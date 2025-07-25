/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "fvMatrix.H"
#include "oversetFvPatchField.H"
#include "calculatedProcessorFvPatchField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "processorFvPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetFvMeshBase::scaleConnection
(
    Field<Type>& coeffs,
    const labelUList& types,
    const scalarList& factor,
    const bool setHoleCellValue,
    const label celli,
    const label facei
) const
{
    const label cType = types[celli];
    const scalar f = factor[celli];

    if (cType == cellCellStencil::INTERPOLATED)
    {
        coeffs[facei] *= 1.0-f;
    }
    else if (cType == cellCellStencil::HOLE)
    {
        // Disconnect hole cell from influence of neighbour
        coeffs[facei] = pTraits<Type>::zero;
    }
    else if (cType == cellCellStencil::SPECIAL)
    {
        if (setHoleCellValue)
        {
            // Behave like hole
            coeffs[facei] = pTraits<Type>::zero;
        }
        else
        {
            // Behave like interpolated
            coeffs[facei] *= 1.0-f;
        }
    }
}

template<class GeoField, class PatchType>
void Foam::oversetFvMeshBase::correctBoundaryConditions
(
    typename GeoField::Boundary& bfld,
    const bool typeOnly
)
{
    // Alternative (C++14)
    //
    // bfld.evaluate_if
    // (
    //     [typeOnly](const auto& pfld) -> bool
    //     {
    //         return (typeOnly == bool(isA<PatchType>(pfld)));
    //     },
    //     UPstream::defaultCommsType
    // );

    const UPstream::commsTypes commsType = UPstream::defaultCommsType;
    const label startOfRequests = UPstream::nRequests();

    for (auto& pfld : bfld)
    {
        if (typeOnly == bool(isA<PatchType>(pfld)))
        {
            pfld.initEvaluate(commsType);
        }
    }

    // Wait for outstanding requests (non-blocking)
    UPstream::waitRequests(startOfRequests);

    for (auto& pfld : bfld)
    {
        if (typeOnly == bool(isA<PatchType>(pfld)))
        {
            pfld.evaluate(commsType);
        }
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::oversetFvMeshBase::normalisation
(
    const fvMatrix<Type>& m
) const
{
    // Determine normalisation. This is normally the original diagonal plus
    // remote contributions. This needs to be stabilised for hole cells
    // which can have a zero diagonal. Assume that if any component has
    // a non-zero diagonal the cell does not need stabilisation.
    tmp<scalarField> tnorm(tmp<scalarField>::New(m.diag()));
    scalarField& norm = tnorm.ref();

    // Add remote coeffs to duplicate behaviour of fvMatrix::addBoundaryDiag
    const FieldField<Field, Type>& internalCoeffs = m.internalCoeffs();
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        forAll(internalCoeffs, patchi)
        {
            const labelUList& fc = mesh_.lduAddr().patchAddr(patchi);
            const Field<Type>& intCoeffs = internalCoeffs[patchi];
            const scalarField cmptCoeffs(intCoeffs.component(cmpt));
            forAll(fc, i)
            {
                norm[fc[i]] += cmptCoeffs[i];
            }
        }
    }

    // Count number of problematic cells
    label nZeroDiag = 0;
    forAll(norm, celli)
    {
        const scalar& n = norm[celli];
        if (magSqr(n) < sqr(SMALL))
        {
            //Pout<< "For field " << m.psi().name()
            //    << " have diagonal " << n << " for cell " << celli
            //    << " at:" << cellCentres()[celli] << endl;
            nZeroDiag++;
        }
    }

    if (debug)
    {
        Pout<< "For field " << m.psi().name() << " have zero diagonals for "
            << returnReduce(nZeroDiag, sumOp<label>()) << " cells" << endl;
    }

    if (returnReduceOr(nZeroDiag))
    {
        // Walk out the norm across hole cells

        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        const cellCellStencilObject& overlap = Stencil::New(mesh_);
        const labelUList& types = overlap.cellTypes();

        // label nHoles = 0;
        scalarField extrapolatedNorm(norm);
        forAll(types, celli)
        {
            if (types[celli] == cellCellStencil::HOLE)
            {
                extrapolatedNorm[celli] = -GREAT;
                // ++nHoles;
            }
        }

        bitSet isFront(mesh_.nFaces());
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label ownType = types[own[facei]];
            label neiType = types[nei[facei]];
            if
            (
                (ownType == cellCellStencil::HOLE)
             != (neiType == cellCellStencil::HOLE)
            )
            {
                isFront.set(facei);
            }
        }
        labelList nbrTypes;
        syncTools::swapBoundaryCellList(mesh_, types, nbrTypes);
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            ++facei
        )
        {
            const label ownType = types[own[facei]];
            const label neiType = nbrTypes[facei-mesh_.nInternalFaces()];
            if
            (
                (ownType == cellCellStencil::HOLE)
             != (neiType == cellCellStencil::HOLE)
            )
            {
                isFront.set(facei);
            }
        }


        while (true)
        {
            scalarField nbrNorm;
            syncTools::swapBoundaryCellList(mesh_, extrapolatedNorm, nbrNorm);

            bitSet newIsFront(mesh_.nFaces());
            scalarField newNorm(extrapolatedNorm);

            label nChanged = 0;
            for (const label facei : isFront)
            {
                if (extrapolatedNorm[own[facei]] == -GREAT)
                {
                    // Average owner cell, add faces to newFront
                    newNorm[own[facei]] = cellAverage
                    (
                        types,
                        nbrTypes,
                        extrapolatedNorm,
                        nbrNorm,
                        own[facei],
                        newIsFront
                    );
                    nChanged++;
                }
                if
                (
                    mesh_.isInternalFace(facei)
                 && extrapolatedNorm[nei[facei]] == -GREAT
                )
                {
                    // Average nei cell, add faces to newFront
                    newNorm[nei[facei]] = cellAverage
                    (
                        types,
                        nbrTypes,
                        extrapolatedNorm,
                        nbrNorm,
                        nei[facei],
                        newIsFront
                    );
                    nChanged++;
                }
            }

            if (!returnReduceOr(nChanged))
            {
                break;
            }

            // Transfer new front
            extrapolatedNorm.transfer(newNorm);
            isFront.transfer(newIsFront);
            syncTools::syncFaceList(mesh_, isFront, maxEqOp<unsigned int>());
        }


        forAll(norm, celli)
        {
            scalar& n = norm[celli];
            if (magSqr(n) < sqr(SMALL))
            {
                //Pout<< "For field " << m.psi().name()
                //    << " for cell " << celli
                //    << " at:" << cellCentres()[celli]
                //    << " have norm " << n
                //    << " have extrapolated norm " << extrapolatedNorm[celli]
                //    << endl;
                // Override the norm
                n = extrapolatedNorm[celli];
            }
        }
    }
    return tnorm;
}


template<class Type>
void Foam::oversetFvMeshBase::addInterpolation
(
    fvMatrix<Type>& m,
    const scalarField& normalisation,
    const bool setHoleCellValue,
    const Type& holeCellValue
) const
{
    const cellCellStencilObject& overlap = Stencil::New(mesh_);
    const List<scalarList>& wghts = overlap.cellInterpolationWeights();
    const labelListList& stencil = overlap.cellStencil();
    const scalarList& factor = overlap.cellInterpolationWeight();
    const labelUList& types = overlap.cellTypes();


    // Force asymmetric matrix (if it wasn't already)
    scalarField& lower = m.lower();
    scalarField& upper = m.upper();
    Field<Type>& source = m.source();
    scalarField& diag = m.diag();


    // Get the addressing. Note that the addressing is now extended with
    // any interpolation faces.
    const lduAddressing& addr = lduAddr();
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const lduInterfacePtrsList& interfaces = allInterfaces_;

    if (!isA<fvMeshPrimitiveLduAddressing>(addr))
    {
        FatalErrorInFunction
            << "Problem : addressing is not fvMeshPrimitiveLduAddressing"
            << exit(FatalError);
    }


    // 1. Adapt lduMatrix for additional faces and new ordering
    upper.setSize(upperAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, upper);
    lower.setSize(lowerAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, lower);


    //const label nOldInterfaces = dynamicMotionSolverFvMesh::interfaces().size();
    const label nOldInterfaces = mesh_.fvMesh::interfaces().size();


    if (interfaces.size() > nOldInterfaces)
    {
        // Extend matrix coefficients
        m.internalCoeffs().setSize(interfaces.size());
        m.boundaryCoeffs().setSize(interfaces.size());

        // 1b. Adapt for additional interfaces
        for
        (
            label patchi = nOldInterfaces;
            patchi < interfaces.size();
            patchi++
        )
        {
            const labelUList& fc = interfaces[patchi].faceCells();
            m.internalCoeffs().set(patchi, new Field<Type>(fc.size(), Zero));
            m.boundaryCoeffs().set(patchi, new Field<Type>(fc.size(), Zero));
        }

        // 1c. Adapt field for additional interfaceFields (note: solver uses
        //     GeometricField::scalarInterfaces() to get hold of interfaces)
        typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

        typename GeoField::Boundary& bfld =
            const_cast<GeoField&>(m.psi()).boundaryFieldRef();

        bfld.setSize(interfaces.size());


        // This gets quite interesting: we do not want to add additional
        // fvPatches (since direct correspondence to polyMesh) so instead
        // add a reference to an existing processor patch
        label addPatchi = 0;
        for (label patchi = 0; patchi < nOldInterfaces; patchi++)
        {
            if (isA<processorFvPatch>(bfld[patchi].patch()))
            {
                addPatchi = patchi;
                break;
            }
        }

        for
        (
            label patchi = nOldInterfaces;
            patchi < interfaces.size();
            patchi++
        )
        {
            bfld.set
            (
                patchi,
                new calculatedProcessorFvPatchField<Type>
                (
                    interfaces[patchi],
                    bfld[addPatchi].patch(),    // dummy processorFvPatch
                    m.psi()
                )
            );
        }
    }


    // 2. Adapt fvMatrix level: faceFluxCorrectionPtr
    // Question: do we need to do this?
    // This seems to be set/used only by the gaussLaplacianScheme and
    // fvMatrix:correction, both of which are outside the linear solver.


    // Clear out existing connections on cells to be interpolated
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: could avoid doing the zeroing of the new faces since these
    //       are set to zero anyway.

    forAll(upperAddr, facei)
    {
        const label l = lowerAddr[facei];
        scaleConnection(upper, types, factor, setHoleCellValue, l, facei);
        const label u = upperAddr[facei];
        scaleConnection(lower, types, factor, setHoleCellValue, u, facei);
        /*
        if (types[upperAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect upper from lower
            label celli = upperAddr[facei];
            lower[facei] *= 1.0-factor[celli];
        }
        if (types[lowerAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect lower from upper
            label celli = lowerAddr[facei];
            upper[facei] *= 1.0-factor[celli];
        }
        */
    }

    for (label patchi = 0; patchi < nOldInterfaces; ++patchi)
    {
        const labelUList& fc = addr.patchAddr(patchi);
        Field<Type>& bCoeffs = m.boundaryCoeffs()[patchi];
        Field<Type>& iCoeffs = m.internalCoeffs()[patchi];
        forAll(fc, i)
        {
            scaleConnection(bCoeffs, types, factor, setHoleCellValue, fc[i], i);

            scaleConnection(iCoeffs, types, factor, setHoleCellValue, fc[i], i);
        }
        /*
        const labelUList& fc = addr.patchAddr(patchi);
        Field<Type>& intCoeffs = m.internalCoeffs()[patchi];
        Field<Type>& bouCoeffs = m.boundaryCoeffs()[patchi];
        forAll(fc, i)
        {
            label celli = fc[i];
            {
                if (types[celli] == cellCellStencil::INTERPOLATED)
                {
                    scalar f = factor[celli];
                    intCoeffs[i] *= 1.0-f;
                    bouCoeffs[i] *= 1.0-f;
                }
                else if (types[celli] == cellCellStencil::HOLE)
                {
                    intCoeffs[i] = pTraits<Type>::zero;
                    bouCoeffs[i] = pTraits<Type>::zero;
                }
            }
        }
        */
    }



    // Modify matrix
    // ~~~~~~~~~~~~~

    // Do hole cells. Note: maybe put into interpolationCells() loop above?
    forAll(types, celli)
    {
        if
        (
            types[celli] == cellCellStencil::HOLE
         || (setHoleCellValue && types[celli] == cellCellStencil::SPECIAL)
        )
        {
            const Type wantedValue
            (
                setHoleCellValue
              ? holeCellValue
              : m.psi()[celli]
            );
            diag[celli] = normalisation[celli];
            source[celli] = normalisation[celli]*wantedValue;
        }
        else if
        (
            types[celli] == cellCellStencil::INTERPOLATED
         || (!setHoleCellValue && types[celli] == cellCellStencil::SPECIAL)
        )
        {
            const scalar f = factor[celli];
            const scalarList& w = wghts[celli];
            const labelList& nbrs = stencil[celli];
            const labelList& nbrFaces = stencilFaces_[celli];
            const labelList& nbrPatches = stencilPatches_[celli];

            diag[celli] *= (1.0-f);
            source[celli] *= (1.0-f);

            forAll(nbrs, nbri)
            {
                const label patchi = nbrPatches[nbri];
                const label facei = nbrFaces[nbri];

                if (patchi == -1)
                {
                    const label nbrCelli = nbrs[nbri];
                    // Add the coefficients
                    const scalar s = normalisation[celli]*f*w[nbri];

                    scalar& u = upper[facei];
                    scalar& l = lower[facei];
                    if (celli < nbrCelli)
                    {
                        diag[celli] += s;
                        u += -s;
                    }
                    else
                    {
                        diag[celli] += s;
                        l += -s;
                    }
                }
                else
                {
                    // Patch face. Store in boundaryCoeffs. Note sign change.
                    //const label globalCelli = globalCellIDs[nbrs[nbri]];
                    //const label proci =
                    //    globalNumbering.whichProcID(globalCelli);
                    //const label remoteCelli =
                    //    globalNumbering.toLocal(proci, globalCelli);
                    //
                    //Pout<< "for cell:" << celli
                    //    << " need weight from remote slot:" << nbrs[nbri]
                    //    << " proc:" << proci << " remote cell:" << remoteCelli
                    //    << " patch:" << patchi
                    //    << " patchFace:" << facei
                    //    << " weight:" << w[nbri]
                    //    << endl;

                    const scalar s = normalisation[celli]*f*w[nbri];
                    m.boundaryCoeffs()[patchi][facei] += pTraits<Type>::one*s;
                    m.internalCoeffs()[patchi][facei] += pTraits<Type>::one*s;

                    // Note: do NOT add to diagonal - this is in the
                    //       internalCoeffs and gets added to the diagonal
                    //       inside fvMatrix::solve
                }
            }
        }
            /*
            label startLabel = ownerStartAddr[celli];
            label endLabel = ownerStartAddr[celli + 1];

            for (label facei = startLabel; facei < endLabel; facei++)
            {
                upper[facei] = 0.0;
            }

            startLabel = addr.losortStartAddr()[celli];
            endLabel = addr.losortStartAddr()[celli + 1];

            for (label i = startLabel; i < endLabel; i++)
            {
                label facei = losortAddr[i];
                lower[facei] = 0.0;
            }

            diag[celli] = normalisation[celli];
            source[celli] = normalisation[celli]*m.psi()[celli];
            */
    }


    //const globalIndex globalNumbering(V().size());
    //labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    //forAll(V(), cellI)
    //{
    //    globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    //}
    //overlap.cellInterpolationMap().distribute(globalCellIDs);

/*
    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];

        const scalar f = factor[celli];
        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const labelList& nbrFaces = stencilFaces_[celli];
        const labelList& nbrPatches = stencilPatches_[celli];

        if (types[celli] == cellCellStencil::HOLE)
        {
            FatalErrorInFunction << "Found HOLE cell " << celli
                << " at:" << C()[celli]
                << " . Should this be in interpolationCells()????"
                << abort(FatalError);
        }
        else
        {
            // Create interpolation stencil

            diag[celli] *= (1.0-f);
            source[celli] *= (1.0-f);

            forAll(nbrs, nbri)
            {
                label patchi = nbrPatches[nbri];
                label facei = nbrFaces[nbri];

                if (patchi == -1)
                {
                    label nbrCelli = nbrs[nbri];

                    // Add the coefficients
                    const scalar s = normalisation[celli]*f*w[nbri];

                    scalar& u = upper[facei];
                    scalar& l = lower[facei];
                    if (celli < nbrCelli)
                    {
                        diag[celli] += s;
                        u += -s;
                    }
                    else
                    {
                        diag[celli] += s;
                        l += -s;
                    }
                }
                else
                {
                    // Patch face. Store in boundaryCoeffs. Note sign change.
                    //const label globalCelli = globalCellIDs[nbrs[nbri]];
                    //const label proci =
                    //    globalNumbering.whichProcID(globalCelli);
                    //const label remoteCelli =
                    //    globalNumbering.toLocal(proci, globalCelli);
                    //
                    //Pout<< "for cell:" << celli
                    //    << " need weight from remote slot:" << nbrs[nbri]
                    //    << " proc:" << proci << " remote cell:" << remoteCelli
                    //    << " patch:" << patchi
                    //    << " patchFace:" << facei
                    //    << " weight:" << w[nbri]
                    //    << endl;

                    const scalar s = normalisation[celli]*f*w[nbri];
                    m.boundaryCoeffs()[patchi][facei] += pTraits<Type>::one*s;
                    m.internalCoeffs()[patchi][facei] += pTraits<Type>::one*s;

                    // Note: do NOT add to diagonal - this is in the
                    //       internalCoeffs and gets added to the diagonal
                    //       inside fvMatrix::solve
                }
            }

            //if (mag(diag[celli]) < SMALL)
            //{
            //    Pout<< "for cell:" << celli
            //        << " at:" << this->C()[celli]
            //        << " diag:" << diag[celli] << endl;
            //
            //    forAll(nbrs, nbri)
            //    {
            //        label patchi = nbrPatches[nbri];
            //        label facei = nbrFaces[nbri];
            //
            //        const label globalCelli = globalCellIDs[nbrs[nbri]];
            //        const label proci =
            //            globalNumbering.whichProcID(globalCelli);
            //        const label remoteCelli =
            //            globalNumbering.toLocal(proci, globalCelli);
            //
            //        Pout<< " need weight from slot:" << nbrs[nbri]
            //            << " proc:" << proci << " remote cell:"
            //            << remoteCelli
            //            << " patch:" << patchi
            //            << " patchFace:" << facei
            //            << " weight:" << w[nbri]
            //            << endl;
            //    }
            //    Pout<< endl;
            //}
        }
    }
*/
}


template<class Type>
Foam::SolverPerformance<Type> Foam::oversetFvMeshBase::solveOverset
(
    fvMatrix<Type>& m,
    const dictionary& dict
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
    // Check we're running with bcs that can handle implicit matrix manipulation
    typename GeoField::Boundary& bpsi =
        const_cast<GeoField&>(m.psi()).boundaryFieldRef();

    bool hasOverset = false;
    forAll(bpsi, patchi)
    {
        if (isA<oversetFvPatchField<Type>>(bpsi[patchi]))
        {
            hasOverset = true;
            break;
        }
    }

    if (!hasOverset)
    {
        if (debug)
        {
            Pout<< "oversetFvMeshBase::solveOverset() :"
                << " bypassing matrix adjustment for field " << m.psi().name()
                << endl;
        }
        //return dynamicMotionSolverFvMesh::solve(m, dict);
        return mesh_.fvMesh::solve(m, dict);
    }

    if (debug)
    {
        Pout<< "oversetFvMeshBase::solveOverset() :"
            << " adjusting matrix for interpolation for field "
            << m.psi().name() << endl;
    }

    // Calculate stabilised diagonal as normalisation for interpolation
    const scalarField norm(normalisation(m));

    if (debug && mesh_.time().writeTime())
    {
        volScalarField scale
        (
            IOobject
            (
                m.psi().name() + "_normalisation",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        scale.ref().field() = norm;
        correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(scale.boundaryFieldRef(), false);
        scale.write();

        if (debug)
        {
            Pout<< "oversetFvMeshBase::solveOverset() :"
                << " writing matrix normalisation for field " << m.psi().name()
                << " to " << scale.name() << endl;
        }
    }


    // Switch to extended addressing (requires mesh::update() having been
    // called)
    active(true);

    // Adapt matrix
    scalarField oldUpper(m.upper());
    scalarField oldLower(m.lower());
    FieldField<Field, Type> oldInt(m.internalCoeffs());
    FieldField<Field, Type> oldBou(m.boundaryCoeffs());

    Field<Type> oldSource(m.source());
    scalarField oldDiag(m.diag());

    const label oldSize = bpsi.size();

    // Insert the interpolation into the matrix (done inside
    // oversetFvPatchField<Type>::manipulateMatrix)
    m.boundaryManipulate(bpsi);

    // Swap psi values so added patches have patchNeighbourField
    correctBoundaryConditions<GeoField, calculatedProcessorFvPatchField<Type>>
    (
        bpsi,
        true
    );

    // Use lower level solver
    //SolverPerformance<Type> s(dynamicMotionSolverFvMesh::solve(m, dict));
    SolverPerformance<Type> s(mesh_.fvMesh::solve(m, dict));

    // Restore boundary field
    bpsi.setSize(oldSize);

    // Restore matrix
    m.upper().transfer(oldUpper);
    m.lower().transfer(oldLower);

    m.source().transfer(oldSource);
    m.diag().transfer(oldDiag);

    m.internalCoeffs().transfer(oldInt);
    m.boundaryCoeffs().transfer(oldBou);


    const cellCellStencilObject& overlap = Stencil::New(mesh_);
    const labelUList& types = overlap.cellTypes();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    auto& psi =
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>(m.psi());

    // Mirror hole cell values next to calculated cells
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const label ownType = types[own[facei]];
        const label neiType = types[nei[facei]];

        if
        (
            ownType == cellCellStencil::HOLE
         && neiType == cellCellStencil::CALCULATED)
        {
            psi[own[facei]] = psi[nei[facei]];
        }
        else if
        (
            ownType == cellCellStencil::CALCULATED
         && neiType == cellCellStencil::HOLE
        )
        {
             psi[nei[facei]] = psi[own[facei]];
        }
    }

    // Switch to original addressing
    active(false);

    return s;
}


template<class Type>
void Foam::oversetFvMeshBase::write
(
    Ostream& os,
    const fvMatrix<Type>& m,
    const lduAddressing& addr,
    const lduInterfacePtrsList& interfaces
) const
{
    os  << "** Matrix **" << endl;
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& ownerStartAddr = addr.ownerStartAddr();
    const labelUList& losortAddr = addr.losortAddr();

    const scalarField& lower = m.lower();
    const scalarField& upper = m.upper();
    const Field<Type>& source = m.source();
    const scalarField& diag = m.diag();


    // Invert patch addressing
    labelListList cellToPatch(addr.size());
    labelListList cellToPatchFace(addr.size());
    {
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                const labelUList& fc = interfaces[patchi].faceCells();

                forAll(fc, i)
                {
                    cellToPatch[fc[i]].append(patchi);
                    cellToPatchFace[fc[i]].append(i);
                }
            }
        }
    }

    forAll(source, celli)
    {
        os  << "cell:" << celli << " diag:" << diag[celli]
            << " source:" << source[celli] << endl;

        label startLabel = ownerStartAddr[celli];
        label endLabel = ownerStartAddr[celli + 1];

        for (label facei = startLabel; facei < endLabel; facei++)
        {
            os  << "    face:" << facei
                << " nbr:" << upperAddr[facei] << " coeff:" << upper[facei]
                << endl;
        }

        startLabel = addr.losortStartAddr()[celli];
        endLabel = addr.losortStartAddr()[celli + 1];

        for (label i = startLabel; i < endLabel; i++)
        {
            label facei = losortAddr[i];
            os  << "    face:" << facei
                << " nbr:" << lowerAddr[facei] << " coeff:" << lower[facei]
                << endl;
        }

        forAll(cellToPatch[celli], i)
        {
            label patchi = cellToPatch[celli][i];
            label patchFacei = cellToPatchFace[celli][i];

            os  << "    patch:" << patchi
                << " patchface:" << patchFacei
                << " intcoeff:" << m.internalCoeffs()[patchi][patchFacei]
                << " boucoeff:" << m.boundaryCoeffs()[patchi][patchFacei]
                << endl;
        }
    }
    forAll(m.internalCoeffs(), patchi)
    {
        if (m.internalCoeffs().set(patchi))
        {
            const labelUList& fc = addr.patchAddr(patchi);

            os  << "patch:" << patchi
                //<< " type:" << interfaces[patchi].type()
                << " size:" << fc.size() << endl;
            if
            (
                interfaces.set(patchi)
             && isA<processorLduInterface>(interfaces[patchi])
            )
            {
                const processorLduInterface& ppp =
                    refCast<const processorLduInterface>(interfaces[patchi]);
                os  << "(processor with my:" << ppp.myProcNo()
                    << " nbr:" << ppp.neighbProcNo()
                    << ")" << endl;
            }

            forAll(fc, i)
            {
                os  << "    cell:" << fc[i]
                    << " int:" << m.internalCoeffs()[patchi][i]
                    << " boun:" << m.boundaryCoeffs()[patchi][i]
                    << endl;
            }
        }
    }


    lduInterfaceFieldPtrsList interfaceFields =
        m.psi().boundaryField().scalarInterfaces();
    forAll(interfaceFields, inti)
    {
        if (interfaceFields.set(inti))
        {
            os  << "interface:" << inti
                << " if type:" << interfaceFields[inti].interface().type()
                << " fld type:" << interfaceFields[inti].type() << endl;
        }
    }

    os  << "** End of Matrix **" << endl;
}


template<class GeoField>
void Foam::oversetFvMeshBase::correctCoupledBoundaryConditions(GeoField& fld)
{
    auto& bfld = fld.boundaryFieldRef();

    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    const label startOfRequests = UPstream::nRequests();

    for (auto& pfld : bfld)
    {
        if (pfld.coupled())
        {
            //Pout<< "initEval of " << pfld.patch().name() << endl;
            pfld.initEvaluate(commsType);
        }
    }

    // Wait for outstanding requests (non-blocking)
    UPstream::waitRequests(startOfRequests);

    for (auto& pfld : bfld)
    {
        if (pfld.coupled())
        {
            //Pout<< "eval of " << pfld.patch().name() << endl;
            pfld.evaluate(commsType);
        }
    }
}


template<class GeoField>
void Foam::oversetFvMeshBase::checkCoupledBC(const GeoField& fld)
{
    Pout<< "** starting checking of " << fld.name() << endl;

    GeoField fldCorr(fld.name()+"_correct", fld);
    correctCoupledBoundaryConditions(fldCorr);

    const typename GeoField::Boundary& bfld = fld.boundaryField();
    const typename GeoField::Boundary& bfldCorr = fldCorr.boundaryField();

    forAll(bfld, patchi)
    {
        const typename GeoField::Patch& pfld = bfld[patchi];
        const typename GeoField::Patch& pfldCorr = bfldCorr[patchi];

        Pout<< "Patch:" << pfld.patch().name() << endl;

        forAll(pfld, i)
        {
            if (pfld[i] != pfldCorr[i])
            {
                Pout<< "    " << i << "  orig:" << pfld[i]
                    << " corrected:" << pfldCorr[i] << endl;
            }
        }
    }
    Pout<< "** end of " << fld.name() << endl;
}


// ************************************************************************* //
