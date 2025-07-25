/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
#include "cellCellStencil.H"
#include "cellCellStencilObject.H"
#include "oversetFvMeshBase.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    setHoleCellValue_(false),
    fluxCorrection_(false),
    interpolateHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::min),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_(),
    fringeFaces_(),
    zoneId_(-1)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    setHoleCellValue_(ptf.setHoleCellValue_),
    fluxCorrection_(ptf.fluxCorrection_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(ptf.fringeUpperCoeffs_),
    fringeLowerCoeffs_(ptf.fringeLowerCoeffs_),
    fringeFaces_(ptf.fringeFaces_),
    zoneId_(ptf.zoneId_)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    oversetPatch_(refCast<const oversetFvPatch>(p, dict)),
    setHoleCellValue_(dict.getOrDefault("setHoleCellValue", false)),
    fluxCorrection_
    (
        dict.getOrDefaultCompat
        (
            "fluxCorrection",
            {{"massCorrection", 2206}},
            false
        )
    ),
    interpolateHoleCellValue_
    (
        dict.getOrDefault("interpolateHoleCellValue", false)
    ),
    holeCellValue_
    (
        setHoleCellValue_
      ? dict.get<Type>("holeCellValue")
      : pTraits<Type>::min
    ),
    fringeUpperCoeffs_(),
    fringeLowerCoeffs_(),
    fringeFaces_(),
    zoneId_(dict.getOrDefault<label>("zone", -1))
{
    // Use 'value' supplied, or set to internal field
    if (!this->readValueEntry(dict))
    {
        fvPatchField<Type>::extrapolateInternal();
    }
}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    oversetPatch_(ptf.oversetPatch_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    fluxCorrection_(ptf.fluxCorrection_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(ptf.fringeUpperCoeffs_),
    fringeLowerCoeffs_(ptf.fringeLowerCoeffs_),
    fringeFaces_(ptf.fringeFaces_),
    zoneId_(ptf.zoneId_)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    oversetPatch_(ptf.oversetPatch_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    fluxCorrection_(ptf.fluxCorrection_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_),
    fringeUpperCoeffs_(ptf.fringeUpperCoeffs_),
    fringeLowerCoeffs_(ptf.fringeLowerCoeffs_),
    fringeFaces_(ptf.fringeFaces_),
    zoneId_(ptf.zoneId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetFvPatchField<Type>::storeFringeCoefficients
(
    const fvMatrix<Type>& matrix
)
{
    const fvMesh& mesh = this->internalField().mesh();

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelList& zoneID = overlap.zoneID();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    label fringesFaces = 0;

    forAll(own, facei)
    {
        const label zonei = zoneID[own[facei]];

        const label ownType = cellTypes[own[facei]];
        const label neiType = cellTypes[nei[facei]];

        const bool ownCalc =
            (ownType == cellCellStencil::CALCULATED)
         && (neiType == cellCellStencil::INTERPOLATED);

        const bool neiCalc =
            (ownType == cellCellStencil::INTERPOLATED)
         && (neiType == cellCellStencil::CALCULATED);

        const bool ownNei = (ownCalc || neiCalc);

        if
        (
            (ownNei && (zonei == zoneId_)) || (ownNei && (zoneId_ == -1))
        )
        {
            fringesFaces++;
        }
    }

    const fvPatchList& patches = mesh.boundary();

    labelList neiCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, neiCellTypes);
    {
        forAll(patches, patchi)
        {
            const fvPatch& curPatch = patches[patchi];

            const labelUList& fc = curPatch.faceCells();

            const label start = curPatch.start();

            forAll(fc, i)
            {
                const label facei = start + i;
                const label celli = fc[i];
                const label ownType = cellTypes[celli];
                const label neiType = neiCellTypes[facei-mesh.nInternalFaces()];

                const label zonei = zoneID[celli];

                const bool ownCalc =
                    (ownType == cellCellStencil::CALCULATED)
                 && (neiType == cellCellStencil::INTERPOLATED);


                if (ownCalc && (zonei == zoneId_))
                {
                    fringesFaces++;
                }
            }
        }
    }
    fringeUpperCoeffs_.setSize(fringesFaces, Zero);
    fringeLowerCoeffs_.setSize(fringesFaces, Zero);
    fringeFaces_.setSize(fringesFaces, -1);

    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    fringesFaces = 0;
    forAll(own, facei)
    {
        const label zonei = zoneID[own[facei]];

        const label ownType = cellTypes[own[facei]];
        const label neiType = cellTypes[nei[facei]];

        const bool ownCalc =
            (ownType == cellCellStencil::CALCULATED)
         && (neiType == cellCellStencil::INTERPOLATED);

        const bool neiCalc =
            (ownType == cellCellStencil::INTERPOLATED)
         && (neiType == cellCellStencil::CALCULATED);

        const bool ownNei = (ownCalc || neiCalc);

        if
        (
            (ownNei && (zonei == zoneId_)) || (ownNei && (zoneId_ == -1))
        )
        {
            fringeUpperCoeffs_[fringesFaces] = upper[facei];
            fringeLowerCoeffs_[fringesFaces] = lower[facei];
            fringeFaces_[fringesFaces] = facei;
            fringesFaces++;
        }
    }

    forAll(boundaryMesh, patchi)
    {
        const polyPatch& p = boundaryMesh[patchi];

        if (isA<coupledPolyPatch>(p))
        {
            const labelUList& fc = p.faceCells();
            const label start = p.start();

            forAll(fc, i)
            {
                const label facei = start + i;
                const label celli = fc[i];
                const label ownType = cellTypes[celli];
                const label neiType = neiCellTypes[facei-mesh.nInternalFaces()];

                const label zonei = zoneID[celli];

                const bool ownCalc =
                    (ownType == cellCellStencil::CALCULATED)
                 && (neiType == cellCellStencil::INTERPOLATED);

                const bool neiCalc =
                    (neiType == cellCellStencil::CALCULATED)
                 && (ownType == cellCellStencil::INTERPOLATED);

                if ((ownCalc||neiCalc)  && (zonei == zoneId_))
                {
                    fringeLowerCoeffs_[fringesFaces] =
                        component
                        (
                            matrix.internalCoeffs()[patchi][facei],
                            0
                        );

                    fringeUpperCoeffs_[fringesFaces] =
                        component
                        (
                            matrix.boundaryCoeffs()[patchi][facei],
                            0
                        );

                    fringeFaces_[fringesFaces] = facei;

                    fringesFaces++;
                }
            }
        }
    }
}


template<class Type>
void Foam::oversetFvPatchField<Type>::fringeFlux
(
    const fvMatrix<Type>& matrix,
    const surfaceScalarField& phi
) const
{
    scalar massIn = 0;
    scalar phiIn = 0;

    const Field<Type>& psi = matrix.psi();

    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    if (this->oversetPatch_.master())
    {
        const fvMesh& mesh = this->internalField().mesh();
        const cellCellStencilObject& overlap = Stencil::New(mesh);
        const labelList& cellTypes = overlap.cellTypes();
        const labelList& zoneID = overlap.zoneID();

        // Check all faces on the outside of interpolated cells
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        label fringesFaces = 0;
        forAll(own, facei)
        {
            const label zonei = zoneID[own[facei]];

            const label ownType = cellTypes[own[facei]];
            const label neiType = cellTypes[nei[facei]];

            const bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            const bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);

            const bool ownNei = (ownCalc || neiCalc);

            if
            (
                (ownNei && (zonei == zoneId_)) || (ownNei && (zoneId_ == -1))
            )
            {
                const label fringei = fringeFaces_[fringesFaces];

                // Get fringe upper/lower coeffs
                //const scalar& ufc = fringeUpperCoeffs_[fringesFaces];
                //const scalar& lfc = fringeLowerCoeffs_[fringesFaces];

                const scalar& ufc = upper[fringei];
                const scalar& lfc = lower[fringei];

                const Type curFlux =
                    ufc*psi[nei[fringei]] - lfc*psi[own[fringei]];

                if (neiCalc)
                {
                    phiIn -= phi[fringei];
                    massIn -= curFlux;
                }
                else
                {
                    phiIn += phi[fringei];
                    massIn += curFlux;
                }
                fringesFaces++;
            }
        }
    }

    reduce(massIn, sumOp<scalar>());
    reduce(phiIn, sumOp<scalar>());

    Info << " gSum(phi) on fringes " << phiIn << endl;
    Info << " gSum(p.flux) on fringes " << massIn << endl;
}


template<class Type>
void Foam::oversetFvPatchField<Type>::adjustPsi
(
    solveScalarField& psi,
    const lduAddressing& lduAddr,
    solveScalarField& result
) const
{
    const fvMesh& mesh = this->internalField().mesh();

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelList& zoneID = overlap.zoneID();

    // Pass-1: accumulate all fluxes, calculate correction factor

    scalarField interpolatedCoeffs(fringeUpperCoeffs_.size(), Zero);

    // Options for scaling corrections
    scalar massIn = 0;
    scalar offDiagCoeffs = 0;
    labelField facePerCell(cellTypes.size(), 0);

    // Check all faces on the outside of interpolated cells
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    label fringesFaces = 0;
    {
        forAll(own, facei)
        {
            const label zonei = zoneID[own[facei]];

            const label ownType = cellTypes[own[facei]];
            const label neiType = cellTypes[nei[facei]];

            const bool ownCalc =
                (ownType == cellCellStencil::CALCULATED)
             && (neiType == cellCellStencil::INTERPOLATED);

            const bool neiCalc =
                (ownType == cellCellStencil::INTERPOLATED)
             && (neiType == cellCellStencil::CALCULATED);

            const bool ownNei = (ownCalc || neiCalc);

            if
            (
                (ownNei && (zonei == zoneId_)) || (ownNei && (zoneId_ == -1))
            )
            {
                // Get fringe upper/lower coeffs
                const scalar& ufc = fringeUpperCoeffs_[fringesFaces];
                const scalar& lfc = fringeLowerCoeffs_[fringesFaces];

                const scalar curFlux =
                    ufc*psi[nei[facei]] - lfc*psi[own[facei]];

                if (neiCalc) // interpolated is owner
                {
                    massIn -= curFlux;
                    offDiagCoeffs += lfc;
                    facePerCell[own[facei]]++;
                }
                else
                {
                    massIn += curFlux;
                    offDiagCoeffs += ufc;
                    facePerCell[nei[facei]]++;
                }
                fringesFaces++;
            }
        }
    }

    scalarField weights(facePerCell.size(), scalar(1));
    forAll (weights, celli)
    {
        if (facePerCell[celli] > 1)
        {
            weights[celli] = scalar(1)/facePerCell[celli];
        }
    }

    // Check all coupled faces on the outside of interpolated cells
    labelList neiCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, neiCellTypes);

    const fvPatchList& boundaryMesh = mesh.boundary();

    forAll(boundaryMesh, patchi)
    {
        const polyPatch& p = mesh.boundaryMesh()[patchi];

        if (isA<coupledPolyPatch>(p))
        {
            const auto& coupledPatch = refCast<const coupledPolyPatch>(p);

            const labelUList& fc = p.faceCells();
            const label start = p.start();

            forAll(fc, i)
            {
                const label facei = start + i;
                const label celli = fc[i];
                const label ownType = cellTypes[celli];
                const label neiType = neiCellTypes[facei-mesh.nInternalFaces()];

                const label zonei = zoneID[celli];

                const bool ownCalc =
                    (ownType == cellCellStencil::CALCULATED)
                 && (neiType == cellCellStencil::INTERPOLATED);

                const bool neiCalc =
                    (neiType == cellCellStencil::CALCULATED)
                 && (ownType == cellCellStencil::INTERPOLATED);

                if ((ownCalc||neiCalc)  && (zonei == zoneId_))
                {
                    const scalar psiOwn = psi[celli];
                    const scalar& lfc = fringeLowerCoeffs_[fringesFaces];
                    const scalar curFlux = lfc*psiOwn;

                    if (ownCalc)
                    {
                        massIn -= curFlux;

                        if (coupledPatch.owner())
                        {
                            offDiagCoeffs -= lfc;
                        }

                        fringesFaces++;
                    }
                    else
                    {
                        massIn += curFlux;

                        if (coupledPatch.owner())
                        {
                            offDiagCoeffs -= lfc;
                        }

                        fringesFaces++;
                    }
                }
            }
        }
    }

    reduce(massIn, sumOp<scalar>());
    reduce(offDiagCoeffs, sumOp<scalar>());

    scalar psiCorr = -massIn/offDiagCoeffs;

    forAll (cellTypes, celli)
    {
        const bool bInter = (cellTypes[celli] == cellCellStencil::INTERPOLATED);
        const label zonei = zoneID[celli];
        if
        (
            (bInter && (zonei == zoneId_)) ||(bInter && (zoneId_ == -1))
        )
        {
            psi[celli] += psiCorr;
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::oversetFvPatchField<Type>::
patchNeighbourField() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
void Foam::oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (this->oversetPatch_.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing. Flux correction already done
            // in the linear solver (through the initUpdateInterfaceMatrix
            // method below)
            if (debug)
            {
                Info<< "Skipping overset interpolation for solved-for field "
                    << fldName << endl;
            }
        }
        else if (!fvSchemes.found("oversetInterpolation"))
        {
            IOWarningInFunction(fvSchemes)
                << "Missing required dictionary entry"
                << " 'oversetInterpolation'"
                << ". Skipping overset interpolation for field "
                << fldName << endl;
        }
        else if (fvSchemes.found("oversetInterpolationRequired"))
        {
            // Backwards compatibility mode: only interpolate what is
            // explicitly mentioned

            if (fvSchemes.found("oversetInterpolationSuppressed"))
            {
                FatalIOErrorInFunction(fvSchemes)
                    << "Cannot have both dictionary entry"
                    << " 'oversetInterpolationSuppresed' and "
                    << " 'oversetInterpolationRequired' for field "
                    << fldName << exit(FatalIOError);
            }
            const dictionary& intDict = fvSchemes.subDict
            (
                "oversetInterpolationRequired"
            );
            if (intDict.found(fldName))
            {
                if (debug)
                {
                    Info<< "Interpolating field " << fldName << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.
                mesh.interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
            }
            else if (debug)
            {
                Info<< "Skipping overset interpolation for field "
                    << fldName << endl;
            }
        }
        else
        {
            const dictionary* dictPtr
            (
                fvSchemes.findDict("oversetInterpolationSuppressed")
            );

            const wordHashSet& suppress =
                Stencil::New(mesh).nonInterpolatedFields();

            bool skipInterpolate = suppress.found(fldName);

            if (dictPtr)
            {
                skipInterpolate =
                    skipInterpolate
                 || dictPtr->found(fldName);
            }

            if (skipInterpolate)
            {
                if (debug)
                {
                    Info<< "Skipping suppressed overset interpolation"
                        << " for field " << fldName << endl;
                }
            }
            else
            {
                if (debug)
                {
                    Info<< "Interpolating non-suppressed field " << fldName
                        << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.

                const cellCellStencil& overlap = Stencil::New(mesh);

                Field<Type>& fld =
                    const_cast<Field<Type>&>(this->primitiveField());

                // tbd: different weights for different variables ...
                cellCellStencil::interpolate
                (
                    fld,
                    mesh,
                    overlap,
                    overlap.cellInterpolationWeights()
                );

                if (this->setHoleCellValue_)
                {
                    const labelUList& types = overlap.cellTypes();
                    label nConstrained = 0;
                    forAll(types, celli)
                    {
                        const label cType = types[celli];
                        if
                        (
                            cType == cellCellStencil::HOLE
                         || cType == cellCellStencil::SPECIAL
                        )
                        {
                            fld[celli] = this->holeCellValue_;
                            nConstrained++;
                        }
                    }

                    if (debug)
                    {
                        Pout<< FUNCTION_NAME << " field:" << fldName
                            << " patch:" << this->oversetPatch_.name()
                            << " set:" << nConstrained << " cells to:"
                            << this->holeCellValue_ << endl;
                    }
                }
            }
        }
    }

    coupledFvPatchField<Type>::initEvaluate(commsType);
}


template<class Type>
void Foam::oversetFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    if (this->manipulatedMatrix())
    {
        return;
    }

    const oversetFvPatch& ovp = this->oversetPatch_;

    if (ovp.master())
    {
        if (fluxCorrection_ || (debug & 2))
        {
            storeFringeCoefficients(matrix);
        }

        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const word& fldName = this->internalField().name();

        // Try to find out if the solve routine comes from the unadapted mesh
        // TBD. This should be cleaner.
        if (&mesh.lduAddr() == &mesh.fvMesh::lduAddr())
        {
            return;
        }

        if (debug)
        {
            Pout<< FUNCTION_NAME << " field:" << fldName
                << " patch:" << ovp.name() << endl;
        }


        // Calculate stabilised diagonal as normalisation for interpolation
        const scalarField norm
        (
            dynamic_cast<const oversetFvMeshBase&>(mesh).normalisation(matrix)
        );

        const cellCellStencil& overlap = Stencil::New(mesh);
        const labelUList& types = overlap.cellTypes();
        const labelListList& stencil = overlap.cellStencil();

        dynamic_cast<const oversetFvMeshBase&>(mesh).addInterpolation
        (
            matrix,
            norm,
            this->setHoleCellValue_,
            this->holeCellValue_
        );

        if (debug)
        {
            pointField allCs(mesh.cellCentres());
            const mapDistribute& map = overlap.cellInterpolationMap();
            map.distribute(allCs, false, UPstream::msgType()+1);

            // Make sure we don't interpolate from a hole

            scalarField marker(this->primitiveField().size(), 0);
            forAll(types, celli)
            {
                if (types[celli] == cellCellStencil::HOLE)
                {
                    marker[celli] = 1.0;
                }
            }
            cellCellStencil::interpolate
            (
                marker,
                mesh,
                overlap,
                overlap.cellInterpolationWeights()
            );

            forAll(marker, celli)
            {
                if
                (
                    types[celli] == cellCellStencil::INTERPOLATED
                 && marker[celli] > SMALL
                )
                {
                    WarningInFunction
                        << " field:" << fldName
                        << " patch:" << ovp.name()
                        << " found:" << celli
                        << " at:" << mesh.cellCentres()[celli]
                        << " donorSlots:" << stencil[celli]
                        << " at:"
                        << UIndirectList<point>(allCs, stencil[celli])
                        << " amount-of-hole:" << marker[celli]
                        << endl;
                }
            }

            // Make sure we don't have matrix coefficients for interpolated
            // or hole cells

            const lduAddressing& addr = mesh.lduAddr();
            const labelUList& upperAddr = addr.upperAddr();
            const labelUList& lowerAddr = addr.lowerAddr();
            const scalarField& lower = matrix.lower();
            const scalarField& upper = matrix.upper();

            forAll(lowerAddr, facei)
            {
                const label l = lowerAddr[facei];
                const bool lHole = (types[l] == cellCellStencil::HOLE);
                const label u = upperAddr[facei];
                const bool uHole = (types[u] == cellCellStencil::HOLE);

                if
                (
                    (lHole && upper[facei] != 0.0)
                 || (uHole && lower[facei] != 0.0)
                )
                {
                    FatalErrorInFunction
                        << "Hole-neighbouring face:" << facei
                        << " lower:" << l
                        << " type:" << types[l]
                        << " coeff:" << lower[facei]
                        << " upper:" << upperAddr[facei]
                        << " type:" << types[u]
                        << " coeff:" << upper[facei]
                        << exit(FatalError);
                }


                // Empty donor list: treat like hole but still allow to
                // influence neighbouring domains
                const bool lEmpty =
                (
                    types[l] == cellCellStencil::INTERPOLATED
                 && stencil[l].empty()
                );
                const bool uEmpty =
                (
                    types[u] == cellCellStencil::INTERPOLATED
                 && stencil[u].empty()
                );

                if
                (
                    (lEmpty && upper[facei] != 0.0)
                 || (uEmpty && lower[facei] != 0.0)
                )
                {
                    FatalErrorInFunction
                        << "Still connected face:" << facei << " lower:" << l
                        << " type:" << types[l]
                        << " coeff:" << lower[facei]
                        << " upper:" << u
                        << " type:" << types[u]
                        << " coeff:" << upper[facei]
                        << exit(FatalError);
                }
            }

            forAll(matrix.internalCoeffs(), patchi)
            {
                const labelUList& fc = addr.patchAddr(patchi);
                const Field<Type>& bouCoeffs = matrix.boundaryCoeffs()[patchi];

                forAll(fc, i)
                {
                    const label celli = fc[i];
                    const bool lHole = (types[celli] == cellCellStencil::HOLE);

                    if (lHole && bouCoeffs[i] != pTraits<Type>::zero)
                    {
                        FatalErrorInFunction
                            << "Patch:" << patchi
                            << " patchFace:" << i
                            << " lower:" << celli
                            << " type:" << types[celli]
                            << " bouCoeff:" << bouCoeffs[i]
                            << exit(FatalError);
                    }

                    // Check whether this is influenced by neighbouring domains
                    const bool lEmpty =
                    (
                        types[celli] == cellCellStencil::INTERPOLATED
                     && stencil[celli].empty()
                    );

                    if (lEmpty && bouCoeffs[i] != pTraits<Type>::zero)
                    {
                        FatalErrorInFunction
                            << "Patch:" << patchi
                            << " patchFace:" << i
                            << " lower:" << celli
                            << " type:" << types[celli]
                            << " bouCoeff:" << bouCoeffs[i]
                            << exit(FatalError);
                    }
                }
            }


            // Make sure that diagonal is non-zero. Note: should add
            // boundaryCoeff ...
            const FieldField<Field, Type>& internalCoeffs =
                matrix.internalCoeffs();

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                // Replacement for m.addBoundaryDiag(norm, cmpt);
                scalarField diag(matrix.diag());
                forAll(internalCoeffs, patchi)
                {
                    const labelUList& fc = addr.patchAddr(patchi);
                    const Field<Type>& intCoeffs = internalCoeffs[patchi];
                    const scalarField cmptCoeffs(intCoeffs.component(cmpt));
                    forAll(fc, i)
                    {
                        diag[fc[i]] += cmptCoeffs[i];
                    }
                }

                forAll(diag, celli)
                {
                    if (mag(diag[celli]) < SMALL)
                    {
                        FatalErrorInFunction
                            << "Patch:" << ovp.name()
                            << " cell:" << celli
                            << " at:" << mesh.cellCentres()[celli]
                            << " diag:" << diag[celli]
                            << exit(FatalError);
                    }
                }
            }
        }
    }

     fvPatchField<Type>::manipulateMatrix(matrix);
}


template<class Type>
void Foam::oversetFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label interfacei,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
}


template<class Type>
void Foam::oversetFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    auto& psi = const_cast<solveScalarField&>(psiInternal);

    if (fluxCorrection_ && this->oversetPatch_.master())
    {
        adjustPsi(psi, lduAddr, result);
    }
}


template<class Type>
void Foam::oversetFvPatchField<Type>::write(Ostream& os) const
{
    coupledFvPatchField<Type>::write(os);

    if (this->setHoleCellValue_)
    {
        os.writeEntry("setHoleCellValue", setHoleCellValue_);
        os.writeEntry("holeCellValue", holeCellValue_);
        os.writeEntryIfDifferent
        (
            "interpolateHoleCellValue",
            false,
            interpolateHoleCellValue_
        );
    }
    os.writeEntryIfDifferent
    (
        "fluxCorrection",
        false,
        fluxCorrection_
    );
    os.writeEntryIfDifferent
    (
        "zone",
        label(-1),
        zoneId_
    );
}


// ************************************************************************* //
