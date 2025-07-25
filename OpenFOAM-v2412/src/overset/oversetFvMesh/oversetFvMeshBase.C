/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022,2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "oversetFvMeshBase.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"
#include "zeroGradientFvPatchFields.H"
#include "lduPrimitiveProcessorInterface.H"
#include "globalIndex.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetFvMeshBase, 0);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

bool Foam::oversetFvMeshBase::updateAddressing() const
{
    const cellCellStencilObject& overlap = Stencil::New(mesh_);

    // The (processor-local part of the) stencil determines the local
    // faces to add to the matrix. tbd: parallel
    const labelListList& stencil = overlap.cellStencil();

    // Get the base addressing
    //const lduAddressing& baseAddr = dynamicMotionSolverFvMesh::lduAddr();
    const lduAddressing& baseAddr = mesh_.fvMesh::lduAddr();

    // Add to the base addressing
    labelList lowerAddr;
    labelList upperAddr;
    label nExtraFaces;

    const globalIndex globalNumbering(baseAddr.size());
    labelListList localFaceCells;
    labelListList remoteFaceCells;

    labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    forAll(baseAddr, cellI)
    {
        globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    }
    overlap.cellInterpolationMap().distribute(globalCellIDs);


    reverseFaceMap_ = fvMeshPrimitiveLduAddressing::addAddressing
    (
        baseAddr,
        stencil,
        nExtraFaces,
        lowerAddr,
        upperAddr,
        stencilFaces_,
        globalNumbering,
        globalCellIDs,
        localFaceCells,
        remoteFaceCells
    );

    if (debug)
    {
        Pout<< "oversetFvMeshBase::update() : extended addressing from"
            << " nFaces:" << baseAddr.lowerAddr().size()
            << " to nFaces:" << lowerAddr.size()
            << " nExtraFaces:" << nExtraFaces << endl;
    }

    // Now for the tricky bits. We want to hand out processor faces according
    // to the localFaceCells/remoteFaceCells. Ultimately we need
    // per entry in stencil:
    // - the patch (or -1 for internal faces)
    // - the face (is either an internal face index or a patch face index)

    stencilPatches_.setSize(stencilFaces_.size());

    // Per processor to owner (local)/neighbour (remote)
    List<DynamicList<label>> procOwner(Pstream::nProcs());
    List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());
    forAll(stencil, celli)
    {
        const labelList& nbrs = stencil[celli];
        stencilPatches_[celli].setSize(nbrs.size());
        stencilPatches_[celli] = -1;

        forAll(nbrs, nbri)
        {
            if (stencilFaces_[celli][nbri] == -1)
            {
                const label nbrCelli = nbrs[nbri];
                label globalNbr = globalCellIDs[nbrCelli];
                label proci = globalNumbering.whichProcID(globalNbr);
                label remoteCelli = globalNumbering.toLocal(proci, globalNbr);

                // Overwrite the face to be a patch face
                stencilFaces_[celli][nbri] = procOwner[proci].size();
                stencilPatches_[celli][nbri] = proci;
                procOwner[proci].append(celli);
                dynProcNeighbour[proci].append(remoteCelli);

                //Pout<< "From neighbour proc:" << proci
                //    << " allocating patchFace:" << stencilFaces_[celli][nbri]
                //    << " to get remote cell " << remoteCelli
                //    << " onto local cell " << celli << endl;
            }
        }
    }
    labelListList procNeighbour(dynProcNeighbour.size());
    forAll(procNeighbour, i)
    {
        procNeighbour[i] = std::move(dynProcNeighbour[i]);
    }
    labelListList mySendCells;
    Pstream::exchange<labelList, label>(procNeighbour, mySendCells);

    label nbri = 0;
    forAll(procOwner, proci)
    {
        if (procOwner[proci].size())
        {
            nbri++;
        }
        if (mySendCells[proci].size())
        {
            nbri++;
        }
    }
    remoteStencilInterfaces_.setSize(nbri);
    nbri = 0;

    // E.g. if proc1 needs some data from proc2 and proc2 needs some data from
    //      proc1. We first want the pair : proc1 receive and proc2 send
    //      and then the pair : proc1 send, proc2 receive

    labelList procToInterface(Pstream::nProcs(), -1);

    forAll(procOwner, proci)
    {
        if (proci < Pstream::myProcNo() && procOwner[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci].size()
                    << " from " << proci << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
        else if (proci > Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci].size()
                    << " to " << proci << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
    }
    forAll(procOwner, proci)
    {
        if (proci > Pstream::myProcNo() && procOwner[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci].size()
                    << " from " << proci << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
        else if (proci < Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (debug)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci].size()
                    << " to " << proci << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
    }


    // Rewrite stencilPatches now we know the actual interface (procToInterface)
    for (auto& patches : stencilPatches_)
    {
        for (auto& interface : patches)
        {
            if (interface != -1)
            {
                interface = procToInterface[interface]+mesh_.boundary().size();
            }
        }
    }

    // Get addressing and interfaces of all interfaces

    UPtrList<const labelUList> patchAddr;
    {
        const fvBoundaryMesh& fvp = mesh_.boundary();

        patchAddr.setSize(fvp.size() + remoteStencilInterfaces_.size());

        //allInterfaces_ = dynamicMotionSolverFvMesh::interfaces();
        allInterfaces_ = mesh_.fvMesh::interfaces();
        allInterfaces_.setSize(patchAddr.size());

        forAll(fvp, patchi)
        {
            patchAddr.set(patchi, &fvp[patchi].faceCells());
        }

        forAll(remoteStencilInterfaces_, i)
        {
            label patchi = fvp.size()+i;
            const lduPrimitiveProcessorInterface& pp =
                remoteStencilInterfaces_[i];

            //Pout<< "at patch:" << patchi
            //    << " have procPatch:" << pp.type()
            //    << " from:" << pp.myProcNo()
            //    << " to:" << pp.neighbProcNo()
            //    << " with fc:" << pp.faceCells().size() << endl;

            patchAddr.set(patchi, &pp.faceCells());
            allInterfaces_.set(patchi, &pp);
        }
    }

    const lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            allInterfaces_
        )
    );

    lduPtr_.reset
    (
        new fvMeshPrimitiveLduAddressing
        (
            mesh_.nCells(),
            std::move(lowerAddr),
            std::move(upperAddr),
            patchAddr,
            ps
        )
   );


    // Check
    if (debug)
    {
        const lduAddressing& addr = lduPtr_();  //this->lduAddr();

        Pout<< "Adapted addressing:"
            << " lower:" << addr.lowerAddr().size()
            << " upper:" << addr.upperAddr().size() << endl;

        // Using lduAddressing::patch
        forAll(patchAddr, patchi)
        {
            Pout<< "    " << patchi << "\tpatchAddr:"
                << addr.patchAddr(patchi).size()
                << endl;
        }

        // Using interfaces
        const lduInterfacePtrsList& iFaces = mesh_.interfaces();
        Pout<< "Adapted interFaces:" << iFaces.size() << endl;
        forAll(iFaces, patchi)
        {
            if (iFaces.set(patchi))
            {
                Pout<< "    " << patchi << "\tinterface:"
                    << iFaces[patchi].type() << endl;
            }
        }
    }

    return true;
}


Foam::scalar Foam::oversetFvMeshBase::cellAverage
(
    const labelList& types,
    const labelList& nbrTypes,
    const scalarField& norm,
    const scalarField& nbrNorm,
    const label celli,
    bitSet& isFront
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const cell& cFaces = mesh_.cells()[celli];

    scalar avg = 0.0;
    label nTotal = 0;
    for (const label facei : cFaces)
    {
        if (mesh_.isInternalFace(facei))
        {
            label nbrCelli = (own[facei] == celli ? nei[facei] : own[facei]);
            if (norm[nbrCelli] == -GREAT)
            {
                // Invalid neighbour. Add to front
                isFront.set(facei);
            }
            else
            {
                // Valid neighbour. Add to average
                avg += norm[nbrCelli];
                ++nTotal;
            }
        }
        else
        {
            if (nbrNorm[facei-mesh_.nInternalFaces()] == -GREAT)
            {
                isFront.set(facei);
            }
            else
            {
                avg += nbrNorm[facei-mesh_.nInternalFaces()];
                ++nTotal;
            }
        }
    }

    if (nTotal)
    {
        return avg/nTotal;
    }
    else
    {
        return norm[celli];
    }
}


void Foam::oversetFvMeshBase::writeAgglomeration
(
    const GAMGAgglomeration& agglom
) const
{
    labelList cellToCoarse(identity(mesh_.nCells()));
    labelListList coarseToCell(invertOneToMany(mesh_.nCells(), cellToCoarse));

    // Write initial agglomeration
    {
        volScalarField scalarAgglomeration
        (
            IOobject
            (
                "agglomeration",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        scalarField& fld = scalarAgglomeration.primitiveFieldRef();
        forAll(fld, celli)
        {
            fld[celli] = cellToCoarse[celli];
        }
        fld /= max(fld);
        correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >(scalarAgglomeration.boundaryFieldRef(), false);
        scalarAgglomeration.write();

        Info<< "Writing initial cell distribution to "
            << mesh_.time().timeName() << endl;
    }


    for (label level = 0; level < agglom.size(); level++)
    {
        const labelList& addr = agglom.restrictAddressing(level);
        label coarseSize = max(addr)+1;

        Info<< "Level : " << level << endl
            << "    current size      : "
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;

        forAll(addr, fineI)
        {
            const labelList& cellLabels = coarseToCell[fineI];
            forAll(cellLabels, i)
            {
                cellToCoarse[cellLabels[i]] = addr[fineI];
            }
        }
        coarseToCell = invertOneToMany(coarseSize, cellToCoarse);

        // Write agglomeration
        {
            volScalarField scalarAgglomeration
            (
                IOobject
                (
                    "agglomeration_" + Foam::name(level),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            );
            scalarField& fld = scalarAgglomeration.primitiveFieldRef();
            forAll(fld, celli)
            {
                fld[celli] = cellToCoarse[celli];
            }
            //if (normalise)
            //{
            //    fld /= max(fld);
            //}
            correctBoundaryConditions
            <
                volScalarField,
                oversetFvPatchField<scalar>
            >(scalarAgglomeration.boundaryFieldRef(), false);
            scalarAgglomeration.write();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetFvMeshBase::oversetFvMeshBase(const fvMesh& mesh, bool doInit)
:
    mesh_(mesh),
    active_(false)
{
    // Load stencil (but do not update)
    (void)Stencil::New(mesh_, false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetFvMeshBase::~oversetFvMeshBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduAddressing& Foam::oversetFvMeshBase::lduAddr() const
{
    if (!active_)
    {
        //return dynamicMotionSolverFvMesh::lduAddr();
        return mesh_.fvMesh::lduAddr();
    }
    if (!lduPtr_)
    {
        // Build extended addressing
        updateAddressing();
    }
    return *lduPtr_;
}


Foam::lduInterfacePtrsList Foam::oversetFvMeshBase::interfaces() const
{
    if (!active_)
    {
        return mesh_.fvMesh::interfaces();
    }
    if (!lduPtr_)
    {
        // Build extended addressing
        updateAddressing();
    }
    return allInterfaces_;
}


const Foam::fvMeshPrimitiveLduAddressing&
Foam::oversetFvMeshBase::primitiveLduAddr() const
{
    if (!lduPtr_)
    {
        FatalErrorInFunction
            << "Extended addressing not allocated" << abort(FatalError);
    }

    return *lduPtr_;
}


void Foam::oversetFvMeshBase::clearOut()
{
    // Cell-cell stencil is already mesh object. Clear out local
    // addressing to force rebuilding addressing
    lduPtr_.clear();
}


bool Foam::oversetFvMeshBase::update()
{
    {
        // Calculate the local extra faces for the interpolation. Note: could
        // let demand-driven lduAddr() trigger it but just to make sure.
        updateAddressing();

        // Addressing and/or weights have changed. Make interpolated cells
        // up to date with donors
        interpolateFields();

        return true;
    }

    return false;
}


bool Foam::oversetFvMeshBase::interpolateFields()
{
    const cellCellStencilObject& overlap = Stencil::New(mesh_);

    // Add the stencil suppression list
    wordHashSet suppressed(overlap.nonInterpolatedFields());

    // Use whatever the solver has set up as suppression list
    const dictionary* dictPtr
    (
        mesh_.schemesDict().findDict("oversetInterpolationSuppressed")
    );
    if (dictPtr)
    {
        suppressed.insert(dictPtr->toc());
    }

    overlap.interpolate<volScalarField>(mesh_, suppressed);
    overlap.interpolate<volVectorField>(mesh_, suppressed);
    overlap.interpolate<volSphericalTensorField>(mesh_, suppressed);
    overlap.interpolate<volSymmTensorField>(mesh_, suppressed);
    overlap.interpolate<volTensorField>(mesh_, suppressed);

    return true;
}


bool Foam::oversetFvMeshBase::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // For postprocessing : write cellTypes and zoneID
    bool ok = true;
    {
        const cellCellStencilObject& overlap = Stencil::New(mesh_);

        const labelUList& cellTypes = overlap.cellTypes();

        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        forAll(volTypes.internalField(), cellI)
        {
            volTypes[cellI] = cellTypes[cellI];
        }
        volTypes.correctBoundaryConditions();
        volTypes.writeObject(streamOpt, writeOnProc);
    }
    {
        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        const cellCellStencilObject& overlap = Stencil::New(mesh_);
        const labelIOList& zoneID = overlap.zoneID();

        forAll(zoneID, cellI)
        {
            volZoneID[cellI] = zoneID[cellI];
        }
        volZoneID.correctBoundaryConditions();
        volZoneID.writeObject(streamOpt, writeOnProc);
    }
    if (debug)
    {
        const cellCellStencilObject& overlap = Stencil::New(mesh_);
        const labelIOList& zoneID = overlap.zoneID();
        const labelListList& cellStencil = overlap.cellStencil();

        // Get remote zones
        labelList donorZoneID(zoneID);
        overlap.cellInterpolationMap().distribute(donorZoneID);

        // Get remote cellCentres
        pointField cc(mesh_.C());
        overlap.cellInterpolationMap().distribute(cc);

        volScalarField volDonorZoneID
        (
            IOobject
            (
                "donorZoneID",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            scalar(-1),
            dimless,
            fvPatchFieldBase::zeroGradientType()
        );

        forAll(cellStencil, cellI)
        {
            const labelList& stencil = cellStencil[cellI];
            if (stencil.size())
            {
                volDonorZoneID[cellI] = donorZoneID[stencil[0]];
                for (label i = 1; i < stencil.size(); i++)
                {
                    if (donorZoneID[stencil[i]] != volDonorZoneID[cellI])
                    {
                        WarningInFunction << "Mixed donor meshes for cell "
                            << cellI << " at " << mesh_.C()[cellI]
                            << " donors:" << UIndirectList<point>(cc, stencil)
                            << endl;
                        volDonorZoneID[cellI] = -2;
                    }
                }
            }
        }
        //- Do not correctBoundaryConditions since re-interpolates!
        //volDonorZoneID.correctBoundaryConditions();
        cellCellStencil::correctBoundaryConditions
        <
            volScalarField,
            oversetFvPatchField<scalar>
        >
        (
            volDonorZoneID
        );
        ok = volDonorZoneID.writeObject(streamOpt, writeOnProc);
    }

    return ok;
}


// ************************************************************************* //
