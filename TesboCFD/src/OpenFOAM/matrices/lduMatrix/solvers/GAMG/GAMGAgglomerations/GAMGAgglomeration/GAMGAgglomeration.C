/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "Time.H"
#include "GAMGInterface.H"
#include "GAMGProcAgglomeration.H"
#include "pairGAMGAgglomeration.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
    defineRunTimeSelectionTable(GAMGAgglomeration, geometry);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels
(
    const label nCreatedLevels,
    const bool doProcessorAgglomerate
)
{
    nCells_.setSize(nCreatedLevels, 0);
    restrictAddressing_.setSize(nCreatedLevels);
    nFaces_.setSize(nCreatedLevels, 0);
    faceRestrictAddressing_.setSize(nCreatedLevels);
    faceFlipMap_.setSize(nCreatedLevels);
    nPatchFaces_.setSize(nCreatedLevels);
    patchFaceRestrictAddressing_.setSize(nCreatedLevels);
    meshLevels_.setSize(nCreatedLevels);

    // Have procCommunicator_ always, even if not procAgglomerating
    procCommunicator_.setSize(nCreatedLevels + 1, -1);
    if (doProcessorAgglomerate && processorAgglomerate())
    {
        procAgglomMap_.setSize(nCreatedLevels);
        agglomProcIDs_.setSize(nCreatedLevels);
        procAgglomCommunicator_.setSize(nCreatedLevels);
        procCellOffsets_.setSize(nCreatedLevels);
        procFaceMap_.setSize(nCreatedLevels);
        procBoundaryMap_.setSize(nCreatedLevels);
        procBoundaryFaceMap_.setSize(nCreatedLevels);

        procAgglomeratorPtr_().agglomerate();
    }
}


void Foam::GAMGAgglomeration::printLevels() const
{
    Info<< "GAMGAgglomeration:" << nl
        << "    local agglomerator     : " << type() << nl;
    if (processorAgglomerate())
    {
        Info<< "    processor agglomerator : "
            << procAgglomeratorPtr_().type() << nl
            << nl;
    }

    Info<< setw(36) << "nCells"
        << setw(20) << "nFaces/nCells"
        << setw(20) << "nInterfaces"
        << setw(20) << "nIntFaces/nCells"
        << setw(12) << "profile"
        << nl
        << setw(8) << "Level"
        << setw(8) << "nProcs"
        << "    "
        << setw(8) << "avg"
        << setw(8) << "max"
        << "    "
        << setw(8) << "avg"
        << setw(8) << "max"
        << "    "
        << setw(8) << "avg"
        << setw(8) << "max"
        << "    "
        << setw(8) << "avg"
        << setw(8) << "max"
        //<< "    "
        << setw(12) << "avg"
        << nl
        << setw(8) << "-----"
        << setw(8) << "------"
        << "    "
        << setw(8) << "---"
        << setw(8) << "---"
        << "    "
        << setw(8) << "---"
        << setw(8) << "---"
        << "    "
        << setw(8) << "---"
        << setw(8) << "---"
        << "    "
        << setw(8) << "---"
        << setw(8) << "---"
        //<< "    "
        << setw(12) << "---"
        //<< "    "
        << nl;

    const label maxSize = returnReduce(size(), maxOp<label>());

    for (label levelI = 0; levelI <= maxSize; levelI++)
    {
        label nProcs = 0;
        label nCells = 0;
        scalar faceCellRatio = 0;
        label nInterfaces = 0;
        label nIntFaces = 0;
        scalar ratio = 0.0;
        scalar profile = 0.0;

        if (hasMeshLevel(levelI))
        {
            nProcs = 1;

            const lduMesh& fineMesh = meshLevel(levelI);
            nCells = fineMesh.lduAddr().size();
            faceCellRatio =
                scalar(fineMesh.lduAddr().lowerAddr().size())/nCells;

            const lduInterfacePtrsList interfaces =
                fineMesh.interfaces();
            forAll(interfaces, i)
            {
                if (interfaces.set(i))
                {
                    nInterfaces++;
                    nIntFaces += interfaces[i].faceCells().size();
                }
            }
            ratio = scalar(nIntFaces)/nCells;

            profile = fineMesh.lduAddr().band().second();
        }

        label totNprocs = returnReduce(nProcs, sumOp<label>());

        label maxNCells = returnReduce(nCells, maxOp<label>());
        label totNCells = returnReduce(nCells, sumOp<label>());

        scalar maxFaceCellRatio =
            returnReduce(faceCellRatio, maxOp<scalar>());
        scalar totFaceCellRatio =
            returnReduce(faceCellRatio, sumOp<scalar>());

        label maxNInt = returnReduce(nInterfaces, maxOp<label>());
        label totNInt = returnReduce(nInterfaces, sumOp<label>());

        scalar maxRatio = returnReduce(ratio, maxOp<scalar>());
        scalar totRatio = returnReduce(ratio, sumOp<scalar>());

        scalar totProfile = returnReduce(profile, sumOp<scalar>());

        const int oldPrecision = Info.stream().precision(4);

        Info<< setw(8) << levelI
            << setw(8) << totNprocs
            << "    "
            << setw(8) << totNCells/totNprocs
            << setw(8) << maxNCells
            << "    "
            << setw(8) << totFaceCellRatio/totNprocs
            << setw(8) << maxFaceCellRatio
            << "    "
            << setw(8) << scalar(totNInt)/totNprocs
            << setw(8) << maxNInt
            << "    "
            << setw(8) << totRatio/totNprocs
            << setw(8) << maxRatio
            << setw(12) << totProfile/totNprocs
            << nl;

        Info.stream().precision(oldPrecision);
    }
    Info<< endl;
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
    const label nCellsInCoarsestLevel,
    const label nFineCells,
    const label nCoarseCells,
    const label comm
) const
{
    const label nTotalCoarseCells =
        returnReduce(nCoarseCells, sumOp<label>(), UPstream::msgType(), comm);
    if (nTotalCoarseCells < Pstream::nProcs(comm)*nCellsInCoarsestLevel)
    {
        return false;
    }
    else
    {
        const label nTotalFineCells =
            returnReduce(nFineCells, sumOp<label>(), UPstream::msgType(), comm);
        return nTotalCoarseCells < nTotalFineCells;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    MeshObject_type(mesh),

    maxLevels_(50),

    updateInterval_(controlDict.getOrDefault<label>("updateInterval", 1)),
    requireUpdate_(false),

    nCellsInCoarsestLevel_
    (
        controlDict.getOrDefault<label>("nCellsInCoarsestLevel", 10)
    ),
    meshInterfaces_(mesh.interfaces()),
    procAgglomeratorPtr_
    (
        (
            (UPstream::nProcs(mesh.comm()) > 1)
         && controlDict.found("processorAgglomerator")
        )
      ? GAMGProcAgglomeration::New
        (
            controlDict.get<word>("processorAgglomerator"),
            *this,
            controlDict
        )
      : autoPtr<GAMGProcAgglomeration>()
    ),

    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    nFaces_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),
    faceFlipMap_(maxLevels_),
    nPatchFaces_(maxLevels_),
    patchFaceRestrictAddressing_(maxLevels_),

    meshLevels_(maxLevels_)
{
    // Limit the cells in the coarsest level based on the local number of
    // cells.  Note: 2 for pair-wise
    nCellsInCoarsestLevel_ =
        max(1, min(mesh.lduAddr().size()/2, nCellsInCoarsestLevel_));

    // Ensure all procs see the same nCellsInCoarsestLevel_
    reduce(nCellsInCoarsestLevel_, minOp<label>());

    procCommunicator_.setSize(maxLevels_ + 1, -1);
    if (processorAgglomerate())
    {
        procAgglomMap_.setSize(maxLevels_);
        agglomProcIDs_.setSize(maxLevels_);
        procAgglomCommunicator_.setSize(maxLevels_);
        procCellOffsets_.setSize(maxLevels_);
        procFaceMap_.setSize(maxLevels_);
        procBoundaryMap_.setSize(maxLevels_);
        procBoundaryFaceMap_.setSize(maxLevels_);
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
{
    const GAMGAgglomeration* agglomPtr =
        mesh.thisDb().cfindObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );

    if (agglomPtr)
    {
        if (agglomPtr->requireUpdate_)
        {
            mesh.thisDb().checkOut(const_cast<GAMGAgglomeration*>(agglomPtr));
            return GAMGAgglomeration::New(mesh, controlDict);
        }

        return *agglomPtr;
    }

    {
        const word agglomeratorType
        (
            controlDict.getOrDefault<word>("agglomerator", "faceAreaPair")
        );

        mesh.thisDb().time().libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            lduMeshConstructorTablePtr_
        );

        auto* ctorPtr = lduMeshConstructorTable(agglomeratorType);

        if (!ctorPtr)
        {
            FatalErrorInFunction
                << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid matrix GAMGAgglomeration types :"
                << lduMatrixConstructorTablePtr_->sortedToc() << endl
                << "Valid geometric GAMGAgglomeration types :"
                << lduMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        auto agglomPtr(ctorPtr(mesh, controlDict));
        if (debug)
        {
            agglomPtr().printLevels();
        }
        return regIOobject::store(agglomPtr);
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
{
    const lduMesh& mesh = matrix.mesh();

    const GAMGAgglomeration* agglomPtr =
        mesh.thisDb().cfindObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );

    if (agglomPtr)
    {
        if (agglomPtr->requireUpdate_)
        {
            mesh.thisDb().checkOut(const_cast<GAMGAgglomeration*>(agglomPtr));
            return GAMGAgglomeration::New(matrix, controlDict);
        }

        return *agglomPtr;
    }

    {
        const word agglomeratorType
        (
            controlDict.getOrDefault<word>("agglomerator", "faceAreaPair")
        );

        mesh.thisDb().time().libs().open
        (
            controlDict,
            "algebraicGAMGAgglomerationLibs",
            lduMatrixConstructorTablePtr_
        );

        auto* ctorPtr = lduMatrixConstructorTable(agglomeratorType);

        if (!ctorPtr)
        {
            return New(mesh, controlDict);
        }
        else
        {
            auto agglomPtr(ctorPtr(matrix, controlDict));
            if (debug)
            {
                agglomPtr().printLevels();
            }
            return regIOobject::store(agglomPtr);
        }
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
{

    const GAMGAgglomeration* agglomPtr =
        mesh.thisDb().cfindObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );

    if (agglomPtr)
    {
        if (agglomPtr->requireUpdate_)
        {
            mesh.thisDb().checkOut(const_cast<GAMGAgglomeration*>(agglomPtr));
            return GAMGAgglomeration::New(mesh, controlDict);
        }

        return *agglomPtr;
    }

    {
        const word agglomeratorType
        (
            controlDict.getOrDefault<word>("agglomerator", "faceAreaPair")
        );

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            geometryConstructorTablePtr_
        );

        auto* ctorPtr = geometryConstructorTable(agglomeratorType);

        if (!ctorPtr)
        {
            FatalErrorInFunction
                << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid geometric GAMGAgglomeration types :"
                << geometryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        auto agglomPtr
        (
            ctorPtr
            (
                mesh,
                cellVolumes,
                faceAreas,
                controlDict
            )
        );
        if (debug)
        {
            agglomPtr().printLevels();
        }
        return regIOobject::store(agglomPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::GAMGAgglomeration::requiresUpdate() const
{
    return
        (updateInterval_ > 0)
     && ((mesh_.thisDb().time().timeIndex() % updateInterval_) == 0);
}


bool Foam::GAMGAgglomeration::movePoints()
{
    if (requiresUpdate())
    {
        requireUpdate_ = true;
    }
    // What to return?
    return requireUpdate_;
}


const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return mesh_;
    }
    else
    {
        return meshLevels_[i - 1];
    }
}


bool Foam::GAMGAgglomeration::hasMeshLevel(const label i) const
{
    if (i == 0)
    {
        return true;
    }
    else
    {
        return meshLevels_.set(i - 1);
    }
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return meshInterfaces_;
    }
    else
    {
        return meshLevels_[i - 1].rawInterfaces();
    }
}


void Foam::GAMGAgglomeration::clearLevel(const label i)
{
    if (hasMeshLevel(i))
    {
        meshLevels_.set(i - 1, nullptr);

        if (i < nCells_.size())
        {
            nCells_[i] = -555;
            restrictAddressing_.set(i, nullptr);
            nFaces_[i] = -666;
            faceRestrictAddressing_.set(i, nullptr);
            faceFlipMap_.set(i, nullptr);
            nPatchFaces_.set(i, nullptr);
            patchFaceRestrictAddressing_.set(i, nullptr);
        }
    }
}


const Foam::labelList& Foam::GAMGAgglomeration::procAgglomMap
(
    const label leveli
) const
{
    return procAgglomMap_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::agglomProcIDs
(
    const label leveli
) const
{
    return agglomProcIDs_[leveli];
}


bool Foam::GAMGAgglomeration::hasProcMesh(const label leveli) const
{
    return procCommunicator_[leveli] != -1;
}


Foam::label Foam::GAMGAgglomeration::procCommunicator(const label leveli) const
{
    return procCommunicator_[leveli];
}


Foam::label Foam::GAMGAgglomeration::agglomCommunicator(const label leveli) const
{
    return procAgglomCommunicator_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::cellOffsets
(
    const label leveli
) const
{
    return procCellOffsets_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::faceMap
(
    const label leveli
) const
{
    return procFaceMap_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::boundaryMap
(
    const label leveli
) const
{
    return procBoundaryMap_[leveli];
}


const Foam::labelListListList& Foam::GAMGAgglomeration::boundaryFaceMap
(
    const label leveli
) const
{
    return procBoundaryFaceMap_[leveli];
}


bool Foam::GAMGAgglomeration::checkRestriction
(
    labelList& newRestrict,
    label& nNewCoarse,
    const lduAddressing& fineAddressing,
    const labelUList& restriction,
    const label nCoarse
)
{
    if (fineAddressing.size() != restriction.size())
    {
        FatalErrorInFunction
            << "nCells:" << fineAddressing.size()
            << " agglom:" << restriction.size()
            << abort(FatalError);
    }

    // Seed (master) for every region
    labelList master(identity(fineAddressing.size()));

    // Now loop and transport master through region
    const labelUList& lower = fineAddressing.lowerAddr();
    const labelUList& upper = fineAddressing.upperAddr();

    while (true)
    {
        label nChanged = 0;

        forAll(lower, facei)
        {
            const label own = lower[facei];
            const label nei = upper[facei];

            if (restriction[own] == restriction[nei])
            {
                // coarse-mesh-internal face

                if (master[own] < master[nei])
                {
                    master[nei] = master[own];
                    nChanged++;
                }
                else if (master[own] > master[nei])
                {
                    master[own] = master[nei];
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }


    // Count number of regions/masters per coarse cell
    labelListList coarseToMasters(nCoarse);
    nNewCoarse = 0;
    forAll(restriction, celli)
    {
        labelList& masters = coarseToMasters[restriction[celli]];

        if (!masters.found(master[celli]))
        {
            masters.append(master[celli]);
            nNewCoarse++;
        }
    }

    if (nNewCoarse > nCoarse)
    {
        //WarningInFunction
        //    << "Have " << nCoarse
        //    << " agglomerated cells but " << nNewCoarse
        //    << " disconnected regions" << endl;

        // Keep coarseToMasters[0] the original coarse, allocate new ones
        // for the others
        labelListList coarseToNewCoarse(coarseToMasters.size());

        nNewCoarse = nCoarse;

        forAll(coarseToMasters, coarseI)
        {
            const labelList& masters = coarseToMasters[coarseI];

            labelList& newCoarse = coarseToNewCoarse[coarseI];
            newCoarse.setSize(masters.size());
            newCoarse[0] = coarseI;
            for (label i=1; i<newCoarse.size(); i++)
            {
                newCoarse[i] = nNewCoarse++;
            }
        }

        newRestrict.setSize(fineAddressing.size());
        forAll(restriction, celli)
        {
            const label coarseI = restriction[celli];

            const label index = coarseToMasters[coarseI].find(master[celli]);
            newRestrict[celli] = coarseToNewCoarse[coarseI][index];
        }

        return false;
    }

    return true;
}


// ************************************************************************* //
