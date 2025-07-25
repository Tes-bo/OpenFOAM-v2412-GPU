/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "parFvFieldDistributor.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parFvFieldDistributor::verbose_ = 1;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parFvFieldDistributor::createPatchFaceMaps()
{
    const fvBoundaryMesh& fvb = srcMesh_.boundary();

    patchFaceMaps_.resize(fvb.size());

    forAll(fvb, patchi)
    {
        if (!isA<processorFvPatch>(fvb[patchi]))
        {
            // Create compact map for patch faces only
            // - compact for used faces only (destination patch faces)
            labelList oldToNewSub;
            labelList oldToNewConstruct;

            // Copy face map
            patchFaceMaps_.set
            (
                patchi,
                new mapDistributeBase(distMap_.faceMap())
            );

            patchFaceMaps_[patchi].compactRemoteData
            (
                bitSet(tgtMesh_.boundaryMesh()[patchi].range()),
                oldToNewSub,
                oldToNewConstruct,
                srcMesh_.nFaces(),   // max index of subMap
                UPstream::msgType()
            );
            //Pout<< "patchMap:" << patchFaceMaps_[patchi] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parFvFieldDistributor::parFvFieldDistributor
(
    const fvMesh& srcMesh,
    fvMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    const bool isWriteProc
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap),
    dummyHandler_(fileOperation::null()),
    writeHandler_(dummyHandler_),
    isWriteProc_(isWriteProc)
{
    createPatchFaceMaps();
}


Foam::parFvFieldDistributor::parFvFieldDistributor
(
    const fvMesh& srcMesh,
    fvMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    refPtr<fileOperation>& writeHandler
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap),
    dummyHandler_(nullptr),
    writeHandler_(writeHandler),
    isWriteProc_(Switch::INVALID)
{
    createPatchFaceMaps();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parFvFieldDistributor::reconstructPoints()
{
    // Reconstruct the points for moving mesh cases and write them out
    distributedFieldMapper mapper
    (
        labelUList::null(),
        distMap_.pointMap()
    );

    pointField newPoints(srcMesh_.points(), mapper);
    tgtMesh_.movePoints(newPoints);

    if (isWriteProc_.good())
    {
        if (UPstream::master())
        {
            tgtMesh_.write();
        }
    }
    else if (writeHandler_ && writeHandler_->good())
    {
        auto oldHandler = fileOperation::fileHandler(writeHandler_);
        const label oldComm = UPstream::commWorld(fileHandler().comm());

        tgtMesh_.write();

        writeHandler_  = fileOperation::fileHandler(oldHandler);
        UPstream::commWorld(oldComm);
    }
}


// ************************************************************************* //
