/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "parPointFieldDistributor.H"
#include "processorPointPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parPointFieldDistributor::verbose_ = 1;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parPointFieldDistributor::parPointFieldDistributor
(
    const pointMesh& srcMesh,
    const bool savePoints,
    const bool isWriteProc
)
:
    srcMesh_(srcMesh),
    nOldPoints_(srcMesh.size()),
    patchMeshPoints_(),
    tgtMeshRef_(nullptr),
    distMapRef_(nullptr),
    patchPointMaps_(),
    dummyHandler_(fileOperation::null()),
    writeHandler_(dummyHandler_),
    isWriteProc_(isWriteProc)
{
    if (savePoints)
    {
        saveMeshPoints();
    }
}


Foam::parPointFieldDistributor::parPointFieldDistributor
(
    const pointMesh& srcMesh,
    const bool savePoints,
    refPtr<fileOperation>& writeHandler
)
:
    srcMesh_(srcMesh),
    nOldPoints_(srcMesh.size()),
    patchMeshPoints_(),
    tgtMeshRef_(nullptr),
    distMapRef_(nullptr),
    patchPointMaps_(),
    dummyHandler_(nullptr),
    writeHandler_(writeHandler),
    isWriteProc_(Switch::INVALID)
{
    if (savePoints)
    {
        saveMeshPoints();
    }
}


Foam::parPointFieldDistributor::parPointFieldDistributor
(
    const pointMesh& srcMesh,
    const pointMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    const bool savePoints,
    const bool isWriteProc
)
:
    srcMesh_(srcMesh),
    nOldPoints_(srcMesh.size()),
    patchMeshPoints_(),
    tgtMeshRef_(tgtMesh),
    distMapRef_(distMap),
    patchPointMaps_(),
    dummyHandler_(fileOperation::null()),
    writeHandler_(dummyHandler_),
    isWriteProc_(isWriteProc)
{
    if (savePoints)
    {
        saveMeshPoints();
    }
}


Foam::parPointFieldDistributor::parPointFieldDistributor
(
    const pointMesh& srcMesh,
    const pointMesh& tgtMesh,
    const mapDistributePolyMesh& distMap,
    const bool savePoints,
    refPtr<fileOperation>& writeHandler
)
:
    srcMesh_(srcMesh),
    nOldPoints_(srcMesh.size()),
    patchMeshPoints_(),
    tgtMeshRef_(tgtMesh),
    distMapRef_(distMap),
    patchPointMaps_(),
    dummyHandler_(nullptr),
    writeHandler_(writeHandler),
    isWriteProc_(Switch::INVALID)
{
    if (savePoints)
    {
        saveMeshPoints();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::parPointFieldDistributor::hasMeshPoints() const
{
    return !patchMeshPoints_.empty();
}


bool Foam::parPointFieldDistributor::hasPatchPointMaps() const
{
    return !patchPointMaps_.empty();
}


bool Foam::parPointFieldDistributor::hasTarget() const
{
    return (tgtMeshRef_ && distMapRef_);
}


void Foam::parPointFieldDistributor::clearMeshPoints()
{
    patchMeshPoints_.clear();
}

void Foam::parPointFieldDistributor::clearPatchPointMaps()
{
    patchPointMaps_.clear();
}


void Foam::parPointFieldDistributor::saveMeshPoints()
{
    const pointBoundaryMesh& patches = srcMesh_.boundary();

    patchMeshPoints_.clear();
    patchMeshPoints_.resize(patches.size());

    forAll(patches, patchi)
    {
        if (!isA<processorPointPatch>(patches[patchi]))
        {
            // Copy meshPoints
            patchMeshPoints_.set
            (
                patchi,
                new labelList(patches[patchi].meshPoints())
            );
        }
    }
}


void Foam::parPointFieldDistributor::createPatchPointMaps()
{
    if (!tgtMeshRef_ || !distMapRef_)
    {
        FatalErrorInFunction
            << "Cannot create maps without target mesh and/or distribution!"
            << abort(FatalError);
    }

    const auto& tgtMesh = tgtMeshRef_();
    const auto& distMap = distMapRef_();

    const auto& newPatches = tgtMesh.boundary();
    const auto& oldPatches = srcMesh_.boundary();

    patchPointMaps_.clear();
    patchPointMaps_.resize(oldPatches.size());

    // if (patchPointMaps_.size() != patchMeshPoints_.size())
    // {
    //     // Warn?
    // }

    forAll(oldPatches, patchi)
    {
        if (!isA<processorPointPatch>(oldPatches[patchi]))
        {
            // Create map for patch points only
            labelList oldToNewSub;
            labelList oldToNewConstruct;

            // Copy point map
            patchPointMaps_.set
            (
                patchi,
                new mapDistributeBase(distMap.pointMap())
            );

            const labelList& oldMeshPoints =
            (
                patchMeshPoints_.test(patchi)
              ? patchMeshPoints_[patchi]
              : oldPatches[patchi].meshPoints()  // <- Questionable!
            );

            patchPointMaps_[patchi].compactData
            (
                oldMeshPoints,
                newPatches[patchi].meshPoints(),
                oldToNewSub,
                oldToNewConstruct,
                nOldPoints_,
                UPstream::msgType()
            );
        }
    }
}


void Foam::parPointFieldDistributor::resetTarget()
{
    tgtMeshRef_.reset(nullptr);
    distMapRef_.reset(nullptr);

    // Old maps are now invalid
    clearPatchPointMaps();
}


void Foam::parPointFieldDistributor::resetTarget
(
    const pointMesh& tgtMesh,
    const mapDistributePolyMesh& distMap
)
{
    tgtMeshRef_.cref(tgtMesh);
    distMapRef_.cref(distMap);

    // Old maps are now invalid
    clearPatchPointMaps();
}


Foam::label Foam::parPointFieldDistributor::distributeAllFields
(
    const IOobjectList& objects,
    const wordRes& selected
) const
{
    label nTotal = 0;

    nTotal += distributePointFields<scalar>(objects, selected);
    nTotal += distributePointFields<vector>(objects, selected);
    nTotal += distributePointFields<symmTensor>(objects, selected);
    nTotal += distributePointFields<sphericalTensor>(objects, selected);
    nTotal += distributePointFields<tensor>(objects, selected);

    return nTotal;
}


// ************************************************************************* //
