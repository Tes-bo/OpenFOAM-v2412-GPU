/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "mapDistributePolyMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::calcPatchSizes()
{
    oldPatchSizes_.resize_nocopy(oldPatchStarts_.size());

    if (oldPatchStarts_.size())
    {
        // Calculate old patch sizes
        for (label patchi = 0; patchi < oldPatchStarts_.size() - 1; patchi++)
        {
            oldPatchSizes_[patchi] =
                oldPatchStarts_[patchi + 1] - oldPatchStarts_[patchi];
        }

        // Set the last one by hand
        const label lastPatchID = oldPatchStarts_.size() - 1;

        oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

        if (min(oldPatchSizes_) < 0)
        {
            FatalErrorInFunction
                << "Calculated negative old patch size:" << oldPatchSizes_ << nl
                << "Error in mapping data" << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributePolyMesh::mapDistributePolyMesh()
:
    mapDistributePolyMesh(UPstream::worldComm)
{}


Foam::mapDistributePolyMesh::mapDistributePolyMesh(const label comm)
:
    nOldPoints_(0),
    nOldFaces_(0),
    nOldCells_(0),
    oldPatchSizes_(),
    oldPatchStarts_(),
    oldPatchNMeshPoints_(),
    pointMap_(comm),
    faceMap_(comm),
    cellMap_(comm),
    patchMap_(comm)
{}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    const mapDistributePolyMesh& map
)
:
    mapDistributePolyMesh()
{
    deepCopy(map);
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    mapDistributePolyMesh&& map
)
:
    mapDistributePolyMesh()
{
    transfer(map);
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    mapDistributePolyMesh& map,
    bool reuse
)
:
    mapDistributePolyMesh()
{
    if (reuse)
    {
        transfer(map);
    }
    else
    {
        deepCopy(map);
    }
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    const polyMesh& mesh,

    // mesh before changes
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    labelList&& oldPatchStarts,
    labelList&& oldPatchNMeshPoints,

    // how to subset pieces of mesh to send across
    labelListList&& subPointMap,
    labelListList&& subFaceMap,
    labelListList&& subCellMap,
    labelListList&& subPatchMap,

    // how to reconstruct received mesh
    labelListList&& constructPointMap,
    labelListList&& constructFaceMap,
    labelListList&& constructCellMap,
    labelListList&& constructPatchMap,

    const bool subFaceHasFlip,
    const bool constructFaceHasFlip
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(),
    oldPatchStarts_(std::move(oldPatchStarts)),
    oldPatchNMeshPoints_(std::move(oldPatchNMeshPoints)),
    pointMap_
    (
        mesh.nPoints(),
        std::move(subPointMap),
        std::move(constructPointMap)
    ),
    faceMap_
    (
        mesh.nFaces(),
        std::move(subFaceMap),
        std::move(constructFaceMap),
        subFaceHasFlip,
        constructFaceHasFlip
    ),
    cellMap_
    (
        mesh.nCells(),
        std::move(subCellMap),
        std::move(constructCellMap)
    ),
    patchMap_
    (
        mesh.boundaryMesh().size(),
        std::move(subPatchMap),
        std::move(constructPatchMap)
    )
{
    calcPatchSizes();
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    // mesh before changes
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    labelList&& oldPatchStarts,
    labelList&& oldPatchNMeshPoints,

    // how to transfer pieces of mesh
    mapDistribute&& pointMap,
    mapDistribute&& faceMap,
    mapDistribute&& cellMap,
    mapDistribute&& patchMap
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    oldPatchSizes_(),
    oldPatchStarts_(std::move(oldPatchStarts)),
    oldPatchNMeshPoints_(std::move(oldPatchNMeshPoints)),
    pointMap_(std::move(pointMap)),
    faceMap_(std::move(faceMap)),
    cellMap_(std::move(cellMap)),
    patchMap_(std::move(patchMap))
{
    calcPatchSizes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::clear()
{
    nOldPoints_ = 0;
    nOldFaces_ = 0;
    nOldCells_ = 0;
    oldPatchSizes_.clear();
    oldPatchStarts_.clear();
    oldPatchNMeshPoints_.clear();
    pointMap_.clear();
    faceMap_.clear();
    cellMap_.clear();
    patchMap_.clear();
}


void Foam::mapDistributePolyMesh::deepCopy(const mapDistributePolyMesh& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
        return;
    }

    nOldPoints_ = rhs.nOldPoints_;
    nOldFaces_ = rhs.nOldFaces_;
    nOldCells_ = rhs.nOldCells_;
    oldPatchSizes_ = rhs.oldPatchSizes_;
    oldPatchStarts_ = rhs.oldPatchStarts_;
    oldPatchNMeshPoints_ = rhs.oldPatchNMeshPoints_;
    pointMap_ = rhs.pointMap_;
    faceMap_ = rhs.faceMap_;
    cellMap_ = rhs.cellMap_;
    patchMap_ = rhs.patchMap_;
}


void Foam::mapDistributePolyMesh::transfer(mapDistributePolyMesh& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
        return;
    }

    nOldPoints_ = rhs.nOldPoints_;
    nOldFaces_ = rhs.nOldFaces_;
    nOldCells_ = rhs.nOldCells_;
    oldPatchSizes_.transfer(rhs.oldPatchSizes_);
    oldPatchStarts_.transfer(rhs.oldPatchStarts_);
    oldPatchNMeshPoints_.transfer(rhs.oldPatchNMeshPoints_);
    pointMap_.transfer(rhs.pointMap_);
    faceMap_.transfer(rhs.faceMap_);
    cellMap_.transfer(rhs.cellMap_);
    patchMap_.transfer(rhs.patchMap_);

    rhs.nOldPoints_ = 0;
    rhs.nOldFaces_ = 0;
    rhs.nOldCells_ = 0;
}


void Foam::mapDistributePolyMesh::distributePointIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        ListOps::createWithValue<bool>(nOldPoints(), lst, true, false)
    );

    // Distribute
    distributePointData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::mapDistributePolyMesh::distributeFaceIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        ListOps::createWithValue<bool>(nOldFaces(), lst, true, false)
    );

    // Distribute
    distributeFaceData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::mapDistributePolyMesh::distributeCellIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        ListOps::createWithValue<bool>(nOldCells(), lst, true, false)
    );

    // Distribute
    distributeCellData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


void Foam::mapDistributePolyMesh::distributePatchIndices(labelList& lst) const
{
    // Construct boolList from selected elements
    boolList isSelected
    (
        ListOps::createWithValue<bool>
        (
            oldPatchStarts().size(), // nOldPatches
            lst,
            true, // set value
            false // default value
        )
    );

    // Distribute
    distributePatchData(isSelected);

    // Collect selected elements
    lst = findIndices(isSelected, true);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::operator=(const mapDistributePolyMesh& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
        return;
    }

    nOldPoints_ = rhs.nOldPoints_;
    nOldFaces_ = rhs.nOldFaces_;
    nOldCells_ = rhs.nOldCells_;
    oldPatchSizes_ = rhs.oldPatchSizes_;
    oldPatchStarts_ = rhs.oldPatchStarts_;
    oldPatchNMeshPoints_ = rhs.oldPatchNMeshPoints_;
    pointMap_ = rhs.pointMap_;
    faceMap_ = rhs.faceMap_;
    cellMap_ = rhs.cellMap_;
    patchMap_ = rhs.patchMap_;
}


void Foam::mapDistributePolyMesh::operator=(mapDistributePolyMesh&& rhs)
{
    if (this != &rhs)
    {
        // Avoid self assignment
        transfer(rhs);
    }
}


// ************************************************************************* //
