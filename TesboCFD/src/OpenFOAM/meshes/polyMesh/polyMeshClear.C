/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "globalMeshData.H"
#include "MeshObject.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::removeBoundary()
{
    DebugInFunction << "Removing boundary patches." << endl;

    boundary_.clear();

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::clearGeom()
{
    DebugInFunction << "Clearing geometric data" << endl;

    // Clear all geometric mesh objects
    meshObject::clear<pointMesh, GeometricMeshObject>(*this);
    meshObject::clear<polyMesh, GeometricMeshObject>(*this);

    primitiveMesh::clearGeom();

    boundary_.clearGeom();

    // Reset valid directions (could change with rotation)
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Remove the cell tree
    cellTreePtr_.reset(nullptr);
}


void Foam::polyMesh::updateGeomPoints
(
    pointIOField&& newPoints,
    autoPtr<labelIOList>& newTetBasePtIsPtr
)
{
    DebugInFunction
        << "Updating geometric data with newPoints:"
        << newPoints.size() << " newTetBasePtIs:"
        << bool(newTetBasePtIsPtr) << endl;

    if (points_.size() != 0 && points_.size() != newPoints.size())
    {
        FatalErrorInFunction
            << "Point motion detected but number of points "
            << newPoints.size() << " in "
            << newPoints.objectPath() << " does not correspond to "
            << " current " << points_.size()
            << exit(FatalError);
    }

    // Clear all geometric mesh objects that are not 'movable'
    meshObject::clearUpto
    <
        pointMesh,
        TopologicalMeshObject,
        MoveableMeshObject
    >
    (
        *this
    );
    meshObject::clearUpto
    <
        polyMesh,
        TopologicalMeshObject,
        MoveableMeshObject
    >
    (
        *this
    );

    primitiveMesh::clearGeom();

    boundary_.clearGeom();

    // Reset valid directions (could change with rotation)
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Remove the cell tree
    cellTreePtr_.reset(nullptr);

    // Update local data
    points_.instance() = newPoints.instance();
    points_.transfer(newPoints);

    // Optional new tet base points
    if (newTetBasePtIsPtr)
    {
        tetBasePtIsPtr_ = std::move(newTetBasePtIsPtr);
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    // Derived info
    bounds_ = boundBox(points_);

    // Rotation can cause direction vector to change
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Update all 'movable' objects
    meshObject::movePoints<polyMesh>(*this);
    meshObject::movePoints<pointMesh>(*this);
}


void Foam::polyMesh::clearAddressing(const bool isMeshUpdate)
{
    DebugInFunction
        << "isMeshUpdate:" << isMeshUpdate << endl;

    if (isMeshUpdate)
    {
        // Part of a mesh update.
        // Keep meshObjects that have an updateMesh callback
        meshObject::clearUpto
        <
            pointMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
        <
            polyMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObject::clear<pointMesh, TopologicalMeshObject>(*this);
        meshObject::clear<polyMesh, TopologicalMeshObject>(*this);
    }

    primitiveMesh::clearAddressing();

    // parallelData depends on the processorPatch ordering so force
    // recalculation
    globalMeshDataPtr_.reset(nullptr);

    // Reset valid directions
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.reset(nullptr);

    // Remove the cell tree
    cellTreePtr_.reset(nullptr);
}


void Foam::polyMesh::clearPrimitives()
{
    resetMotion();

    points_.clear();
    faces_.clear();
    owner_.clear();
    neighbour_.clear();

    clearedPrimitives_ = true;
}


void Foam::polyMesh::clearOut(const bool isMeshUpdate)
{
    clearGeom();  // not implementable?  isMeshUpdate
    clearAddressing(isMeshUpdate);
}


void Foam::polyMesh::clearTetBasePtIs()
{
    DebugInFunction << "Clearing tet base points" << endl;

    tetBasePtIsPtr_.reset(nullptr);
}


void Foam::polyMesh::clearCellTree()
{
    DebugInFunction << endl;

    cellTreePtr_.reset(nullptr);
}


// ************************************************************************* //
