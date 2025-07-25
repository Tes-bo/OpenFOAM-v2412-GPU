/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
#include "HashOps.H"
#include "emptyPolyPatch.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::polyTopoChange::reorder
(
    const labelUList& oldToNew,
    DynamicList<Type>& lst
)
{
    // Create copy
    DynamicList<Type> oldLst(lst);

    forAll(oldLst, i)
    {
        const label newIdx = oldToNew[i];

        if (newIdx >= 0)
        {
            lst[newIdx] = std::move(oldLst[i]);
        }
    }
}


template<class Type>
void Foam::polyTopoChange::renumberKey
(
    const labelUList& oldToNew,
    Map<Type>& map
)
{
    Map<Type> newMap(map.capacity());

    forAllConstIters(map, iter)
    {
        const label newKey = oldToNew[iter.key()];

        if (newKey >= 0)
        {
            newMap.insert(newKey, iter.val());
        }
    }

    map.transfer(newMap);
}


template<class Type>
Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::makeMesh
(
    autoPtr<Type>& newMeshPtr,
    const IOobject& io,
    const polyMesh& mesh,
    const labelUList& patchMap,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (debug)
    {
        Pout<< "polyTopoChange::changeMesh"
            << "(autoPtr<fvMesh>&, const IOobject&, const fvMesh&"
            << ", const bool, const bool, const bool)"
            << endl;
    }

    if (debug)
    {
        Pout<< "Old mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // new mesh points
    pointField newPoints;
    // number of internal points
    label nInternalPoints;
    // patch slicing
    labelList patchSizes;
    labelList patchStarts;
    // inflate maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromPoints;
    List<objectMap> facesFromEdges;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromPoints;
    List<objectMap> cellsFromEdges;
    List<objectMap> cellsFromFaces;
    List<objectMap> cellsFromCells;

    // old mesh info
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label>> oldFaceZoneMeshPointMaps;

    // Compact, reorder patch faces and calculate mesh/patch maps.
    compactAndReorder
    (
        mesh,
        patchMap,      // from new to old patch
        syncParallel,
        orderCells,
        orderPoints,

        nInternalPoints,
        newPoints,
        patchSizes,
        patchStarts,
        pointsFromPoints,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchStarts,
        oldFaceZoneMeshPointMaps
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Create the mesh
    // ~~~~~~~~~~~~~~~

    //IOobject noReadIO(io);
    //noReadIO.readOpt(IOobject::NO_READ);
    //noReadIO.writeOpt(IOobject::AUTO_WRITE);
    newMeshPtr.reset
    (
        new Type
        (
            io,         //noReadIO
            std::move(newPoints),
            std::move(faces_),
            std::move(faceOwner_),
            std::move(faceNeighbour_)
        )
    );
    Type& newMesh = *newMeshPtr;

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nAdd, nInflate, nMerge, nRemove;
        countMap(pointMap_, reversePointMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Points:"
            << "  added(from point):" << nAdd
            << "  added(from nothing):" << nInflate
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other cell):" << nMerge
            << "  removed:" << nRemove
            << nl
            << endl;
    }


    {
        const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

        polyPatchList newBoundary(patchMap.size());

        forAll(patchMap, patchi)
        {
            const label oldPatchi = patchMap[patchi];

            if (oldPatchi != -1)
            {
                newBoundary.set
                (
                    patchi,
                    oldPatches[oldPatchi].clone
                    (
                        newMesh.boundaryMesh(),
                        patchi,
                        patchSizes[patchi],
                        patchStarts[patchi]
                    )
                );
            }
            else
            {
                // Added patch
                newBoundary.set
                (
                    patchi,
                    new emptyPolyPatch
                    (
                        "patch" + Foam::name(patchi),
                        patchSizes[patchi],
                        patchStarts[patchi],
                        patchi,
                        newMesh.boundaryMesh(),
                        word::null
                    )
                );
            }
        }
        newMesh.addFvPatches(newBoundary);
    }


    // Zones
    // ~~~~~

    // Start off from empty zones.
    const pointZoneMesh& oldPointZones = mesh.pointZones();
    List<pointZone*> pZonePtrs(oldPointZones.size());
    {
        forAll(oldPointZones, i)
        {
            pZonePtrs[i] = new pointZone
            (
                oldPointZones[i].name(),
                i,
                newMesh.pointZones()
            );
        }
    }

    const faceZoneMesh& oldFaceZones = mesh.faceZones();
    List<faceZone*> fZonePtrs(oldFaceZones.size());
    {
        forAll(oldFaceZones, i)
        {
            fZonePtrs[i] = new faceZone
            (
                oldFaceZones[i].name(),
                i,
                newMesh.faceZones()
            );
        }
    }

    const cellZoneMesh& oldCellZones = mesh.cellZones();
    List<cellZone*> cZonePtrs(oldCellZones.size());
    {
        forAll(oldCellZones, i)
        {
            cZonePtrs[i] = new cellZone
            (
                oldCellZones[i].name(),
                i,
                newMesh.cellZones()
            );
        }
    }

    newMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    // Inverse of point/face/cell zone addressing.
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, newMesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    // Clear zone info
    {
        pointZone_.clearStorage();
        faceZone_.clearStorage();
        faceZoneFlip_.clearStorage();
        cellZone_.clearStorage();
    }

    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(newMesh.boundaryMesh().size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        patchMap,
        newMesh.boundaryMesh(),
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(newMesh.faceZones().size());
    calcFaceZonePointMap(newMesh, oldFaceZoneMeshPointMaps, faceZonePointMap);

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(newMesh, Pout);
    }

    labelHashSet flipFaceFluxSet(HashSetOps::used(flipFaceFlux_));

    return autoPtr<mapPolyMesh>::New
    (
        newMesh,
        nOldPoints,
        nOldFaces,
        nOldCells,

        pointMap_,
        pointsFromPoints,

        faceMap_,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,

        cellMap_,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,

        reversePointMap_,
        reverseFaceMap_,
        reverseCellMap_,

        flipFaceFluxSet,

        patchPointMap,

        pointZoneMap,

        faceZonePointMap,
        faceZoneFaceMap,
        cellZoneMap,

        newPoints,          // if empty signals no inflation.
        oldPatchStarts,
        oldPatchNMeshPoints,
        oldCellVolumes,
        true                // steal storage.
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


template<class Type>
Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::makeMesh
(
    autoPtr<Type>& newMeshPtr,
    const IOobject& io,
    const polyMesh& mesh,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    return makeMesh
    (
        newMeshPtr,
        io,
        mesh,
        identity(mesh.boundaryMesh().size()),
        syncParallel,
        orderCells,
        orderPoints
    );
}


// ************************************************************************* //
