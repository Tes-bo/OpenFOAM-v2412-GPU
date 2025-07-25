/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "PatchTools.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PatchTools::gatherAndMerge
(
    const scalar mergeDist,
    const PrimitivePatch<FaceList, PointField>& pp,
    Field
    <
        typename PrimitivePatch<FaceList, PointField>::point_type
    >& mergedPoints,
    List
    <
        typename PrimitivePatch<FaceList, PointField>::face_type
    >& mergedFaces,
    globalIndex& pointAddr,
    globalIndex& faceAddr,
    labelList& pointMergeMap,
    const bool useLocal
)
{
    typedef typename PrimitivePatch<FaceList, PointField>::face_type FaceType;

    // Faces from all ranks
    faceAddr.reset(globalIndex::gatherOnly{}, pp.size());

    // Points from all ranks
    pointAddr.reset
    (
        globalIndex::gatherOnly{},
        (useLocal ? pp.localPoints().size() : pp.points().size())
    );

    if (useLocal)
    {
        faceAddr.gather(pp.localFaces(), mergedFaces);
        pointAddr.gather(pp.localPoints(), mergedPoints);
    }
    else
    {
        faceAddr.gather(pp, mergedFaces);
        pointAddr.gather(pp.points(), mergedPoints);
    }

    // Relabel faces according to global point offsets
    for (const label proci : faceAddr.subProcs())
    {
        SubList<FaceType> slot(mergedFaces, faceAddr.range(proci));

        for (auto& f : slot)
        {
            pointAddr.inplaceToGlobal(proci, f);
        }
    }


    // Merging points
    label nPointsChanged(0);

    labelList boundaryPoints;

    if (UPstream::parRun())
    {
        const globalIndex localPointAddr
        (
            globalIndex::gatherOnly{},
            pp.localPoints().size()
        );

        const globalIndex bndPointAddr
        (
            globalIndex::gatherOnly{},
            pp.boundaryPoints().size()
        );

        bndPointAddr.gather(pp.boundaryPoints(), boundaryPoints);

        // Relabel according to global point offsets
        for (const label proci : localPointAddr.subProcs())
        {
            SubList<label> slot(boundaryPoints, bndPointAddr.range(proci));
            localPointAddr.inplaceToGlobal(proci, slot);
        }
    }


    if (UPstream::parRun() && UPstream::master())
    {
        labelList pointToUnique;

        nPointsChanged = Foam::inplaceMergePoints
        (
            mergedPoints,
            boundaryPoints,  // selection of points to merge
            mergeDist,
            false,           // verbose = false
            pointToUnique
        );

        if (nPointsChanged)
        {
            // Renumber faces to use unique point numbers
            for (auto& f : mergedFaces)
            {
                inplaceRenumber(pointToUnique, f);
            }

            // Store point mapping
            if (notNull(pointMergeMap))
            {
                pointMergeMap.transfer(pointToUnique);
            }
        }
    }

    if (!nPointsChanged && notNull(pointMergeMap))
    {
        // Safety
        pointMergeMap = identity(mergedPoints.size());
    }
}


template<class FaceList, class PointField>
void Foam::PatchTools::gatherAndMerge
(
    const scalar mergeDist,
    const PrimitivePatch<FaceList, PointField>& pp,
    Field
    <
        typename PrimitivePatch<FaceList, PointField>::point_type
    >& mergedPoints,
    List
    <
        typename PrimitivePatch<FaceList, PointField>::face_type
    >& mergedFaces,
    labelList& pointMergeMap,
    const bool useLocal
)
{
    globalIndex pointAddr;
    globalIndex faceAddr;

    PatchTools::gatherAndMerge<FaceList, PointField>
    (
        mergeDist,
        pp,
        mergedPoints,
        mergedFaces,
        pointAddr,
        faceAddr,
        pointMergeMap,
        useLocal
    );
}


template<class FaceList>
void Foam::PatchTools::gatherAndMerge
(
    const polyMesh& mesh,
    const FaceList& localFaces,
    const labelList& meshPoints,
    const Map<label>& meshPointMap,

    labelList& pointToGlobal,
    labelList& uniqueMeshPointLabels,
    autoPtr<globalIndex>& globalPointsPtr,
    autoPtr<globalIndex>& globalFacesPtr,
    List<typename FaceList::value_type>& mergedFaces,
    pointField& mergedPoints
)
{
    typedef typename FaceList::value_type FaceType;

    if (UPstream::parRun())
    {
        // Renumber the points/faces into unique points
        globalPointsPtr = mesh.globalData().mergePoints
        (
            meshPoints,
            meshPointMap,
            pointToGlobal,
            uniqueMeshPointLabels
        );

        globalFacesPtr.reset(new globalIndex(localFaces.size()));

        // Renumber faces locally
        List<FaceType> myFaces(localFaces);
        for (auto& f : myFaces)
        {
            inplaceRenumber(pointToGlobal, f);
        }

        // Can also use
        //     UIndirectList<point>(mesh.points(), uniqueMeshPointLabels)
        // but favour communication over local memory use
        globalPointsPtr().gather
        (
            pointField(mesh.points(), uniqueMeshPointLabels),
            mergedPoints
        );
        globalFacesPtr().gather(myFaces, mergedFaces);
    }
    else
    {
        pointToGlobal = identity(meshPoints.size());
        uniqueMeshPointLabels = pointToGlobal;

        globalPointsPtr.reset(new globalIndex(meshPoints.size()));
        globalFacesPtr.reset(new globalIndex(localFaces.size()));

        mergedFaces = localFaces;
        mergedPoints = pointField(mesh.points(), meshPoints);
    }
}


// ************************************************************************* //
