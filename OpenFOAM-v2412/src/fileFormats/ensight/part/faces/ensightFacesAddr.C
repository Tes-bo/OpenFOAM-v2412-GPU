/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "ensightFaces.H"
#include "ensightOutput.H"

#include "polyMesh.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::ensightFaces::uniqueMeshPoints
(
    const polyMesh& mesh,
    labelList& uniqueMeshPointLabels,
    bool parallel
) const
{
    const ensightFaces& part = *this;

    parallel = parallel && Pstream::parRun();

    // Renumber the patch points/faces into unique points
    label nPoints = 0;  // Total number of points
    labelList pointToGlobal;  // local point to unique global index

    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();


    // Use the properly sorted faceIds (ensightFaces) and do NOT use
    // the faceZone or anything else directly, otherwise the
    // point-maps will not correspond.
    // - perform face-flipping later

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(faces, part.faceIds()),
        points
    );

    if (parallel)
    {
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        nPoints = globalPointsPtr().totalSize();  // nPoints (global)
    }
    else
    {
        // Non-parallel
        // - all information already available from PrimitivePatch

        nPoints = pp.meshPoints().size();
        uniqueMeshPointLabels = pp.meshPoints();

        // Not needed: pointToGlobal
    }

    return nPoints;
}


// ************************************************************************* //
