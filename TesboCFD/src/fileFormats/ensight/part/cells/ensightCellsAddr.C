/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "ensightCells.H"
#include "polyMesh.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "ListOps.H"
#include "manifoldCellsMeshObject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Map<Foam::label>
Foam::ensightCells::meshPointMap(const polyMesh& mesh) const
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();
    const faceList& meshFaces = mesh.faces();

    Map<label> pointMap(8*this->size());

    // Pass 1: markup used points from cells

    for (const label celli : this->cellIds())
    {
        for (const label facei : meshCells[celli])
        {
            for (const label pointi : meshFaces[facei])
            {
                pointMap.insert(pointi, 0);
            }
        }
    }

    // Compact point numbering, preserves the original order
    label nPoints = 0;
    for (const label pointi : pointMap.sortedToc())
    {
        pointMap(pointi) = nPoints++;
    }

    return pointMap;
}


Foam::label Foam::ensightCells::meshPointMapppings
(
    const polyMesh& mesh,
    labelList& pointToGlobalRequest,
    labelList& uniqueMeshPointLabels,
    bool parallel
) const
{
    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();

    labelList pointToGlobal;

    const bool rewritePointMap = notNull(pointToGlobalRequest);

    if (notNull(pointToGlobalRequest))
    {
        pointToGlobal.transfer(pointToGlobalRequest);
    }


    const ensightCells& part = *this;

    parallel = parallel && Pstream::parRun();

    // Renumber the points/faces into unique points
    label nPoints = 0;  // Total number of points

    bool allCells = (part.size() == mesh.nCells());

    if (parallel)
    {
        Pstream::reduceAnd(allCells);

        if (allCells)
        {
            // All cells used, and thus all points

            autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

            nPoints = globalPointsPtr().totalSize();   // nPoints (global)
        }
        else
        {
            // Map mesh point index to local (compact) point index

            Map<label> meshPointMap(part.meshPointMap(mesh));

            labelList meshPoints(meshPointMap.sortedToc());

            autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    meshPoints,
                    meshPointMap,
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

            nPoints = globalPointsPtr().totalSize();   // nPoints (global)

            meshPointMap.clear();

            // The mergePoints returns pointToGlobal under the
            // assumption of local addressing
            // (eg, patch localFaces).
            // Recast as original mesh points to new global points

            if (rewritePointMap)
            {
                labelList oldToNew(mesh.nPoints(), -1);

                forAll(meshPoints, i)
                {
                    const label orig = meshPoints[i];
                    const label glob = pointToGlobal[i];

                    oldToNew[orig] = glob;
                }

                pointToGlobal.transfer(oldToNew);
            }
        }
    }
    else
    {
        // Non-parallel

        nPoints = mesh.nPoints();
        pointToGlobal.resize_nocopy(nPoints);

        if (allCells)
        {
            // All cells used, and thus all points

            uniqueMeshPointLabels.resize_nocopy(nPoints);

            Foam::identity(pointToGlobal);
            Foam::identity(uniqueMeshPointLabels);
        }
        else
        {
            // Mark up with -1 for unused entries
            pointToGlobal = -1;

            nPoints = 0;

            // Pass 1: markup used points from cells

            for (const label celli : this->cellIds())
            {
                for (const label facei : meshCells[celli])
                {
                    for (const label pointi : mesh.faces()[facei])
                    {
                        if (pointToGlobal[pointi] == -1)
                        {
                            pointToGlobal[pointi] = nPoints++;
                        }
                    }
                }
            }

            // Pass 2:
            //
            // Compact point numbering, preserving original point order
            uniqueMeshPointLabels.resize(nPoints);

            nPoints = 0;
            forAll(pointToGlobal, pointi)
            {
                if (pointToGlobal[pointi] != -1)
                {
                    pointToGlobal[pointi] = nPoints;

                    uniqueMeshPointLabels[nPoints] = pointi;

                    ++nPoints;
                }
            }
        }
    }

    if (notNull(pointToGlobalRequest))
    {
        pointToGlobalRequest.transfer(pointToGlobal);
    }

    return nPoints;
}


Foam::label Foam::ensightCells::uniqueMeshPoints
(
    const polyMesh& mesh,
    labelList& uniqueMeshPointLabels,
    bool parallel
) const
{
    return meshPointMapppings
    (
        mesh,
        const_cast<labelList&>(labelList::null()),  // Ignore pointToGlobal
        uniqueMeshPointLabels,
        parallel
    );
}


// ************************************************************************* //
