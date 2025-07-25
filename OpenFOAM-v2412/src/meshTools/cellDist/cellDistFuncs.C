/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020,2024 OpenCFD Ltd.
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

#include "cellDistFuncs.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"
#include "uindirectPrimitivePatch.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellDistFuncs, 0);
}

bool Foam::cellDistFuncs::useCombinedWallPatch = true;

registerOptSwitch
(
    "useCombinedWallPatch",
    bool,
    Foam::cellDistFuncs::useCombinedWallPatch
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellDistFuncs::cellDistFuncs(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelHashSet Foam::cellDistFuncs::getPatchIDs
(
    const UList<wordRe>& patchNames
) const
{
    return mesh().boundaryMesh().patchSet(patchNames, false);
}


// size of largest patch (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::maxPatchSize
(
    const labelHashSet& patchIDs
) const
{
    label maxSize = 0;

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            maxSize = Foam::max(maxSize, patch.size());
        }
    }
    return maxSize;
}


// sum of patch sizes (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::sumPatchSize
(
    const labelHashSet& patchIDs
)
const
{
    label sum = 0;

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            sum += patch.size();
        }
    }
    return sum;
}


// Gets nearest wall for cells next to wall
void Foam::cellDistFuncs::correctBoundaryFaceCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    const auto& pbm = mesh().boundaryMesh();

    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(patchIDs));

    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh().cellCentres();
    const labelList& faceOwner = mesh().faceOwner();

    forAll(pbm, patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = pbm[patchi];
            const auto areaFraction(patch.areaFraction());

            // Check cells with face on wall
            forAll(patch, patchFacei)
            {
                if (areaFraction && (areaFraction()[patchFacei] <= 0.5))
                {
                    // For cyclicACMI: more cyclic than wall
                    continue;
                }

                getPointNeighbours(patch, patchFacei, neighbours);

                label celli = faceOwner[patch.start() + patchFacei];

                label minFacei = -1;

                wallDistCorrected[celli] = smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minFacei
                );

                // Store wallCell and its nearest neighbour
                nearestFace.insert(celli, patch.start()+minFacei);
            }
        }
    }
}


// Correct all cells connected to wall (via point) and not in nearestFace
void Foam::cellDistFuncs::correctBoundaryPointCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    // Correct all (non-visited) cells with point on wall

    const auto& pbm = mesh().boundaryMesh();
    const vectorField& cellCentres = mesh().cellCentres();

    forAll(pbm, patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = pbm[patchi];
            const auto& localFaces = patch.localFaces();
            const labelList& meshPoints = patch.meshPoints();
            const labelListList& pointFaces = patch.pointFaces();

            bitSet isWallPoint(meshPoints.size(), true);
            {
                const auto areaFraction(patch.areaFraction());

                // Check cells with face on wall
                forAll(patch, patchFacei)
                {
                    if (areaFraction && (areaFraction()[patchFacei] <= 0.5))
                    {
                        // For cyclicACMI: more cyclic than wall
                        isWallPoint.unset(localFaces[patchFacei]);
                    }
                }
            }


            forAll(meshPoints, patchPointi)
            {
                const label vertI = meshPoints[patchPointi];

                if (!isWallPoint[patchPointi])
                {
                    continue;
                }

                const labelList& neighbours = mesh().pointCells(vertI);

                for (const label celli : neighbours)
                {
                    if (!nearestFace.found(celli))
                    {
                        const labelList& wallFaces = pointFaces[patchPointi];

                        label minFacei = -1;

                        wallDistCorrected[celli] = smallestDist
                        (
                            cellCentres[celli],
                            patch,
                            wallFaces,
                            minFacei
                        );

                        // Store wallCell and its nearest neighbour
                        nearestFace.insert(celli, patch.start()+minFacei);
                    }
                }
            }
        }
    }
}


void Foam::cellDistFuncs::correctBoundaryCells
(
    const labelList& patchIDs,
    const bool doPointCells,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    label nWalls = 0;
    {
        for (const label patchi : patchIDs)
        {
            nWalls += mesh().boundaryMesh()[patchi].size();
        }
    }


    DynamicList<label> faceLabels(nWalls);
    {
        for (const label patchi : patchIDs)
        {
            const auto& patch = mesh().boundaryMesh()[patchi];
            forAll(patch, i)
            {
                faceLabels.append(patch.start()+i);
            }
        }
    }

    const uindirectPrimitivePatch wallPatch
    (
        UIndirectList<face>(mesh().faces(), faceLabels),
        mesh_.points()
    );


    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh().cellCentres();

    DynamicList<label> neighbours;

    nWalls = 0;
    for (const label patchi : patchIDs)
    {
        const auto& patch = mesh().boundaryMesh()[patchi];
        const auto areaFraction(patch.areaFraction());
        const labelUList& faceCells = patch.faceCells();

        // Check cells with face on wall
        forAll(patch, patchFacei)
        {
            if (areaFraction && (areaFraction()[patchFacei] <= 0.5))
            {
                // For cyclicACMI: more cyclic than wall
            }
            else
            {
                getPointNeighbours(wallPatch, nWalls, neighbours);

                const label celli = faceCells[patchFacei];

                label minFacei = -1;
                wallDistCorrected[celli] = smallestDist
                (
                    cellCentres[celli],
                    wallPatch,
                    neighbours,
                    minFacei
                );

                // Store wallCell and its nearest neighbour
                nearestFace.insert(celli, nWalls+minFacei);
            }

            nWalls++;
        }
    }

    // Correct all cells with a point on the wall
    if (doPointCells)
    {
        const auto& meshPoints = wallPatch.meshPoints();
        const auto& localFaces = wallPatch.localFaces();

        bitSet isWallPoint(meshPoints.size(), true);

        nWalls = 0;
        for (const label patchi : patchIDs)
        {
            const auto& patch = mesh().boundaryMesh()[patchi];
            const auto areaFraction(patch.areaFraction());

            // Check cells with face on wall
            forAll(patch, patchFacei)
            {
                if (areaFraction && (areaFraction()[patchFacei] <= 0.5))
                {
                    // For cyclicACMI: more cyclic than wall
                    isWallPoint.unset(localFaces[nWalls]);
                }

                nWalls++;
            }
        }

        const auto& pointFaces = wallPatch.pointFaces();

        for (const label patchPointi : isWallPoint)
        {
            const label verti = meshPoints[patchPointi];

            const labelList& neighbours = mesh().pointCells(verti);

            for (const label celli : neighbours)
            {
                if (!nearestFace.found(celli))
                {
                    const labelList& wallFaces = pointFaces[patchPointi];

                    label minFacei = -1;

                    wallDistCorrected[celli] = smallestDist
                    (
                        cellCentres[celli],
                        wallPatch,
                        wallFaces,
                        minFacei
                    );

                    // Store wallCell and its nearest neighbour
                    nearestFace.insert(celli, wallPatch.addressing()[minFacei]);
                }
            }
        }
    }
}


// ************************************************************************* //
