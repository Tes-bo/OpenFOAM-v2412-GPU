/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022,2024 OpenCFD Ltd.
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

#include "regionSplit2D.H"
#include "polyMesh.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static constexpr Foam::label UNASSIGNED = -1;
static constexpr Foam::label BLOCKED = -2;


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::regionSplit2D::regionSplit2D
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const boolList& blockedFaces,
    const label offset
)
:
    labelList(patch.size(), UNASSIGNED),
    nRegions_(0)
{
    globalIndex globalFaces(blockedFaces.size());
    label regioni = globalFaces.toGlobal(0);
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());
    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedRegions;
    label nBlockedFaces = 0;
    forAll(blockedFaces, facei)
    {
        if (blockedFaces.test(facei))
        {
            for (const label edgei : patch.faceEdges()[facei])
            {
                changedEdges.append(edgei);

                // Append globally unique value
                changedRegions.append(regioni);
            }
            nBlockedFaces++;
            regioni++;
        }
        else
        {
            // Block all non-seeded faces from the walk
            allFaceInfo[facei] = BLOCKED;
        }
    }

    // Early exit if there are no blocked faces
    if (!returnReduceOr(nBlockedFaces))
    {
        return;
    }

    PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceRegion
    >
    (
        mesh,
        patch,
        changedEdges,
        changedRegions,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(patch.nEdges(), sumOp<label>())
    );


    // Map from regions to local compact indexing
    // - only for regions that originate from this processor
    Map<label> regionToCompactAddr(changedRegions.size());
    label compactRegioni = 0;
    forAll(allFaceInfo, facei)
    {
        const label regioni = allFaceInfo[facei].region();
        if
        (
            globalFaces.isLocal(regioni)
         && regionToCompactAddr.insert(regioni, compactRegioni)
        )
        {
            compactRegioni++;
        }
    }

    // In-place renumber the local regionI to global (compact) regioni
    {
        const label myProcOffset = globalIndex::calcOffset(compactRegioni);
        forAllIters(regionToCompactAddr, iter)
        {
            iter.val() += myProcOffset;
        }
    }


    // Ensure regionToCompactAddr consistent across all processors
    // - not concerned about the op (keys are unique)
    // - map size will be the number of regions in the set of faces
    Pstream::mapCombineReduce(regionToCompactAddr, minEqOp<label>());

    nRegions_ = regionToCompactAddr.size();

    // Set the region index per face
    forAll(allFaceInfo, facei)
    {
        const label regioni = allFaceInfo[facei].region();
        if (regioni >= 0)
        {
            this->operator[](facei) = regionToCompactAddr[regioni] + offset;
        }
    }
}


// ************************************************************************* //
