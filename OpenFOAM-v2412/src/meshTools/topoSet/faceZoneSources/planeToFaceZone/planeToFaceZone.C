/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "planeToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(planeToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, planeToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, planeToFaceZone, istream);

    addToRunTimeSelectionTable(topoSetFaceZoneSource, planeToFaceZone, word);
    addToRunTimeSelectionTable(topoSetFaceZoneSource, planeToFaceZone, istream);

    addNamedToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        planeToFaceZone,
        word,
        plane
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        planeToFaceZone,
        istream,
        plane
    );
}


Foam::topoSetSource::addToUsageTable Foam::planeToFaceZone::usage_
(
    planeToFaceZone::typeName,
    "\n    Usage: planeToFaceZone (px py pz) (nx ny nz) include\n\n"
    "    Select faces for which the adjacent cell centres lie on opposite "
    " of a plane\n\n"
);


const Foam::Enum
<
    Foam::planeToFaceZone::faceAction
>
Foam::planeToFaceZone::faceActionNames_
({
    { faceAction::ALL, "all" },
    { faceAction::CLOSEST, "closest" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::planeToFaceZone::combine(faceZoneSet& fzSet, const bool add) const
{
    // Mark all cells with centres above the plane
    bitSet cellIsAbovePlane(mesh_.nCells());
    forAll(mesh_.cells(), celli)
    {
        if (((mesh_.cellCentres()[celli] - point_) & normal_) > 0)
        {
            cellIsAbovePlane.set(celli);
        }
    }

    // Mark all faces that sit between cells above and below the plane
    bitSet faceIsOnPlane(mesh_.nFaces());
    forAll(mesh_.faceNeighbour(), facei)
    {
        if
        (
            cellIsAbovePlane[mesh_.faceOwner()[facei]]
         != cellIsAbovePlane[mesh_.faceNeighbour()[facei]]
        )
        {
            faceIsOnPlane.set(facei);
        }
    }
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        forAll(patch, patchFacei)
        {
            const label facei = patch.start() + patchFacei;
            if (patch.coupled() && cellIsAbovePlane[mesh_.faceOwner()[facei]])
            {
                faceIsOnPlane.set(facei);
            }
        }
    }
    syncTools::syncFaceList(mesh_, faceIsOnPlane, xorEqOp<unsigned int>());

    // Convert marked faces to a list of indices
    labelList newSetFaces(faceIsOnPlane.sortedToc());
    faceIsOnPlane.clear();

    // If constructing a single contiguous set, remove all faces except those
    // connected to the contiguous region closest to the specified point
    if (option_ == faceAction::CLOSEST)
    {
        // Step 1: Get locally contiguous regions for the new face set and the
        // total number of regions across all processors.
        labelList newSetFaceRegions(newSetFaces.size(), -1);
        label nRegions = -1;
        {
            // Create a patch of the set faces
            const uindirectPrimitivePatch newSetPatch
            (
                UIndirectList<face>(mesh_.faces(), newSetFaces),
                mesh_.points()
            );

            // Get the region ID-s and store the total number of regions on
            // each processor
            labelList procNRegions(Pstream::nProcs(), -1);
            procNRegions[Pstream::myProcNo()] =
                PatchTools::markZones
                (
                    newSetPatch,
                    bitSet(),  // No border edges
                    newSetFaceRegions
                );
            Pstream::allGatherList(procNRegions);

            // Cumulative sum the number of regions on each processor to get an
            // offset which makes the local region ID-s globally unique
            labelList procRegionOffset(Pstream::nProcs(), Zero);
            for (label proci = 1; proci < Pstream::nProcs(); ++proci)
            {
                procRegionOffset[proci] =
                    procRegionOffset[proci - 1]
                  + procNRegions[proci - 1];
            }

            // Apply the offset
            for (label& regioni : newSetFaceRegions)
            {
                regioni += procRegionOffset[Pstream::myProcNo()];
            }

            // Store the total number of regions across all processors
            nRegions = procRegionOffset.last() + procNRegions.last();
        }

        // Step 2: Create a region map which combines regions which are
        // connected across coupled interfaces
        labelList regionMap(identity(nRegions));
        {
            // Put region labels on connected boundary edges and synchronise to
            // create a list of all regions connected to a given edge
            labelListList meshEdgeRegions(mesh_.nEdges(), labelList());
            forAll(newSetFaces, newSetFacei)
            {
                const label facei = newSetFaces[newSetFacei];
                const label regioni = newSetFaceRegions[newSetFacei];
                for (const label edgei : mesh_.faceEdges()[facei])
                {
                    meshEdgeRegions[edgei] = labelList(one{}, regioni);
                }
            }
            syncTools::syncEdgeList
            (
                mesh_,
                meshEdgeRegions,
                ListOps::appendEqOp<label>(),
                labelList()
            );

            // Combine edge regions to create a list of what regions a given
            // region is connected to
            List<bitSet> regionRegions(nRegions);
            forAll(newSetFaces, newSetFacei)
            {
                const label facei = newSetFaces[newSetFacei];
                const label regioni = newSetFaceRegions[newSetFacei];
                for (const label edgei : mesh_.faceEdges()[facei])
                {
                    // Includes self region (removed below)
                    regionRegions[regioni].set(meshEdgeRegions[edgei]);
                }
                forAll(regionRegions, regioni)
                {
                    // Remove self region
                    regionRegions[regioni].unset(regioni);
                }
            }
            Pstream::listCombineReduce(regionRegions, bitOrEqOp<bitSet>());

            // Collapse the region connections into a map between each region
            // and the lowest numbered region that it connects to
            forAll(regionRegions, regioni)
            {
                for (const label regi : regionRegions[regioni])
                {
                    // minEqOp<label>()
                    regionMap[regi] = min(regionMap[regi], regionMap[regioni]);
                }
            }
        }

        // Step 3: Combine connected regions
        labelList regionNFaces;
        {
            // Remove duplicates from the region map
            label regioni0 = 0;
            forAll(regionMap, regioni)
            {
                if (regionMap[regioni] > regioni0)
                {
                    ++ regioni0;
                    regionMap[regioni] = regioni0;
                }
            }

            // Recompute the number of regions
            nRegions = regioni0 + 1;

            // Renumber the face region ID-s
            newSetFaceRegions =
                IndirectList<label>(regionMap, newSetFaceRegions);

            // Report the final number and size of the regions
            regionNFaces = labelList(nRegions, Zero);
            for (const label regioni : newSetFaceRegions)
            {
                ++ regionNFaces[regioni];
            }
            Pstream::listCombineReduce(regionNFaces, plusEqOp<label>());

            Info<< "    Found " << nRegions << " contiguous regions with "
                << regionNFaces << " faces" << endl;
        }

        // Step 4: Choose the closest region to output
        label selectedRegioni = -1;
        {
            // Compute the region centres
            scalarField regionWeights(nRegions, Zero);
            pointField regionCentres(nRegions, Zero);
            forAll(newSetFaces, newSetFacei)
            {
                const label facei = newSetFaces[newSetFacei];
                const label regioni = newSetFaceRegions[newSetFacei];

                const scalar w = mag(mesh_.faceAreas()[facei]);
                const point& c = mesh_.faceCentres()[facei];

                regionWeights[regioni] += w;
                regionCentres[regioni] += w*c;
            }
            Pstream::listCombineGather(regionWeights, plusEqOp<scalar>());
            Pstream::listCombineGather(regionCentres, plusEqOp<point>());

            if (Pstream::master())
            {
                regionCentres /= regionWeights;
            }
            //Pstream::broadcast(regionWeights);
            Pstream::broadcast(regionCentres);

            // Find the region centroid closest to the reference point
            selectedRegioni =
                returnReduce
                (
                    findMin(mag(regionCentres - point_)()),
                    minOp<label>()
                );

            // Report the selection
            Info<< "    Selecting region " << selectedRegioni << " with "
                << regionNFaces[selectedRegioni]
                << " faces as the closest to point " << point_ << endl;
        }

        // Step 5: Remove any faces from the set list not in the selected region
        {
            // Remove faces from the list by shuffling up and resizing
            label newSetFacei0 = 0;
            forAll(newSetFaces, newSetFacei)
            {
                newSetFaces[newSetFacei0] = newSetFaces[newSetFacei];

                if (newSetFaceRegions[newSetFacei] == selectedRegioni)
                {
                    ++ newSetFacei0;
                }
            }
            newSetFaces.resize(newSetFacei0);
        }
    }

    // Modify the face zone set
    DynamicList<label> newAddressing;
    DynamicList<bool> newFlipMap;
    if (add)
    {
        // Start from copy
        newAddressing = fzSet.addressing();
        newFlipMap = fzSet.flipMap();

        // Add anything from the new set that is not already in the zone set
        const auto& exclude = fzSet;
        for (const label facei : newSetFaces)
        {
            if (!exclude.found(facei))
            {
                newAddressing.append(facei);
                newFlipMap.append(cellIsAbovePlane[mesh_.faceOwner()[facei]]);
            }
        }
    }
    else
    {
        // Start from empty
        newAddressing.reserve(fzSet.addressing().size());
        newFlipMap.reserve(newAddressing.capacity());

        // Add everything from the zone set that is not also in the new set
        bitSet exclude(newSetFaces);
        for (const label facei : fzSet.addressing())
        {
            if (!exclude.found(facei))
            {
                newAddressing.append(facei);
                newFlipMap.append(cellIsAbovePlane[mesh_.faceOwner()[facei]]);
            }
        }
    }
    fzSet.addressing().transfer(newAddressing);
    fzSet.flipMap().transfer(newFlipMap);
    fzSet.updateSet();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::planeToFaceZone::planeToFaceZone
(
    const polyMesh& mesh,
    const point& basePoint,
    const vector& normal,
    const faceAction action
)
:
    topoSetFaceZoneSource(mesh),
    point_(basePoint),
    normal_(normal),
    option_(action)
{}


Foam::planeToFaceZone::planeToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    point_(dict.get<vector>("point")),
    normal_(dict.get<vector>("normal")),
    option_(faceActionNames_.getOrDefault("option", dict, faceAction::ALL))
{}


Foam::planeToFaceZone::planeToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    point_(checkIs(is)),
    normal_(checkIs(is)),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::planeToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
        return;
    }

    faceZoneSet& zoneSet = refCast<faceZoneSet>(set);

    if (action == topoSetSource::NEW || action == topoSetSource::ADD)
    {
        if (verbose_)
        {
            Info<< "    Adding faces that form a plane at "
                << point_ << " with normal " << normal_ << endl;
        }

        combine(zoneSet, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces that form a plane at "
                << point_ << " with normal " << normal_ << endl;
        }

        combine(zoneSet, false);
    }
}


// ************************************************************************* //
