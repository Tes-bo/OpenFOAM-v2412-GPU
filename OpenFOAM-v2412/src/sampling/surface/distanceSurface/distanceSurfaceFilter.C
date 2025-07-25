/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "distanceSurface.H"
#include "regionSplit.H"
#include "syncTools.H"
#include "ListOps.H"
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::distanceSurface::refineBlockedCells
(
    bitSet& ignoreCells,
    const isoSurfaceBase& isoCutter
) const
{
    // With the cell/point distance fields we can take a second pass at
    // pre-filtering.
    // This duplicates how cut detection is determined in the cell/topo
    // algorithms but is fairly inexpensive (creates no geometry)

    bool changed = false;

    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        if (ignoreCells.test(celli))
        {
            continue;
        }

        auto cut = isoCutter.getCellCutType(celli);
        if (!(cut & isoSurfaceBase::ANYCUT))
        {
            ignoreCells.set(celli);
            changed = true;
        }
    }

    return changed;
}


Foam::bitSet Foam::distanceSurface::filterPrepareRegionSplit
(
    const bitSet& ignoreCells
) const
{
    // Prepare for region split

    bitSet blockedFaces(mesh_.nFaces());

    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();

    // Could be more efficient
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        // If only one cell is blocked, the face corresponds
        // to an exposed subMesh face

        if
        (
            ignoreCells.test(faceOwn[facei])
         != ignoreCells.test(faceNei[facei])
        )
        {
            blockedFaces.set(facei);
        }
    }

    for (const polyPatch& patch : mesh_.boundaryMesh())
    {
        if (!patch.coupled())
        {
            continue;
        }

        forAll(patch, patchFacei)
        {
            const label facei = patch.start() + patchFacei;
            if (ignoreCells.test(faceOwn[facei]))
            {
                blockedFaces.set(facei);
            }
        }
    }

    syncTools::syncFaceList(mesh_, blockedFaces, xorEqOp<unsigned int>());

    return blockedFaces;
}


void Foam::distanceSurface::filterKeepLargestRegion
(
    bitSet& ignoreCells
) const
{
    // For region split
    bitSet blockedFaces(filterPrepareRegionSplit(ignoreCells));

    // Split region
    regionSplit rs(mesh_, blockedFaces);
    blockedFaces.clearStorage();

    const labelList& regionColour = rs;

    // Identical number of regions on all processors
    labelList nCutsPerRegion(rs.nRegions(), Zero);

    // Count cut cells (ie, unblocked)
    forAll(regionColour, celli)
    {
        if (!ignoreCells.test(celli))
        {
            ++nCutsPerRegion[regionColour[celli]];
        }
    }

    // Sum totals from all processors
    Pstream::listCombineGather(nCutsPerRegion, plusEqOp<label>());


    // Define which regions to keep
    boolList keepRegion(rs.nRegions(), false);

    if (Pstream::master())
    {
        const label largest = findMax(nCutsPerRegion);

        if (largest == -1)
        {
            // Should not happen
            keepRegion = true;
        }
        else
        {
            keepRegion[largest] = true;
        }

        if (debug)
        {
            Info<< "Had " << sum(nCutsPerRegion) << " cuts, in "
                << nCutsPerRegion.size() << " regions, largest is "
                << largest <<  ": " << flatOutput(nCutsPerRegion) << nl;
        }
    }

    Pstream::broadcast(keepRegion);

    forAll(regionColour, celli)
    {
        if (!keepRegion.test(regionColour[celli]))
        {
            ignoreCells.set(celli);
        }
    }
}


void Foam::distanceSurface::filterKeepNearestRegions
(
    bitSet& ignoreCells
) const
{
    if (nearestPoints_.empty())
    {
        WarningInFunction
            << "Ignoring nearestPoints - no points provided" << nl
            << endl;
        return;
    }

    // For region split
    bitSet blockedFaces(filterPrepareRegionSplit(ignoreCells));

    // Split region
    regionSplit rs(mesh_, blockedFaces);
    blockedFaces.clearStorage();

    const labelList& regionColour = rs;

    const pointField& cc = mesh_.cellCentres();
    const pointField& nearPts = nearestPoints_;

    // The magSqr distance and region
    typedef Tuple2<scalar, label> nearInfo;
    List<nearInfo> nearest(nearPts.size(), nearInfo(GREAT, -1));

    // Also collect cuts per region, may be useful for rejecting
    // small regions. Code as per filterKeepLargestRegion
    labelList nCutsPerRegion(rs.nRegions(), Zero);

    forAll(cc, celli)
    {
        if (ignoreCells.test(celli))
        {
            continue;
        }

        const point& pt = cc[celli];
        const label regioni = regionColour[celli];

        ++nCutsPerRegion[regioni];

        label pointi = 0;
        for (nearInfo& near : nearest)
        {
            const scalar distSqr = magSqr(nearPts[pointi] - pt);
            ++pointi;

            if (distSqr < near.first())
            {
                near.first() = distSqr;
                near.second() = regioni;
            }
        }
    }

    // Sum totals from all processors
    Pstream::listCombineGather(nCutsPerRegion, plusEqOp<label>());

    // Get nearest
    Pstream::listCombineGather(nearest, minFirstEqOp<scalar>());


    // Define which regions to keep

    boolList keepRegion(rs.nRegions(), false);

    if (Pstream::master())
    {
        const label largest = findMax(nCutsPerRegion);

        for (const nearInfo& near : nearest)
        {
            const scalar distSqr = near.first();
            const label regioni = near.second();

            if (regioni != -1 && distSqr < maxDistanceSqr_)
            {
                keepRegion[regioni] = true;
            }
        }

        if (debug)
        {
            Info<< "Had " << sum(nCutsPerRegion) << " cuts, in "
                << nCutsPerRegion.size() << " regions, largest is "
                << largest <<  ": " << flatOutput(nCutsPerRegion) << nl;

            Info<< "nearestPoints (max distance = "
                << sqrt(maxDistanceSqr_) << ")" << nl;

            forAll(nearPts, pointi)
            {
                const scalar dist = sqrt(nearest[pointi].first());
                const label regioni = nearest[pointi].second();

                Info<< "    " << nearPts[pointi] << " region "
                    << regioni << " distance "
                    << dist;

                if (!keepRegion.test(regioni))
                {
                    Info<< " too far";
                }
                Info<< nl;
            }
        }
    }

    Pstream::broadcast(keepRegion);

    forAll(regionColour, celli)
    {
        if (!keepRegion.test(regionColour[celli]))
        {
            ignoreCells.set(celli);
        }
    }
}


void Foam::distanceSurface::filterRegionProximity
(
    bitSet& ignoreCells
) const
{
    const searchableSurface& geom = geometryPtr_();

    // For face distances
    const pointField& fc = surface_.faceCentres();

    // For region split
    bitSet blockedFaces(filterPrepareRegionSplit(ignoreCells));

    // Split region
    regionSplit rs(mesh_, blockedFaces);
    blockedFaces.clearStorage();

    const labelList& regionColour = rs;

    // For each face
    scalarField faceDistance(fc.size(), GREAT);

    {
        List<pointIndexHit> nearest;
        geom.findNearest
        (
            fc,
            // Use initialized field (GREAT) to limit search too
            faceDistance,
            nearest
        );
        calcAbsoluteDistance(faceDistance, fc, nearest);
    }

    // Identical number of regions on all processors
    scalarField areaRegion(rs.nRegions(), Zero);
    scalarField distRegion(rs.nRegions(), Zero);

    forAll(meshCells_, facei)
    {
        const label celli = meshCells_[facei];
        const label regioni = regionColour[celli];

        const scalar faceArea = surface_[facei].mag(surface_.points());
        distRegion[regioni] += (faceDistance[facei] * faceArea);
        areaRegion[regioni] += (faceArea);
    }

    Pstream::listCombineGather(distRegion, plusEqOp<scalar>());
    Pstream::listCombineGather(areaRegion, plusEqOp<scalar>());

    if (Pstream::master())
    {
        forAll(distRegion, regioni)
        {
            distRegion[regioni] /= (areaRegion[regioni] + VSMALL);
        }
    }

    Pstream::broadcast(distRegion);


    // Define the per-face acceptance based on the region average distance

    bitSet acceptFaces(fc.size());
    bool prune(false);

    forAll(meshCells_, facei)
    {
        const label celli = meshCells_[facei];
        const label regioni = regionColour[celli];

        // NB: do not filter by individual faces as well since this
        // has been reported to cause minor holes for surfaces with
        // high curvature! (2021-06-10)

        if (absProximity_ < distRegion[regioni])
        {
            prune = true;
        }
        else
        {
            acceptFaces.set(facei);
        }
    }

    // Heavier debugging
    if (debug & 4)
    {
        const fileName outputName(surfaceName() + "-region-proximity-filter");

        Info<< "Writing debug surface: " << outputName << nl;

        surfaceWriters::vtkWriter writer
        (
            surface_.points(),
            surface_,  // faces
            outputName
        );

        writer.write("absolute-distance", faceDistance);

        // Region segmentation
        labelField faceRegion
        (
            ListOps::create<label>
            (
                meshCells_,
                [&](const label celli){ return regionColour[celli]; }
            )
        );
        writer.write("face-region", faceRegion);

        // Region-wise filter state
        labelField faceFilterState
        (
            ListOps::createWithValue<label>(surface_.size(), acceptFaces, 1, 0)
        );
        writer.write("filter-state", faceFilterState);
    }


    // If filtering with region proximity results in zero faces,
    // revert to face-only proximity filter
    if (returnReduceAnd(acceptFaces.none()))
    {
        acceptFaces.reset();
        prune = false;

        // Consider the absolute proximity of the face centres
        forAll(faceDistance, facei)
        {
            if (absProximity_ < faceDistance[facei])
            {
                prune = true;
            }
            else
            {
                acceptFaces.set(facei);
            }
        }
    }

    if (prune)
    {
        labelList pointMap, faceMap;
        meshedSurface filtered
        (
            surface_.subsetMesh(acceptFaces, pointMap, faceMap)
        );
        surface_.transfer(filtered);

        meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
    }
}


void Foam::distanceSurface::filterFaceProximity()
{
    const searchableSurface& geom = geometryPtr_();

    // For face distances
    const pointField& fc = surface_.faceCentres();

    // For each face
    scalarField faceDistance(fc.size(), GREAT);
    scalarField faceNormalDistance;  // Debugging only
    {
        List<pointIndexHit> nearest;
        geom.findNearest
        (
            fc,
            // Use initialized field (GREAT) to limit search too
            faceDistance,
            nearest
        );
        calcAbsoluteDistance(faceDistance, fc, nearest);

        // Heavier debugging
        if (debug & 4)
        {
            vectorField norms;
            geom.getNormal(nearest, norms);

            faceNormalDistance.resize(fc.size());

            forAll(nearest, i)
            {
                const vector diff(fc[i] - nearest[i].point());

                faceNormalDistance[i] = Foam::mag((diff & norms[i]));
            }
        }
    }

    // Note
    // Proximity checks using the face points (nearest hit) to establish
    // a length scale are too fragile. Can easily have stretched faces
    // where the centre is less than say 0.3-0.5 of the centre-point distance
    // but they are still outside.

    // Using the absolute proximity of the face centres is more robust.

    bitSet acceptFaces(fc.size());
    bool prune(false);

    // Consider the absolute proximity of the face centres
    forAll(faceDistance, facei)
    {
        if (absProximity_ < faceDistance[facei])
        {
            prune = true;
            if (debug & 2)
            {
                Pout<< "trim reject: "
                    << faceDistance[facei] << nl;
            }
        }
        else
        {
            acceptFaces.set(facei);
        }
    }


    // Heavier debugging
    if (debug & 4)
    {
        const fileName outputName(surfaceName() + "-face-proximity-filter");

        Info<< "Writing debug surface: " << outputName << nl;

        surfaceWriters::vtkWriter writer
        (
            surface_.points(),
            surface_,  // faces
            outputName
        );

        writer.write("absolute-distance", faceDistance);
        writer.write("normal-distance", faceNormalDistance);

        // Region-wise filter state
        labelField faceFilterState
        (
            ListOps::createWithValue<label>(surface_.size(), acceptFaces, 1, 0)
        );
        writer.write("filter-state", faceFilterState);
    }


    if (prune)
    {
        labelList pointMap, faceMap;
        meshedSurface filtered
        (
            surface_.subsetMesh(acceptFaces, pointMap, faceMap)
        );
        surface_.transfer(filtered);

        meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
    }
}


// ************************************************************************* //
