/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "bitSet.H"
#include "volFields.H"

// Leak path
#include "shortestPathSet.H"
#include "meshSearch.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"
#include "removeCells.H"
#include "regionSplit.H"

#include "volFields.H"
#include "wallPoints.H"
#include "searchableSurfaces.H"
#include "distributedTriSurfaceMesh.H"

#include "holeToFace.H"
#include "refinementParameters.H"
#include "indirectPrimitivePatch.H"
#include "OBJstream.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::label Foam::meshRefinement::markFakeGapRefinement
//(
//    const scalar planarCos,
//
//    const label nAllowRefine,
//    const labelList& neiLevel,
//    const pointField& neiCc,
//
//    labelList& refineCell,
//    label& nRefine
//) const
//{
//    label oldNRefine = nRefine;
//
//
//    // Collect candidate faces (i.e. intersecting any surface and
//    // owner/neighbour not yet refined.
//    const labelList testFaces(getRefineCandidateFaces(refineCell));
//
//    // Collect segments
//    pointField start(testFaces.size());
//    pointField end(testFaces.size());
//    labelList minLevel(testFaces.size());
//
//    calcCellCellRays
//    (
//        neiCc,
//        neiLevel,
//        testFaces,
//        start,
//        end,
//        minLevel
//    );
//
//
//    // Re-use the gap shooting methods. This needs:
//    //  - shell gapLevel : faked. Set to 0,labelMax
//    //  - surface gapLevel : faked by overwriting
//
//
//    List<FixedList<label, 3>>& surfGapLevel = const_cast
//    <
//        List<FixedList<label, 3>>&
//    >(surfaces_.extendedGapLevel());
//
//    List<volumeType>& surfGapMode = const_cast
//    <
//        List<volumeType>&
//    >(surfaces_.extendedGapMode());
//
//    const List<FixedList<label, 3>> surfOldLevel(surfGapLevel);
//    const List<volumeType> surfOldMode(surfGapMode);
//
//    // Set the extended gap levels
//    forAll(surfaces_.gapLevel(), regioni)
//    {
//        surfGapLevel[regioni] = FixedList<label, 3>
//        ({
//            3,
//           -1,
//            surfaces_.gapLevel()[regioni]+1
//        });
//    }
//    surfGapMode = volumeType::MIXED;
//
//Pout<< "gapLevel was:" << surfOldLevel << endl;
//Pout<< "gapLevel now:" << surfGapLevel << endl;
//Pout<< "gapMode was:" << surfOldMode << endl;
//Pout<< "gapMode now:" << surfGapMode << endl;
//Pout<< "nRefine was:" << oldNRefine << endl;
//
//
//
//    List<List<FixedList<label, 3>>>& shellGapLevel = const_cast
//    <
//        List<List<FixedList<label, 3>>>&
//    >(shells_.extendedGapLevel());
//
//    List<List<volumeType>>& shellGapMode = const_cast
//    <
//        List<List<volumeType>>&
//    >(shells_.extendedGapMode());
//
//    const List<List<FixedList<label, 3>>> shellOldLevel(shellGapLevel);
//    const List<List<volumeType>> shellOldMode(shellGapMode);
//
//    // Set the extended gap levels
//    forAll(shellGapLevel, shelli)
//    {
//        shellGapLevel[shelli] =  FixedList<label, 3>({3, -1, labelMax});
//        shellGapMode[shelli] = volumeType::MIXED;
//    }
//Pout<< "shellLevel was:" << shellOldLevel << endl;
//Pout<< "shellLevel now:" << shellGapLevel << endl;
//
//    const label nAdditionalRefined = markSurfaceGapRefinement
//    (
//        planarCos,
//
//        nAllowRefine,
//        neiLevel,
//        neiCc,
//
//        refineCell,
//        nRefine
//    );
//
//Pout<< "nRefine now:" << nRefine << endl;
//
//    // Restore
//    surfGapLevel = surfOldLevel;
//    surfGapMode = surfOldMode;
//    shellGapLevel = shellOldLevel;
//    shellGapMode = shellOldMode;
//
//    return nAdditionalRefined;
//}


void Foam::meshRefinement::markOutsideFaces
(
    const labelList& cellLevel,
    const labelList& neiLevel,
    const labelList& refineCell,
    bitSet& isOutsideFace
) const
{
    // Get faces:
    // - on outside of cell set
    // - inbetween same cell level (i.e. quads)

    isOutsideFace.setSize(mesh_.nFaces());
    isOutsideFace = Zero;

    forAll(mesh_.faceNeighbour(), facei)
    {
        label own = mesh_.faceOwner()[facei];
        label nei = mesh_.faceNeighbour()[facei];
        if
        (
            (cellLevel[own] == cellLevel[nei])
         && (
                (refineCell[own] != -1)
             != (refineCell[nei] != -1)
            )
        )
        {
            isOutsideFace.set(facei);
        }
    }
    {

        const label nBnd = mesh_.nBoundaryFaces();

        labelList neiRefineCell(nBnd);
        syncTools::swapBoundaryCellList(mesh_, refineCell, neiRefineCell);
        for (label bFacei = 0; bFacei < nBnd; ++bFacei)
        {
            label facei = mesh_.nInternalFaces()+bFacei;
            label own = mesh_.faceOwner()[facei];

            if
            (
                (cellLevel[own] == neiLevel[bFacei])
             && (
                    (refineCell[own] != -1)
                 != (neiRefineCell[bFacei] != -1)
                )
            )
            {
                isOutsideFace.set(facei);
            }
        }
    }
}


Foam::label Foam::meshRefinement::countFaceDirs
(
    const bitSet& isOutsideFace,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    const vectorField& faceAreas = mesh_.faceAreas();

    Vector<bool> haveDirs(vector::uniform(false));

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        if (isOutsideFace[facei])
        {
            const vector& n = faceAreas[facei];
            scalar magSqrN = magSqr(n);

            if (magSqrN > ROOTVSMALL)
            {
                for
                (
                    direction dir = 0;
                    dir < pTraits<vector>::nComponents;
                    dir++
                )
                {
                    if (Foam::sqr(n[dir]) > 0.99*magSqrN)
                    {
                        haveDirs[dir] = true;
                        break;
                    }
                }
            }
        }
    }

    label nDirs = 0;
    forAll(haveDirs, dir)
    {
        if (haveDirs[dir])
        {
            nDirs++;
        }
    }
    return nDirs;
}


void Foam::meshRefinement::growSet
(
    const labelList& neiLevel,
    const bitSet& isOutsideFace,
    labelList& refineCell,
    label& nRefine
) const
{
    // Get cells with three or more outside faces
    const cellList& cells = mesh_.cells();
    forAll(cells, celli)
    {
        if (refineCell[celli] == -1)
        {
            if (countFaceDirs(isOutsideFace, celli) == 3)
            {
                // Mark cell with any value
                refineCell[celli] = 0;
                nRefine++;
            }
        }
    }
}


//void Foam::meshRefinement::markMultiRegionCell
//(
//    const label celli,
//    const FixedList<label, 3>& surface,
//
//    Map<FixedList<label, 3>>& cellToRegions,
//    bitSet& isMultiRegion
//) const
//{
//    if (!isMultiRegion[celli])
//    {
//        Map<FixedList<label, 3>>::iterator fnd = cellToRegions.find(celli);
//        if (!fnd.good())
//        {
//            cellToRegions.insert(celli, surface);
//        }
//        else if (fnd() != surface)
//        {
//            // Found multiple intersections on cell
//            isMultiRegion.set(celli);
//        }
//    }
//}


//void Foam::meshRefinement::detectMultiRegionCells
//(
//    const labelListList& faceZones,
//    const labelList& testFaces,
//
//    const labelList& surface1,
//    const List<pointIndexHit>& hit1,
//    const labelList& region1,
//
//    const labelList& surface2,
//    const List<pointIndexHit>& hit2,
//    const labelList& region2,
//
//    bitSet& isMultiRegion
//) const
//{
//    isMultiRegion.clear();
//    isMultiRegion.setSize(mesh_.nCells());
//
//    Map<FixedList<label, 3>> cellToRegions(testFaces.size());
//
//    forAll(testFaces, i)
//    {
//        const pointIndexHit& info1 = hit1[i];
//        if (info1.hit())
//        {
//            const label facei = testFaces[i];
//            const labelList& fz1 = faceZones[surface1[i]];
//            const FixedList<label, 3> surfaceInfo1
//            ({
//                surface1[i],
//                region1[i],
//                (fz1.size() ? fz1[info1.index()] : region1[i])
//            });
//
//            markMultiRegionCell
//            (
//                mesh_.faceOwner()[facei],
//                surfaceInfo1,
//                cellToRegions,
//                isMultiRegion
//            );
//            if (mesh_.isInternalFace(facei))
//            {
//                markMultiRegionCell
//                (
//                    mesh_.faceNeighbour()[facei],
//                    surfaceInfo1,
//                    cellToRegions,
//                    isMultiRegion
//                );
//            }
//
//            const pointIndexHit& info2 = hit2[i];
//
//            if (info2.hit() && info1 != info2)
//            {
//                const labelList& fz2 = faceZones[surface2[i]];
//                const FixedList<label, 3> surfaceInfo2
//                ({
//                    surface2[i],
//                    region2[i],
//                    (fz2.size() ? fz2[info2.index()] : region2[i])
//                });
//
//                markMultiRegionCell
//                (
//                    mesh_.faceOwner()[facei],
//                    surfaceInfo2,
//                    cellToRegions,
//                    isMultiRegion
//                );
//                if (mesh_.isInternalFace(facei))
//                {
//                    markMultiRegionCell
//                    (
//                        mesh_.faceNeighbour()[facei],
//                        surfaceInfo2,
//                        cellToRegions,
//                        isMultiRegion
//                    );
//                }
//            }
//        }
//    }
//
//
//    if (debug&meshRefinement::MESH)
//    {
//        volScalarField multiCell
//        (
//            IOobject
//            (
//                "multiCell",
//                mesh_.time().timeName(),
//                mesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE,
//                IOobject::NO_REGISTER
//            ),
//            mesh_,
//            dimensionedScalar
//            (
//                "zero",
//                dimensionSet(0, 1, 0, 0, 0),
//                0.0
//            )
//        );
//        forAll(isMultiRegion, celli)
//        {
//            if (isMultiRegion[celli])
//            {
//                multiCell[celli] = 1.0;
//            }
//        }
//
//        Info<< "Writing all multi cells to " << multiCell.name() << endl;
//        multiCell.write();
//    }
//}


void Foam::meshRefinement::markProximityCandidateFaces
(
    const labelList& blockedSurfaces,
    const List<scalarList>& regionToBlockSize,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& testFaces
) const
{
    labelListList faceZones(surfaces_.surfaces().size());
    {
        // If triSurface do additional zoning based on connectivity
        for (const label surfi : blockedSurfaces)
        {
            const label geomi = surfaces_.surfaces()[surfi];
            const searchableSurface& s = surfaces_.geometry()[geomi];
            if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            {
                const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
                const labelListList& edFaces = surf.edgeFaces();
                boolList isOpenEdge(edFaces.size(), false);
                forAll(edFaces, i)
                {
                    if (edFaces[i].size() == 1)
                    {
                        isOpenEdge[i] = true;
                    }
                }

                labelList faceZone;
                const label nZones = surf.markZones(isOpenEdge, faceZone);
                if (nZones > 1)
                {
                    faceZones[surfi].transfer(faceZone);
                }
            }
        }
    }


    // Clever limiting of the number of iterations (= max cells in the channel)
    // since it has too many problematic issues, e.g. with volume refinement
    // and the real check uses the proper distance anyway just disable.
    const label nIters = mesh_.globalData().nTotalCells();


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined)
    //const labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    calcCellCellRays
    (
        neiCc,
        neiLevel,
        testFaces,
        start,
        end,
        minLevel
    );
    // TBD. correct nIters for actual cellLevel (since e.g. volume refinement
    //      might add to cell level). Unfortunately we only have minLevel,
    //      not maxLevel. Problem: what if volume refinement only in middle of
    //      channel? Say channel 1m wide with a 0.1m sphere of refinement
    //      Workaround: have dummy surface with e.g. maxLevel 100 to
    //      force nIters to be high enough.


    // Test for all intersections (with surfaces of higher gap level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    vectorField normal1;

    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    vectorField normal2;

    surfaces_.findNearestIntersection
    (
        blockedSurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,    // local region
        normal1,

        surface2,
        hit2,
        region2,    // local region
        normal2
    );


    // Detect cells that are using multiple surface regions
    //bitSet isMultiRegion;
    //detectMultiRegionCells
    //(
    //    faceZones,
    //    testFaces,
    //
    //    surface1,
    //    hit1,
    //    region1,
    //
    //    surface2,
    //    hit2,
    //    region2,
    //
    //    isMultiRegion
    //);


    label n = 0;
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            n++;
        }
    }

    List<wallPoints> faceDist(n);
    labelList changedFaces(n);
    n = 0;

    DynamicList<point> originLocation(2);
    DynamicList<scalar> originDistSqr(2);
    DynamicList<FixedList<label, 3>> originSurface(2);
    //DynamicList<point> originNormal(2);


    //- To avoid walking through surfaces we mark all faces that have been
    //  intersected. We can either mark only those faces intersecting
    //  blockedSurfaces (i.e. with a 'blockLevel') or mark faces intersecting
    //  any (refinement) surface (this includes e.g. faceZones). This is
    //  much easier since that information is already cached
    //  (meshRefinement::intersectedFaces())

    //bitSet isBlockedFace(mesh_.nFaces());
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            const label facei = testFaces[i];
            //isBlockedFace.set(facei);
            const point& fc = mesh_.faceCentres()[facei];
            const labelList& fz1 = faceZones[surface1[i]];

            originLocation.clear();
            originDistSqr.clear();
            //originNormal.clear();
            originSurface.clear();

            originLocation.append(hit1[i].point());
            originDistSqr.append(hit1[i].point().distSqr(fc));
            //originNormal.append(normal1[i]);
            originSurface.append
            (
                FixedList<label, 3>
                ({
                    surface1[i],
                    region1[i],
                    (fz1.size() ? fz1[hit1[i].index()] : region1[i])
                })
            );

            if (hit2[i].hit() && hit1[i] != hit2[i])
            {
                const labelList& fz2 = faceZones[surface2[i]];
                originLocation.append(hit2[i].point());
                originDistSqr.append(hit2[i].point().distSqr(fc));
                //originNormal.append(normal2[i]);
                originSurface.append
                (
                    FixedList<label, 3>
                    ({
                        surface2[i],
                        region2[i],
                        (fz2.size() ? fz2[hit2[i].index()] : region2[i])
                    })
                );
            }

            // Collect all seed data. Currently walking does not look at
            // surface direction - if so pass in surface normal as well
            faceDist[n] = wallPoints
            (
                originLocation,     // origin
                originDistSqr,      // distance to origin
                originSurface       // surface+region+zone
                //originNormal        // normal at origin
            );
            changedFaces[n] = facei;
            n++;
        }
    }


    // Clear intersection info
    surface1.clear();
    hit1.clear();
    region1.clear();
    normal1.clear();
    surface2.clear();
    hit2.clear();
    region2.clear();
    normal2.clear();


    List<wallPoints> allFaceInfo(mesh_.nFaces());
    List<wallPoints> allCellInfo(mesh_.nCells());

    // Any refinement surface (even a faceZone) should stop the gap walking.
    // This is exactly the information which is cached in the surfaceIndex_
    // field.
    const bitSet isBlockedFace(intersectedFaces());

    wallPoints::trackData td(isBlockedFace, regionToBlockSize);
    FaceCellWave<wallPoints, wallPoints::trackData> wallDistCalc
    (
        mesh_,
        changedFaces,
        faceDist,
        allFaceInfo,
        allCellInfo,
        0,            // max iterations
        td
    );
    wallDistCalc.iterate(nIters);


    // Extract faces of visited cells

    bitSet isProximityFace(mesh_.nFaces(), false);

    // Make sure to include initial intersected ones
    isProximityFace.set(testFaces);

    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].valid(wallDistCalc.data()))
        {
            isProximityFace.set(mesh_.cells()[celli]);
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        isProximityFace,
        orEqOp<unsigned int>(),
        0u
    );

    testFaces = isProximityFace.toc();
}


Foam::label Foam::meshRefinement::markFakeGapRefinement
(
    const scalar planarCos,
    const labelList& blockedSurfaces,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    // Re-work the blockLevel of the blockedSurfaces into a length scale
    // and a number of cells to cross
    List<scalarList> regionToBlockSize(surfaces_.surfaces().size());
    scalar maxBlockSize(-1);
    for (const label surfi : blockedSurfaces)
    {
        const label geomi = surfaces_.surfaces()[surfi];
        const searchableSurface& s = surfaces_.geometry()[geomi];
        const label nRegions = s.regions().size();

        regionToBlockSize[surfi].setSize(nRegions, -1);
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            const label globalRegioni = surfaces_.globalRegion(surfi, regioni);
            const label bLevel = surfaces_.blockLevel()[globalRegioni];

            // If valid blockLevel specified
            if (bLevel != labelMax)
            {
                regionToBlockSize[surfi][regioni] =
                    meshCutter_.level0EdgeLength()/pow(2.0, bLevel);
            }

            maxBlockSize = max(maxBlockSize, regionToBlockSize[surfi][regioni]);

            //const label mLevel = surfaces_.maxLevel()[globalRegioni];
            //// TBD: check for higher cached level of surface due to vol
            ////      refinement. Problem: might still miss refinement bubble
            ////      fully inside thin channel
            //if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            //{
            //    const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
            //}
        }
    }

    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined) and their cells
    labelList cellMap;
    {
        labelList testFaces(getRefineCandidateFaces(refineCell));

        // Extend the set by faceCellFace wave upto given 3*blockLevel size
        markProximityCandidateFaces
        (
            blockedSurfaces,
            regionToBlockSize,
            neiLevel,
            neiCc,
            testFaces
        );


        bitSet isTestCell(mesh_.nCells());
        for (const label facei : testFaces)
        {
            const label own = mesh_.faceOwner()[facei];
            if (refineCell[own] == -1)
            {
                isTestCell.set(own);
            }
            if (mesh_.isInternalFace(facei))
            {
                const label nei = mesh_.faceNeighbour()[facei];
                if (refineCell[nei] == -1)
                {
                    isTestCell.set(nei);
                }
            }
        }
        cellMap = isTestCell.sortedToc();
    }


    // Shoot rays in all 6 directions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Almost exact copy of meshRefinement::markInternalGapRefinement.
    // Difference is only in the length scale (now based on blockLevel)


    // Find per cell centre the nearest point and normal on the surfaces
    List<pointIndexHit> nearInfo;
    vectorField nearNormal;
    labelList nearSurface;
    labelList nearRegion;
    {
        surfaces_.findNearestRegion
        (
            blockedSurfaces,
            pointField(mesh_.cellCentres(), cellMap),
            scalarField(cellMap.size(), maxBlockSize),  // use uniform for now
            nearSurface,
            nearInfo,
            nearRegion,
            nearNormal
        );
    }


    DynamicList<label> map(nearInfo.size());
    DynamicField<point> rayStart(nearInfo.size());
    DynamicField<point> rayEnd(nearInfo.size());
    DynamicField<scalar> gapSize(nearInfo.size());

    DynamicField<point> rayStart2(nearInfo.size());
    DynamicField<point> rayEnd2(nearInfo.size());
    DynamicField<scalar> gapSize2(nearInfo.size());

    label nTestCells = 0;

    forAll(nearInfo, i)
    {
        if (nearInfo[i].hit())
        {
            const label globalRegioni = surfaces_.globalRegion
            (
                nearSurface[i],
                nearRegion[i]
            );

            // Combine info from shell and surface
            const label bLevel = surfaces_.blockLevel()[globalRegioni];
            const FixedList<label, 3> gapInfo
            ({
                100,
                bLevel,
                100
            });
            const volumeType gapMode(volumeType::MIXED);

            // Store wanted number of cells in gap
            const label celli = cellMap[i];
            const label cLevel = meshCutter_.cellLevel()[celli];

            // Construct one or more rays to test for oppositeness
            const label nRays = generateRays
            (
                false,                      // ignore actual surface normal
                nearInfo[i].hitPoint(),
                nearNormal[i],
                gapInfo,
                gapMode,

                mesh_.cellCentres()[celli],
                cLevel,

                rayStart,
                rayEnd,
                gapSize,

                rayStart2,
                rayEnd2,
                gapSize2
            );
            if (nRays > 0)
            {
                nTestCells++;
                for (label j = 0; j < nRays; j++)
                {
                    map.append(i);
                }
            }
        }
    }

    Info<< "Selected " << returnReduce(nTestCells, sumOp<label>())
        << " cells for testing out of "
        << mesh_.globalData().nTotalCells() << endl;
    map.shrink();
    rayStart.shrink();
    rayEnd.shrink();
    gapSize.shrink();

    rayStart2.shrink();
    rayEnd2.shrink();
    gapSize2.shrink();

    cellMap = labelUIndList(cellMap, map)();
    nearNormal = UIndirectList<vector>(nearNormal, map)();
    nearInfo.clear();
    nearSurface.clear();
    nearRegion.clear();


    // Do intersections in pairs
    labelList surf1;
    labelList region1;
    List<pointIndexHit> hit1;
    vectorField normal1;
    surfaces_.findNearestIntersection
    (
        rayStart,
        rayEnd,
        surf1,
        region1,
        hit1,
        normal1
    );

    labelList surf2;
    labelList region2;
    List<pointIndexHit> hit2;
    vectorField normal2;
    surfaces_.findNearestIntersection
    (
        rayStart2,
        rayEnd2,
        surf2,
        region2,
        hit2,
        normal2
    );

    // Extract cell based gap size
    const label oldNRefine = nRefine;
    forAll(surf1, i)
    {
        if (surf1[i] != -1 && surf2[i] != -1)
        {
            // Found intersections with surface. Check for
            // - small gap
            // - coplanar normals

            const label celli = cellMap[i];
            if (celli != -1 && refineCell[celli] == -1)
            {
                const scalar d2 = magSqr(hit1[i].hitPoint()-hit2[i].hitPoint());

                const scalar maxGapSize
                (
                    max
                    (
                        regionToBlockSize[surf1[i]][region1[i]],
                        regionToBlockSize[surf2[i]][region2[i]]
                    )
                );

                if
                (
                    (mag(normal1[i]&normal2[i]) > planarCos)
                 && (d2 < Foam::sqr(maxGapSize))
                )
                {
                    if
                    (
                        !markForRefine
                        (
                            0,                      // mark level
                            nAllowRefine,
                            refineCell[celli],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeGapCells
(
    const scalar planarAngle,
    const labelList& minSurfaceLevel,
    const labelList& globalToMasterPatch,
    const label growIter
)
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList refineCell(mesh_.nCells(), -1);
    label nRefine = 0;
    //markProximityRefinement
    //(
    //    Foam::cos(degToRad(planarAngle)),
    //
    //    minSurfaceLevel,                                // surface min level
    //    labelList(minSurfaceLevel.size(), labelMax),    // surfaceGapLevel
    //
    //    labelMax/Pstream::nProcs(), //nAllowRefine,
    //    neiLevel,
    //    neiCc,
    //
    //    refineCell,
    //    nRefine
    //);


    // Determine minimum blockLevel per surface
    Map<label> surfToBlockLevel;

    forAll(surfaces_.surfaces(), surfi)
    {
        const label geomi = surfaces_.surfaces()[surfi];
        const searchableSurface& s = surfaces_.geometry()[geomi];
        const label nRegions = s.regions().size();

        label minBlockLevel = labelMax;
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            const label globalRegioni = surfaces_.globalRegion(surfi, regioni);
            minBlockLevel = min
            (
                minBlockLevel,
                surfaces_.blockLevel()[globalRegioni]
            );
        }

        if (minBlockLevel < labelMax)
        {
            surfToBlockLevel.insert(surfi, minBlockLevel);
        }
    }


    //markProximityRefinementWave
    //(
    //    Foam::cos(degToRad(planarAngle)),
    //    surfToBlockLevel.sortedToc(),
    //
    //    labelMax/Pstream::nProcs(), //nAllowRefine,
    //    neiLevel,
    //    neiCc,
    //
    //    refineCell,
    //    nRefine
    //);

    markFakeGapRefinement
    (
        Foam::cos(degToRad(planarAngle)),
        surfToBlockLevel.sortedToc(),

        labelMax/Pstream::nProcs(), //nAllowRefine,
        neiLevel,
        neiCc,

        refineCell,
        nRefine
    );


    Info<< "Marked for blocking due to close opposite surfaces         : "
        << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;

    // Remove outliers, i.e. cells with all points exposed
    if (growIter)
    {
        labelList oldRefineCell(refineCell);

        // Pass1: extend the set to fill in gaps
        bitSet isOutsideFace;
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces
            markOutsideFaces
            (
                meshCutter_.cellLevel(),
                neiLevel,
                refineCell,
                isOutsideFace
            );
            // Extend with cells with three outside faces
            growSet(neiLevel, isOutsideFace, refineCell, nRefine);
        }


        // Pass2: erode back to original set if pass1 didn't help
        for (label iter = 0; iter < growIter; iter++)
        {
            // Get outside faces. Ignore cell level.
            markOutsideFaces
            (
                labelList(mesh_.nCells(), 0),
                labelList(neiLevel.size(), 0),
                refineCell,
                isOutsideFace
            );

            // Unmark cells with three or more outside faces
            for (label celli = 0; celli < mesh_.nCells(); celli++)
            {
                if (refineCell[celli] != -1 && oldRefineCell[celli] == -1)
                {
                    if (countFaceDirs(isOutsideFace, celli) >= 3)
                    {
                        refineCell[celli] = -1;
                        --nRefine;
                    }
                }
            }
        }

        Info<< "Marked for blocking after filtering                        : "
            << returnReduce(nRefine, sumOp<label>()) << " cells." << endl;
    }


    // Determine patch for every mesh face
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));
    const label defaultRegion(surfaces_.globalRegion(unnamedSurfaces[0], 0));

    const labelList nearestRegion
    (
        nearestIntersection
        (
            unnamedSurfaces,
            defaultRegion
        )
    );

    // Pack
    labelList cellsToRemove(nRefine);
    nRefine = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell[cellI] != -1)
        {
            cellsToRemove[nRefine++] = cellI;
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    labelList exposedPatches(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        exposedPatches[i] = globalToMasterPatch[nearestRegion[facei]];
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


void Foam::meshRefinement::selectIntersectedFaces
(
    const labelList& selectedSurfaces,
    boolList& isBlockedFace
) const
{
    // Like meshRefinement::selectSeparatedCoupledFaces. tbd: convert to bitSet

    // Check that no connection between inside and outside points
    isBlockedFace.setSize(mesh_.nFaces(), false);

    // Block off separated couples.
    selectSeparatedCoupledFaces(isBlockedFace);

    // Block off intersections with selected surfaces

    // Mark per face (for efficiency)
    boolList isSelectedSurf(surfaces_.surfaces().size(), false);
    UIndirectList<bool>(isSelectedSurf, selectedSurfaces) = true;

    forAll(surfaceIndex_, facei)
    {
        const label surfi = surfaceIndex_[facei];
        if (surfi != -1 && isSelectedSurf[surfi])
        {
            isBlockedFace[facei] = true;
        }
    }
}


//Foam::labelList Foam::meshRefinement::detectLeakCells
//(
//    const boolList& isBlockedFace,
//    const labelList& leakFaces,
//    const labelList& seedCells
//) const
//{
//    int dummyTrackData = 0;
//    List<topoDistanceData<label>> allFaceInfo(mesh_.nFaces());
//    List<topoDistanceData<label>> allCellInfo(mesh_.nCells());
//
//    // Block faces
//    forAll(isBlockedFace, facei)
//    {
//        if (isBlockedFace[facei])
//        {
//            allFaceInfo[facei] = topoDistanceData<label>(labelMax, 123);
//        }
//    }
//    for (const label facei : leakFaces)
//    {
//        allFaceInfo[facei] = topoDistanceData<label>(labelMax, 123);
//    }
//
//    // Walk from inside cell
//    DynamicList<topoDistanceData<label>> faceDist;
//    DynamicList<label> cFaces1;
//    for (const label celli : seedCells)
//    {
//        if (celli != -1)
//        {
//            const labelList& cFaces = mesh_.cells()[celli];
//            faceDist.reserve(cFaces.size());
//            cFaces1.reserve(cFaces.size());
//
//            for (const label facei : cFaces)
//            {
//                if (!allFaceInfo[facei].valid(dummyTrackData))
//                {
//                    cFaces1.append(facei);
//                    faceDist.append(topoDistanceData<label>(0, 123));
//                }
//            }
//        }
//    }
//
//    // Walk through face-cell wave till all cells are reached
//    FaceCellWave<topoDistanceData<label>> wallDistCalc
//    (
//        mesh_,
//        cFaces1,
//        faceDist,
//        allFaceInfo,
//        allCellInfo,
//        mesh_.globalData().nTotalCells()+1   // max iterations
//    );
//
//    label nRemove = 0;
//    for (const label facei : leakFaces)
//    {
//        const label own = mesh_.faceOwner()[facei];
//        if (!allCellInfo[own].valid(dummyTrackData))
//        {
//            nRemove++;
//        }
//        if (mesh_.isInternalFace(facei))
//        {
//            const label nei = mesh_.faceNeighbour()[facei];
//            if (!allCellInfo[nei].valid(dummyTrackData))
//            {
//                nRemove++;
//            }
//        }
//    }
//
//    labelList cellsToRemove(nRemove);
//    nRemove = 0;
//    for (const label facei : leakFaces)
//    {
//        const label own = mesh_.faceOwner()[facei];
//        if (!allCellInfo[own].valid(dummyTrackData))
//        {
//            cellsToRemove[nRemove++] = own;
//        }
//        if (mesh_.isInternalFace(facei))
//        {
//            const label nei = mesh_.faceNeighbour()[facei];
//            if (!allCellInfo[nei].valid(dummyTrackData))
//            {
//                cellsToRemove[nRemove++] = nei;
//            }
//        }
//    }
//
//    if (debug)
//    {
//        volScalarField fld
//        (
//            IOobject
//            (
//                "cellsToKeep",
//                mesh_.time().timeName(),
//                mesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh_,
//            dimensionedScalar(dimless, Zero)
//        );
//        forAll(allCellInfo, celli)
//        {
//            if (allCellInfo[celli].valid(dummyTrackData))
//            {
//                fld[celli] = allCellInfo[celli].distance();
//            }
//        }
//        forAll(fld.boundaryField(), patchi)
//        {
//            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
//            SubList<topoDistanceData<label>> p(pp.patchSlice(allFaceInfo));
//            scalarField pfld
//            (
//                fld.boundaryField()[patchi].size(),
//                Zero
//            );
//            forAll(pfld, i)
//            {
//                if (p[i].valid(dummyTrackData))
//                {
//                    pfld[i] = 1.0*p[i].distance();
//                }
//            }
//            fld.boundaryFieldRef()[patchi] == pfld;
//        }
//        //Note: do not swap cell values so do not do
//        //fld.correctBoundaryConditions();
//        Pout<< "Writing distance field for initial cells "
//            << seedCells << " to " << fld.objectPath() << endl;
//        fld.write();
//    }
//
//    return cellsToRemove;
//}
//
//
//Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeLeakCells
//(
//    const labelList& globalToMasterPatch,
//    const labelList& globalToSlavePatch,
//    const pointField& locationsInMesh,
//    const wordList& zonesInMesh,
//    const pointField& locationsOutsideMesh,
//    const labelList& selectedSurfaces
//)
//{
//    boolList isBlockedFace;
//    selectIntersectedFaces(selectedSurfaces, isBlockedFace);
//
//    // Determine cell regions
//    const regionSplit cellRegion(mesh_, isBlockedFace);
//
//    // Detect locationsInMesh regions
//    labelList insideCells(locationsInMesh.size(), -1);
//    labelList insideRegions(locationsInMesh.size(), -1);
//    forAll(locationsInMesh, i)
//    {
//        insideCells[i] = findCell
//        (
//            mesh_,
//            mergeDistance_*vector::one,   //perturbVec,
//            locationsInMesh[i]
//        );
//        if (insideCells[i] != -1)
//        {
//            insideRegions[i] = cellRegion[insideCells[i]];
//        }
//        reduce(insideRegions[i], maxOp<label>());
//
//        if (insideRegions[i] == -1)
//        {
//            // See if we can perturb a bit
//            insideCells[i] = findCell
//            (
//                mesh_,
//                mergeDistance_*vector::one,   //perturbVec,
//                locationsInMesh[i]+mergeDistance_*vector::one
//            );
//            if (insideCells[i] != -1)
//            {
//                insideRegions[i] = cellRegion[insideCells[i]];
//            }
//            reduce(insideRegions[i], maxOp<label>());
//
//            if (insideRegions[i] == -1)
//            {
//                FatalErrorInFunction
//                    << "Cannot find locationInMesh " << locationsInMesh[i]
//                    << " on any processor" << exit(FatalError);
//            }
//        }
//    }
//
//
//    // Check that all the locations outside the
//    // mesh do not conflict with those inside
//
//    bool haveLeak = false;
//    forAll(locationsOutsideMesh, i)
//    {
//        // Find the region containing the point
//        label regioni = findRegion
//        (
//            mesh_,
//            cellRegion,
//            mergeDistance_*vector::one,   //perturbVec,
//            locationsOutsideMesh[i]
//        );
//
//        if (regioni != -1)
//        {
//            // Check for locationsOutsideMesh overlapping with inside ones
//            if (insideRegions.find(regioni) != -1)
//            {
//                haveLeak = true;
//                WarningInFunction
//                    << "Outside location " << locationsOutsideMesh[i]
//                    << " in region " << regioni
//                    << " is connected to one of the inside points "
//                    << locationsInMesh << endl;
//            }
//        }
//    }
//
//
//    autoPtr<mapPolyMesh> mapPtr;
//    if (returnReduceOr(haveLeak))
//    {
//        // Use shortestPathSet to provide a minimum set of faces needed
//        // to close hole. Tbd: maybe directly use wave?
//        meshSearch searchEngine(mesh_);
//        shortestPathSet leakPath
//        (
//            "leakPath",
//            mesh_,
//            searchEngine,
//            coordSet::coordFormatNames[coordSet::coordFormat::DISTANCE],
//            true,
//            50,    // tbd. Number of iterations. This is the maximum
//                    // number of faces in the leak hole
//
//            //pbm.groupPatchIDs()["wall"],  // patch to grow from
//            meshedPatches(),   // patch to grow from
//
//            locationsInMesh,
//            locationsOutsideMesh,
//            isBlockedFace
//        );
//
//
//        // Use leak path to find minimum set of cells to delete
//        const labelList cellsToRemove
//        (
//            detectLeakCells
//            (
//                isBlockedFace,
//                leakPath.leakFaces(),
//                insideCells
//            )
//        );
//
//        // Re-do intersections to find nearest unnamed surface
//        const label defaultRegion
//        (
//            surfaces().globalRegion
//            (
//                selectedSurfaces[0],
//                0
//            )
//        );
//
//        const labelList nearestRegion
//        (
//            nearestIntersection
//            (
//                selectedSurfaces,
//                defaultRegion
//            )
//        );
//
//
//        // Remove cells
//        removeCells cellRemover(mesh_);
//        const labelList exposedFaces
//        (
//            cellRemover.getExposedFaces(cellsToRemove)
//        );
//
//        labelList exposedPatches(exposedFaces.size());
//        forAll(exposedFaces, i)
//        {
//            label facei = exposedFaces[i];
//            exposedPatches[i] = globalToMasterPatch[nearestRegion[facei]];
//        }
//
//        mapPtr = doRemoveCells
//        (
//            cellsToRemove,
//            exposedFaces,
//            exposedPatches,
//            cellRemover
//        );
//
//
//        // Put the exposed faces into a special faceZone
//        {
//            // Add to "frozenFaces" zone
//            faceZoneMesh& faceZones = mesh_.faceZones();
//
//            // Get current frozen faces (if any)
//            bitSet isFrozenFace(mesh_.nFaces());
//            label zonei = faceZones.findZoneID("frozenFaces");
//            if (zonei != -1)
//            {
//                const bitSet oldSet(mesh_.nFaces(), faceZones[zonei]);
//                isFrozenFace.set(oldSet);
//            }
//
//            // Add newly exposed faces (if not yet in any faceZone!)
//            const labelList exposed
//            (
//                renumber
//                (
//                    mapPtr().reverseFaceMap(),
//                    exposedFaces
//                )
//            );
//            bitSet isZonedFace(mesh_.nFaces(), faceZones.zoneMap().toc());
//            for (const label facei : exposed)
//            {
//                if (!isZonedFace[facei])
//                {
//                    isFrozenFace.set(facei);
//                }
//            }
//
//            syncTools::syncFaceList
//            (
//                mesh_,
//                isFrozenFace,
//                orEqOp<unsigned int>(),
//                0u
//            );
//
//            // Add faceZone if non existing
//            faceZones.clearAddressing();
//            if (zonei == -1)
//            {
//                zonei = faceZones.size();
//
//                faceZones.emplace_back
//                (
//                    "frozenFaces",    // name
//                    zonei,            // index
//                    faceZones
//                );
//            }
//
//            // Update faceZone with new contents
//            labelList frozenFaces(isFrozenFace.toc());
//            boolList frozenFlip(frozenFaces.size(), false);
//
//            faceZones[zonei].resetAddressing
//            (
//                std::move(frozenFaces),
//                std::move(frozenFlip)
//            );
//        }
//
//
//        //// Put the exposed points into a special pointZone
//        //if (false)
//        //{
//        //    const labelList meshFaceIDs
//        //    (
//        //        renumber
//        //        (
//        //            mapPtr().reverseFaceMap(),
//        //            exposedFaces
//        //        )
//        //    );
//        //    const uindirectPrimitivePatch pp
//        //    (
//        //        UIndirectList<face>(mesh_.faces(), meshFaceIDs),
//        //        mesh_.points()
//        //    );
//        //
//        //    // Count number of faces per edge
//        //    const labelListList& edgeFaces = pp.edgeFaces();
//        //    labelList nEdgeFaces(edgeFaces.size());
//        //    forAll(edgeFaces, edgei)
//        //    {
//        //        nEdgeFaces[edgei] = edgeFaces[edgei].size();
//        //    }
//        //
//        //    // Sync across processor patches
//        //    if (Pstream::parRun())
//        //    {
//        //        const globalMeshData& globalData = mesh_.globalData();
//        //        const mapDistribute& map = globalData.globalEdgeSlavesMap();
//        //        const indirectPrimitivePatch& cpp =
//        //            globalData.coupledPatch();
//        //
//        //        // Match pp edges to coupled edges
//        //        labelList patchEdges;
//        //        labelList coupledEdges;
//        //        PackedBoolList sameEdgeOrientation;
//        //        PatchTools::matchEdges
//        //        (
//        //            pp,
//        //            cpp,
//        //            patchEdges,
//        //            coupledEdges,
//        //            sameEdgeOrientation
//        //        );
//        //
//        //
//        //        // Convert patch-edge data into cpp-edge data
//        //        labelList coupledNEdgeFaces(map.constructSize(), Zero);
//        //        UIndirectList<label>(coupledNEdgeFaces, coupledEdges) =
//        //            UIndirectList<label>(nEdgeFaces, patchEdges);
//        //
//        //        // Synchronise
//        //        globalData.syncData
//        //        (
//        //            coupledNEdgeFaces,
//        //            globalData.globalEdgeSlaves(),
//        //            globalData.globalEdgeTransformedSlaves(),
//        //            map,
//        //            plusEqOp<label>()
//        //        );
//        //
//        //        // Convert back from cpp-edge to patch-edge
//        //        UIndirectList<label>(nEdgeFaces, patchEdges) =
//        //            UIndirectList<label>(coupledNEdgeFaces, coupledEdges);
//        //    }
//        //
//        //    // Freeze all internal points
//        //    bitSet isFrozenPoint(mesh_.nPoints());
//        //    forAll(nEdgeFaces, edgei)
//        //    {
//        //        if (nEdgeFaces[edgei] != 1)
//        //        {
//        //            const edge& e = pp.edges()[edgei];
//        //            isFrozenPoint.set(pp.meshPoints()[e[0]]);
//        //            isFrozenPoint.set(pp.meshPoints()[e[1]]);
//        //        }
//        //    }
//        //    // Add to "frozenPoints" zone
//        //    pointZoneMesh& pointZones = mesh_.pointZones();
//        //    pointZones.clearAddressing();
//        //    label zonei = pointZones.findZoneID("frozenPoints");
//        //    if (zonei != -1)
//        //    {
//        //        const bitSet oldSet(mesh_.nPoints(), pointZones[zonei]);
//        //        // Add to isFrozenPoint
//        //        isFrozenPoint.set(oldSet);
//        //    }
//        //
//        //    syncTools::syncPointList
//        //    (
//        //        mesh_,
//        //        isFrozenPoint,
//        //        orEqOp<unsigned int>(),
//        //        0u
//        //    );
//        //
//        //    if (zonei == -1)
//        //    {
//        //        zonei = pointZones.size();
//        //
//        //        pointZones.emplace_back
//        //        (
//        //            "frozenPoints",             // name
//        //            isFrozenPoint.sortedToc(),  // addressing
//        //            zonei,                      // index
//        //            pointZones                  // pointZoneMesh
//        //        );
//        //    }
//        //}
//
//
//        if (debug&meshRefinement::MESH)
//        {
//            const_cast<Time&>(mesh_.time())++;
//
//            Pout<< "Writing current mesh to time "
//                << timeName() << endl;
//            write
//            (
//                meshRefinement::debugType(debug),
//                meshRefinement::writeType
//                (
//                    meshRefinement::writeLevel()
//                  | meshRefinement::WRITEMESH
//                ),
//                mesh_.time().path()/timeName()
//            );
//            Pout<< "Dumped mesh in = "
//                << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
//        }
//    }
//    return mapPtr;
//}


//Foam::autoPtr<Foam::mapDistribute> Foam::meshRefinement::nearestFace
//(
//    const globalIndex& globalSeedFaces,
//    const labelList& seedFaces,
//    const labelList& closureFaces
//) const
//{
//    // Takes a set of faces for which we have information (seedFaces,
//    // globalSeedFaces - these are e.g. intersected faces) and walks out
//    // across nonSeedFace. Returns for
//    // every nonSeedFace the nearest seed face (in global indexing).
//    // Used e.g. in hole closing. Assumes that seedFaces, closureFaces
//    // are a small subset of the master
//    // so are not large - it uses edge addressing on the closureFaces
//
//
//    if (seedFaces.size() != globalSeedFaces.localSize())
//    {
//        FatalErrorInFunction << "problem : seedFaces:" << seedFaces.size()
//            << " globalSeedFaces:" << globalSeedFaces.localSize()
//            << exit(FatalError);
//    }
//
//    //// Mark mesh points that are used by any closureFaces. This is for
//    //// initial filtering
//    //bitSet isNonSeedPoint(mesh.nPoints());
//    //for (const label facei : closureFaces)
//    //{
//    //    isNonSeedPoint.set(mesh_.faces()[facei]);
//    //}
//    //syncTools::syncPointList
//    //(
//    //    mesh_,
//    //    isNonSeedPoint,
//    //    orEqOp<unsigned int>(),
//    //    0u
//    //);
//
//    // Make patch
//    const uindirectPrimitivePatch pp
//    (
//        IndirectList<face>(mesh_.faces(), closureFaces),
//        mesh_.points()
//    );
//    const edgeList& edges = pp.edges();
//    const labelList& mp = pp.meshPoints();
//    const label nBndEdges = pp.nEdges() - pp.nInternalEdges();
//
//    // For all faces in seedFaces mark the edge with a face. No special
//    // handling for multiple faces sharing the edge - first one wins
//    EdgeMap<label> edgeMap(pp.nEdges());
//    for (const label facei : seedFaces)
//    {
//        const label globalFacei = globalSeedFaces.toGlobal(facei);
//        const face& f = mesh_.faces()[facei];
//        forAll(f, fp)
//        {
//            label nextFp = f.fcIndex(fp);
//            edgeMap.insert(edge(f[fp], f[nextFp]), globalFacei);
//        }
//    }
//    syncTools::syncEdgeMap(mesh_, edgeMap, maxEqOp<label>());
//
//
//
//    // Seed
//    DynamicList<label> initialEdges(2*nBndEdges);
//    DynamicList<edgeTopoDistanceData<label, uindirectPrimitivePatch>>
//        initialEdgesInfo(2*nBndEdges);
//    forAll(edges, edgei)
//    {
//        const edge& e = edges[edgei];
//        const edge meshE = edge(mp[e[0]], mp[e[1]]);
//
//        const auto iter = edgeMap.cfind(meshE);
//        if (iter.good())
//        {
//            initialEdges.append(edgei);
//            initialEdgesInfo.append
//            (
//                edgeTopoDistanceData<label, uindirectPrimitivePatch>
//                (
//                    0,          // distance
//                    iter.val()  // globalFacei
//                )
//            );
//        }
//    }
//
//    // Data on all edges and faces
//    List<edgeTopoDistanceData<label, uindirectPrimitivePatch>> allEdgeInfo
//    (
//        pp.nEdges()
//    );
//    List<edgeTopoDistanceData<label, uindirectPrimitivePatch>> allFaceInfo
//    (
//        pp.size()
//    );
//
//    // Walk
//    PatchEdgeFaceWave
//    <
//        uindirectPrimitivePatch,
//        edgeTopoDistanceData<label, uindirectPrimitivePatch>
//    > calc
//    (
//        mesh_,
//        pp,
//        initialEdges,
//        initialEdgesInfo,
//        allEdgeInfo,
//        allFaceInfo,
//        returnReduce(pp.nEdges(), sumOp<label>())
//    );
//
//
//    // Per non-seed face the seed face
//    labelList closureToSeed(pp.size(), -1);
//    forAll(allFaceInfo, facei)
//    {
//        if (allFaceInfo[facei].valid(calc.data()))
//        {
//            closureToSeed[facei] = allFaceInfo[facei].data();
//        }
//    }
//
//    List<Map<label>> compactMap;
//    return autoPtr<mapDistribute>::New
//    (
//        globalSeedFaces,
//        closureToSeed,
//        compactMap
//    );
//}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::blockLeakFaces
(
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh,
    const labelList& selectedSurfaces
)
{
    // Problem: this is based on cached intersection information. This might
    // not include the surface we actually want to use. In which case the
    // surface would not be seen as intersected!
    boolList isBlockedFace;
    selectIntersectedFaces(selectedSurfaces, isBlockedFace);

    // Determine cell regions
    const regionSplit cellRegion(mesh_, isBlockedFace);

    // Detect locationsInMesh regions
    labelList insideCells(locationsInMesh.size(), -1);
    labelList insideRegions(locationsInMesh.size(), -1);
    forAll(locationsInMesh, i)
    {
        insideCells[i] = findCell
        (
            mesh_,
            mergeDistance_*vector::one,   //perturbVec,
            locationsInMesh[i]
        );
        if (insideCells[i] != -1)
        {
            insideRegions[i] = cellRegion[insideCells[i]];
        }
        reduce(insideRegions[i], maxOp<label>());

        if (insideRegions[i] == -1)
        {
            // See if we can perturb a bit
            insideCells[i] = findCell
            (
                mesh_,
                mergeDistance_*vector::one,   //perturbVec,
                locationsInMesh[i]+mergeDistance_*vector::one
            );
            if (insideCells[i] != -1)
            {
                insideRegions[i] = cellRegion[insideCells[i]];
            }
            reduce(insideRegions[i], maxOp<label>());

            if (insideRegions[i] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find locationInMesh " << locationsInMesh[i]
                    << " on any processor" << exit(FatalError);
            }
        }
    }


    // Check that all the locations outside the
    // mesh do not conflict with those inside

    bool haveLeak = false;
    forAll(locationsOutsideMesh, i)
    {
        // Find the region containing the point
        label regioni = findRegion
        (
            mesh_,
            cellRegion,
            mergeDistance_*vector::one,   //perturbVec,
            locationsOutsideMesh[i]
        );

        if (regioni != -1)
        {
            // Check for locationsOutsideMesh overlapping with inside ones
            if (insideRegions.find(regioni) != -1)
            {
                haveLeak = true;
                WarningInFunction
                    << "Outside location " << locationsOutsideMesh[i]
                    << " in region " << regioni
                    << " is connected to one of the inside points "
                    << locationsInMesh << endl;
            }
        }
    }


    autoPtr<mapPolyMesh> mapPtr;
    if (returnReduceOr(haveLeak))
    {
        // Use holeToFace to provide a minimum set of faces needed
        // to close hole.

        const List<pointField> allLocations
        (
            refinementParameters::zonePoints
            (
                locationsInMesh,
                zonesInMesh,
                locationsOutsideMesh
            )
        );

        const labelList blockingFaces(findIndices(isBlockedFace, true));

        labelList closureFaces;
        labelList closureToBlocked;
        autoPtr<mapDistribute> closureMapPtr;
        {
            const globalIndex globalBlockingFaces(blockingFaces.size());

            closureMapPtr = holeToFace::calcClosure
            (
                mesh_,
                allLocations,
                blockingFaces,
                globalBlockingFaces,
                true,                   // allow erosion

                closureFaces,
                closureToBlocked
            );

            if (debug)
            {
                Pout<< "meshRefinement::blockLeakFaces :"
                    << " found closure faces:" << closureFaces.size()
                    << "  map:" << bool(closureMapPtr) << endl;
            }

            if (!closureMapPtr)
            {
                FatalErrorInFunction
                    << "have leak but did not find any closure faces"
                    << exit(FatalError);
            }
        }

        // Baffle faces
        labelList ownPatch(mesh_.nFaces(), -1);
        labelList neiPatch(mesh_.nFaces(), -1);

        // Keep (global) boundary faces in their patch
        {
            const polyBoundaryMesh& patches = mesh_.boundaryMesh();
            for (label patchi = 0; patchi < patches.nNonProcessor(); ++patchi)
            {
                const polyPatch& pp = patches[patchi];

                forAll(pp, i)
                {
                    ownPatch[pp.start()+i] = patchi;
                    neiPatch[pp.start()+i] = patchi;
                }
            }
        }

        const faceZoneMesh& fzs = mesh_.faceZones();

        // Mark zone per face
        labelList faceToZone(mesh_.nFaces(), -1);
        boolList faceToFlip(mesh_.nFaces(), false);
        forAll(fzs, zonei)
        {
            const labelList& addressing = fzs[zonei];
            const boolList& flipMap = fzs[zonei].flipMap();

            forAll(addressing, i)
            {
                faceToZone[addressing[i]] = zonei;
                faceToFlip[addressing[i]] = flipMap[i];
            }
        }


        // Fetch patch and zone info from blockingFaces
        labelList packedOwnPatch(labelUIndList(ownPatch, blockingFaces));
        closureMapPtr->distribute(packedOwnPatch);
        labelList packedNeiPatch(labelUIndList(neiPatch, blockingFaces));
        closureMapPtr->distribute(packedNeiPatch);
        labelList packedZone(labelUIndList(faceToZone, blockingFaces));
        closureMapPtr->distribute(packedZone);
        boolList packedFlip(UIndirectList<bool>(faceToFlip, blockingFaces));
        closureMapPtr->distribute(packedFlip);
        forAll(closureFaces, i)
        {
            const label facei = closureFaces[i];
            const label sloti = closureToBlocked[i];

            if (sloti != -1)
            {
                // TBD. how to orient own/nei patch?
                ownPatch[facei] = packedOwnPatch[sloti];
                neiPatch[facei] = packedNeiPatch[sloti];
                faceToZone[facei] = packedZone[sloti];
                faceToFlip[facei] = packedFlip[sloti];
            }
        }


        // Add faces to faceZone. For now do this outside of createBaffles
        // below
        {
            List<DynamicList<label>> zoneToFaces(fzs.size());
            List<DynamicList<bool>> zoneToFlip(fzs.size());

            // Add any to-be-patched face
            forAll(faceToZone, facei)
            {
                const label zonei = faceToZone[facei];
                if (zonei != -1)
                {
                    zoneToFaces[zonei].append(facei);
                    zoneToFlip[zonei].append(faceToFlip[facei]);
                }
            }

            forAll(zoneToFaces, zonei)
            {
                surfaceZonesInfo::addFaceZone
                (
                    fzs[zonei].name(),
                    zoneToFaces[zonei],
                    zoneToFlip[zonei],
                    mesh_
                );
            }
        }

        // Put the points of closureFaces into a special pointZone
        {
            const uindirectPrimitivePatch pp
            (
                UIndirectList<face>(mesh_.faces(), closureFaces),
                mesh_.points()
            );

            // Count number of faces per edge
            const labelList nEdgeFaces(countEdgeFaces(pp));

            // Freeze all internal points
            bitSet isFrozenPoint(mesh_.nPoints());
            forAll(nEdgeFaces, edgei)
            {
                if (nEdgeFaces[edgei] != 1)
                {
                    const edge& e = pp.edges()[edgei];
                    isFrozenPoint.set(pp.meshPoints()[e[0]]);
                    isFrozenPoint.set(pp.meshPoints()[e[1]]);
                }
            }

            // Lookup/add pointZone and include its points
            pointZoneMesh& pointZones =
                const_cast<pointZoneMesh&>(mesh_.pointZones());
            const label zonei = addPointZone("frozenPoints");
            const bitSet oldSet(mesh_.nPoints(), pointZones[zonei]);
            isFrozenPoint.set(oldSet);

            syncTools::syncPointList
            (
                mesh_,
                isFrozenPoint,
                orEqOp<unsigned int>(),
                0u
            );

            // Override addressing
            pointZones.clearAddressing();
            pointZones[zonei] = isFrozenPoint.sortedToc();
        }



        // Create baffles for faces
        mapPtr = createBaffles(ownPatch, neiPatch);

        //// Put the exposed faces into a special faceZone
        //{
        //    // Add newly exposed faces (if not yet in any faceZone!)
        //    const labelList exposed
        //    (
        //        renumber
        //        (
        //            mapPtr().reverseFaceMap(),
        //            blockingFaces
        //        )
        //    );
        //
        //    surfaceZonesInfo::addFaceZone
        //    (
        //        "frozenFaces",
        //        exposed,
        //        boolList(exposed.size(), false),
        //        mesh_
        //    );
        //}


        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh_.time())++;

            Pout<< "Writing current mesh to time "
                << timeName() << endl;
            write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh_.time().path()/timeName()
            );
            Pout<< "Dumped mesh in = "
                << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
        }
    }
    return mapPtr;
}



/*
Foam::label Foam::meshRefinement::markProximityRefinementWave
(
    const scalar planarCos,
    const labelList& blockedSurfaces,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    labelListList faceZones(surfaces_.surfaces().size());
    {
        // If triSurface do additional zoning based on connectivity
        for (const label surfi : blockedSurfaces)
        {
            const label geomi = surfaces_.surfaces()[surfi];
            const searchableSurface& s = surfaces_.geometry()[geomi];
            if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            {
                const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
                const labelListList& edFaces = surf.edgeFaces();
                boolList isOpenEdge(edFaces.size(), false);
                forAll(edFaces, i)
                {
                    if (edFaces[i].size() == 1)
                    {
                        isOpenEdge[i] = true;
                    }
                }

                labelList faceZone;
                const label nZones = surf.markZones(isOpenEdge, faceZone);
                if (nZones > 1)
                {
                    faceZones[surfi].transfer(faceZone);
                }
            }
        }
    }


    // Re-work the blockLevel of the blockedSurfaces into a length scale
    // and a number of cells to cross
    List<scalarList> regionToBlockSize(surfaces_.surfaces().size());
    for (const label surfi : blockedSurfaces)
    {
        const label geomi = surfaces_.surfaces()[surfi];
        const searchableSurface& s = surfaces_.geometry()[geomi];
        const label nRegions = s.regions().size();
        regionToBlockSize[surfi].setSize(nRegions);
        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            const label globalRegioni = surfaces_.globalRegion(surfi, regioni);
            const label bLevel = surfaces_.blockLevel()[globalRegioni];
            regionToBlockSize[surfi][regioni] =
                meshCutter_.level0EdgeLength()/pow(2.0, bLevel);

            //const label mLevel = surfaces_.maxLevel()[globalRegioni];
            //// TBD: check for higher cached level of surface due to vol
            ////      refinement. Problem: might still miss refinement bubble
            ////      fully inside thin channel
            //if (isA<triSurfaceMesh>(s) && !isA<distributedTriSurfaceMesh>(s))
            //{
            //    const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);
            //}

            //nIters = max(nIters, (2<<(mLevel-bLevel)));
        }
    }

    // Clever limiting of the number of iterations (= max cells in the channel)
    // since it has too many problematic issues, e.g. with volume refinement
    // and the real check uses the proper distance anyway just disable.
    const label nIters = mesh_.globalData().nTotalCells();


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined)
    const labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    calcCellCellRays
    (
        neiCc,
        neiLevel,
        testFaces,
        start,
        end,
        minLevel
    );
    // TBD. correct nIters for actual cellLevel (since e.g. volume refinement
    //      might add to cell level). Unfortunately we only have minLevel,
    //      not maxLevel. Problem: what if volume refinement only in middle of
    //      channel? Say channel 1m wide with a 0.1m sphere of refinement
    //      Workaround: have dummy surface with e.g. maxLevel 100 to
    //      force nIters to be high enough.


    // Test for all intersections (with surfaces of higher gap level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    vectorField normal1;

    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    vectorField normal2;

    surfaces_.findNearestIntersection
    (
        blockedSurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,    // local region
        normal1,

        surface2,
        hit2,
        region2,    // local region
        normal2
    );


    // Detect cells that are using multiple surface regions
    //bitSet isMultiRegion;
    //detectMultiRegionCells
    //(
    //    faceZones,
    //    testFaces,
    //
    //    surface1,
    //    hit1,
    //    region1,
    //
    //    surface2,
    //    hit2,
    //    region2,
    //
    //    isMultiRegion
    //);


    label n = 0;
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            n++;
        }
    }

    List<wallPoints> faceDist(n);
    labelList changedFaces(n);
    n = 0;

    DynamicList<point> originLocation(2);
    DynamicList<scalar> originDistSqr(2);
    DynamicList<FixedList<label, 3>> originSurface(2);
    //DynamicList<point> originNormal(2);


    //- To avoid walking through surfaces we mark all faces that have been
    //  intersected. We can either mark only those faces intersecting
    //  blockedSurfaces (i.e. with a 'blockLevel') or mark faces intersecting
    //  any (refinement) surface (this includes e.g. faceZones). This is
    //  much easier since that information is already cached
    //  (meshRefinement::intersectedFaces())

    //bitSet isBlockedFace(mesh_.nFaces());
    forAll(testFaces, i)
    {
        if (hit1[i].hit())
        {
            const label facei = testFaces[i];
            //isBlockedFace.set(facei);
            const point& fc = mesh_.faceCentres()[facei];
            const labelList& fz1 = faceZones[surface1[i]];

            originLocation.clear();
            originDistSqr.clear();
            //originNormal.clear();
            originSurface.clear();

            originLocation.append(hit1[i].point());
            originDistSqr.append(hit1[i].point().distSqr(fc));
            //originNormal.append(normal1[i]);
            originSurface.append
            (
                FixedList<label, 3>
                ({
                    surface1[i],
                    region1[i],
                    (fz1.size() ? fz1[hit1[i].index()] : region1[i])
                })
            );

            if (hit2[i].hit() && hit1[i] != hit2[i])
            {
                const labelList& fz2 = faceZones[surface2[i]];
                originLocation.append(hit2[i].point());
                originDistSqr.append(hit2[i].point().distSqr(fc));
                //originNormal.append(normal2[i]);
                originSurface.append
                (
                    FixedList<label, 3>
                    ({
                        surface2[i],
                        region2[i],
                        (fz2.size() ? fz2[hit2[i].index()] : region2[i])
                    })
                );
            }

            // Collect all seed data. Currently walking does not look at
            // surface direction - if so pass in surface normal as well
            faceDist[n] = wallPoints
            (
                originLocation,     // origin
                originDistSqr,      // distance to origin
                originSurface       // surface+region+zone
                //originNormal        // normal at origin
            );
            changedFaces[n] = facei;
            n++;
        }
    }


    // Clear intersection info
    surface1.clear();
    hit1.clear();
    region1.clear();
    normal1.clear();
    surface2.clear();
    hit2.clear();
    region2.clear();
    normal2.clear();


    List<wallPoints> allFaceInfo(mesh_.nFaces());
    List<wallPoints> allCellInfo(mesh_.nCells());

    // Any refinement surface (even a faceZone) should stop the gap walking.
    // This is exactly the information which is cached in the surfaceIndex_
    // field.
    const bitSet isBlockedFace(intersectedFaces());

    wallPoints::trackData td(isBlockedFace, regionToBlockSize);
    FaceCellWave<wallPoints, wallPoints::trackData> wallDistCalc
    (
        mesh_,
        changedFaces,
        faceDist,
        allFaceInfo,
        allCellInfo,
        0,            // max iterations
        td
    );
    wallDistCalc.iterate(nIters);


    if (debug&meshRefinement::MESH)
    {
        // Dump current nearest opposite surfaces
        volScalarField distance
        (
            IOobject
            (
                "gapSize",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimLength,  //dimensionSet(0, 1, 0, 0, 0),
                -1
            )
        );

        forAll(allCellInfo, celli)
        {
            if (allCellInfo[celli].valid(wallDistCalc.data()))
            {
                const point& cc = mesh_.cellCentres()[celli];
                // Nearest surface points
                const List<point>& origin = allCellInfo[celli].origin();

                // Find 'opposite' pair with minimum distance
                for (label i = 0; i < origin.size(); i++)
                {
                    for (label j = i + 1; j < origin.size(); j++)
                    {
                        if (((cc-origin[i]) & (cc-origin[j])) < 0)
                        {
                            const scalar d(mag(origin[i]-origin[j]));
                            if (distance[celli] < 0)
                            {
                                distance[celli] = d;
                            }
                            else
                            {
                                distance[celli] = min(distance[celli], d);
                            }
                        }
                    }
                }
            }
        }
        distance.correctBoundaryConditions();

        Info<< "Writing measured gap distance to "
            << distance.name() << endl;
        distance.write();
    }



    // Detect tight gaps:
    // - cell is inbetween the two surfaces
    // - two surfaces are planarish
    // - two surfaces are not too far apart
    //   (number of walking iterations is a too-coarse measure)

    scalarField smallGapDistance(mesh_.nCells(), 0.0);
    label nMulti = 0;
    label nSmallGap = 0;

    //OBJstream str(mesh_.time().timePath()/"multiRegion.obj");


    forAll(allCellInfo, celli)
    {
        if (allCellInfo[celli].valid(wallDistCalc.data()))
        {
            const point& cc = mesh_.cellCentres()[celli];

            const List<point>& origin = allCellInfo[celli].origin();
            const List<FixedList<label, 3>>& surface =
                allCellInfo[celli].surface();

            // Find pair with minimum distance
            for (label i = 0; i < origin.size(); i++)
            {
                for (label j = i + 1; j < origin.size(); j++)
                {
                    //if (isMultiRegion[celli])
                    //{
                    //    // The intersection locations are too inaccurate
                    //    // (since not proper nearest, just a cell-cell ray
                    //    //  intersection) so include always
                    //
                    //    smallGapDistance[celli] =
                    //        max(smallGapDistance[celli], maxDist);
                    //
                    //
                    //    str.writeLine(cc, origin[i]);
                    //    str.writeLine(cc, origin[j]);
                    //
                    //    nMulti++;
                    //}
                    //else
                    if (((cc-origin[i]) & (cc-origin[j])) < 0)
                    {
                        const label surfi = surface[i][0];
                        const label regioni = surface[i][1];

                        const label surfj = surface[j][0];
                        const label regionj = surface[j][1];

                        const scalar maxSize = max
                        (
                            regionToBlockSize[surfi][regioni],
                            regionToBlockSize[surfj][regionj]
                        );

                        if
                        (
                            magSqr(origin[i]-origin[j])
                          < Foam::sqr(2*maxSize)
                        )
                        {
                            const scalar maxDist
                            (
                                max
                                (
                                    mag(cc-origin[i]),
                                    mag(cc-origin[j])
                                )
                            );

                            smallGapDistance[celli] =
                                max(smallGapDistance[celli], maxDist);
                            nSmallGap++;
                        }
                    }
                }
            }
        }
    }


    if (debug)
    {
        Info<< "Marked for blocking due to intersecting multiple surfaces  : "
            << returnReduce(nMulti, sumOp<label>()) << " cells." << endl;
        Info<< "Marked for blocking due to close opposite surfaces         : "
            << returnReduce(nSmallGap, sumOp<label>()) << " cells." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        volScalarField distance
        (
            IOobject
            (
                "smallGapDistance",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        );
        distance.field() = smallGapDistance;
        distance.correctBoundaryConditions();

        Info<< "Writing all small-gap cells to "
            << distance.name() << endl;
        distance.write();
    }


    // Mark refinement
    const label oldNRefine = nRefine;
    forAll(smallGapDistance, celli)
    {
        if (smallGapDistance[celli] > SMALL)
        {
            if
            (
                !markForRefine
                (
                    0,                      // mark level
                    nAllowRefine,
                    refineCell[celli],
                    nRefine
                )
            )
            {
                if (debug)
                {
                    Pout<< "Stopped refining since reaching my cell"
                        << " limit of " << mesh_.nCells()+7*nRefine
                        << endl;
                }
                break;
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}
*/

// ************************************************************************* //
