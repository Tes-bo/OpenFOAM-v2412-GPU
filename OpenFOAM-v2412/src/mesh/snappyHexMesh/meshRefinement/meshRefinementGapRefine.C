/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
#include "Time.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "triSurfaceMesh.H"
#include "treeDataCell.H"
#include "searchableSurfaces.H"
#include "DynamicField.H"
#include "transportData.H"
#include "FaceCellWave.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::markSurfaceGapRefinement
(
    const scalar planarCos,

    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Get the gap level for the shells
    const labelList maxLevel(shells_.maxGapLevel());

    label oldNRefine = nRefine;

    if (max(maxLevel) > 0)
    {
        // Use cached surfaceIndex_ to detect if any intersection. If so
        // re-intersect to determine level wanted.

        // Collect candidate faces
        // ~~~~~~~~~~~~~~~~~~~~~~~

        labelList testFaces(getRefineCandidateFaces(refineCell));

        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(testFaces.size());
        pointField end(testFaces.size());
        {
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
        }


        // Collect cells to test for inside/outside in shell
        labelList cellToCompact(mesh_.nCells(), -1);
        labelList bFaceToCompact(mesh_.nBoundaryFaces(), -1);
        labelList gapShell;
        List<FixedList<label, 3>> shellGapInfo;
        List<volumeType> shellGapMode;
        {
            DynamicField<point> compactToCc(mesh_.nCells()/10);
            DynamicList<label> compactToLevel(compactToCc.capacity());
            forAll(testFaces, i)
            {
                label faceI = testFaces[i];
                label own = mesh_.faceOwner()[faceI];
                if (cellToCompact[own] == -1)
                {
                    cellToCompact[own] = compactToCc.size();
                    compactToCc.append(cellCentres[own]);
                    compactToLevel.append(cellLevel[own]);
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (cellToCompact[nei] == -1)
                    {
                        cellToCompact[nei] = compactToCc.size();
                        compactToCc.append(cellCentres[nei]);
                        compactToLevel.append(cellLevel[nei]);
                    }
                }
                else
                {
                    label bFaceI = faceI - mesh_.nInternalFaces();
                    if (bFaceToCompact[bFaceI] == -1)
                    {
                        bFaceToCompact[bFaceI] = compactToCc.size();
                        compactToCc.append(neiCc[bFaceI]);
                        compactToLevel.append(neiLevel[bFaceI]);
                    }
                }
            }

            shells_.findHigherGapLevel
            (
                compactToCc,
                compactToLevel,

                gapShell,
                shellGapInfo,
                shellGapMode
            );
        }


        //const fileName dir(mesh_.time().path()/timeName());
        //if (debug)
        //{
        //    mkDir(dir);
        //    OBJstream insideStr(dir/"insideShell.obj");
        //    OBJstream outsideStr(dir/"outsideShell.obj");
        //    Pout<< "Writing points to:" << nl
        //        << "    inside : " << insideStr.name() << nl
        //        << "    outside: " << outsideStr.name() << nl
        //        << endl;
        //
        //    forAll(cellToCompact, celli)
        //    {
        //        const label compacti = cellToCompact[celli];
        //
        //        if (compacti != -1)
        //        {
        //            if (gapShell[compacti] != -1)
        //            {
        //                insideStr.write(mesh_.cellCentres()[celli]);
        //            }
        //            else
        //            {
        //                outsideStr.write(mesh_.cellCentres()[celli]);
        //            }
        //        }
        //    }
        //    forAll(bFaceToCompact, bFacei)
        //    {
        //        const label compacti = bFaceToCompact[bFacei];
        //        if (compacti != -1)
        //        {
        //            if (gapShell[compacti] != -1)
        //            {
        //                insideStr.write(neiCc[bFacei]);
        //            }
        //            else
        //            {
        //                outsideStr.write(neiCc[bFacei]);
        //            }
        //        }
        //    }
        //}


        const List<FixedList<label, 3>>& extendedGapLevel =
            surfaces_.extendedGapLevel();
        const List<volumeType>& extendedGapMode =
            surfaces_.extendedGapMode();
        const boolList& extendedGapSelf = surfaces_.gapSelf();

        labelList ccSurface1;
        List<pointIndexHit> ccHit1;
        labelList ccRegion1;
        vectorField ccNormal1;
        {
            labelList ccSurface2;
            List<pointIndexHit> ccHit2;
            labelList ccRegion2;
            vectorField ccNormal2;

            surfaces_.findNearestIntersection
            (
                identity(surfaces_.surfaces().size()),
                start,
                end,

                ccSurface1,
                ccHit1,
                ccRegion1,
                ccNormal1,

                ccSurface2,
                ccHit2,
                ccRegion2,
                ccNormal2
            );
        }

        start.clear();
        end.clear();

        DynamicField<point> rayStart(2*ccSurface1.size());
        DynamicField<point> rayEnd(2*ccSurface1.size());
        DynamicField<scalar> gapSize(2*ccSurface1.size());

        DynamicField<point> rayStart2(2*ccSurface1.size());
        DynamicField<point> rayEnd2(2*ccSurface1.size());
        DynamicField<scalar> gapSize2(2*ccSurface1.size());

        DynamicList<label> cellMap(2*ccSurface1.size());
        DynamicList<label> compactMap(2*ccSurface1.size());

        forAll(ccSurface1, i)
        {
            label surfI = ccSurface1[i];

            if (surfI != -1)
            {
                label globalRegionI =
                    surfaces_.globalRegion(surfI, ccRegion1[i]);

                label faceI = testFaces[i];
                const point& surfPt = ccHit1[i].hitPoint();

                label own = mesh_.faceOwner()[faceI];
                if
                (
                    cellToCompact[own] != -1
                 && shellGapInfo[cellToCompact[own]][2] > 0
                )
                {
                    // Combine info from shell and surface
                    label compactI = cellToCompact[own];
                    FixedList<label, 3> gapInfo;
                    volumeType gapMode;
                    mergeGapInfo
                    (
                        shellGapInfo[compactI],
                        shellGapMode[compactI],
                        extendedGapLevel[globalRegionI],
                        extendedGapMode[globalRegionI],

                        gapInfo,
                        gapMode
                    );

                    const point& cc = cellCentres[own];
                    label nRays = generateRays
                    (
                        false,
                        surfPt,
                        ccNormal1[i],
                        gapInfo,
                        gapMode,
                        surfPt+((cc-surfPt)&ccNormal1[i])*ccNormal1[i],
                        cellLevel[own],

                        rayStart,
                        rayEnd,
                        gapSize,

                        rayStart2,
                        rayEnd2,
                        gapSize2
                    );
                    for (label j = 0; j < nRays; j++)
                    {
                        cellMap.append(own);
                        compactMap.append(i);
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if
                    (
                        cellToCompact[nei] != -1
                     && shellGapInfo[cellToCompact[nei]][2] > 0
                    )
                    {
                        // Combine info from shell and surface
                        label compactI = cellToCompact[nei];
                        FixedList<label, 3> gapInfo;
                        volumeType gapMode;
                        mergeGapInfo
                        (
                            shellGapInfo[compactI],
                            shellGapMode[compactI],
                            extendedGapLevel[globalRegionI],
                            extendedGapMode[globalRegionI],

                            gapInfo,
                            gapMode
                        );

                        const point& cc = cellCentres[nei];
                        label nRays = generateRays
                        (
                            false,
                            surfPt,
                            ccNormal1[i],
                            gapInfo,
                            gapMode,
                            surfPt+((cc-surfPt)&ccNormal1[i])*ccNormal1[i],
                            cellLevel[nei],

                            rayStart,
                            rayEnd,
                            gapSize,

                            rayStart2,
                            rayEnd2,
                            gapSize2
                        );
                        for (label j = 0; j < nRays; j++)
                        {
                            cellMap.append(nei);
                            compactMap.append(i);
                        }
                    }
                }
                else
                {
                    // Note: on coupled face. What cell are we going to
                    // refine? We've got the neighbouring cell centre
                    // and level but we cannot mark it for refinement on
                    // this side...
                    label bFaceI = faceI - mesh_.nInternalFaces();

                    if
                    (
                        bFaceToCompact[bFaceI] != -1
                     && shellGapInfo[bFaceToCompact[bFaceI]][2] > 0
                    )
                    {
                        // Combine info from shell and surface
                        label compactI = bFaceToCompact[bFaceI];
                        FixedList<label, 3> gapInfo;
                        volumeType gapMode;
                        mergeGapInfo
                        (
                            shellGapInfo[compactI],
                            shellGapMode[compactI],
                            extendedGapLevel[globalRegionI],
                            extendedGapMode[globalRegionI],

                            gapInfo,
                            gapMode
                        );

                        const point& cc = neiCc[bFaceI];
                        label nRays = generateRays
                        (
                            false,
                            surfPt,
                            ccNormal1[i],
                            gapInfo,
                            gapMode,
                            surfPt+((cc-surfPt)&ccNormal1[i])*ccNormal1[i],
                            neiLevel[bFaceI],

                            rayStart,
                            rayEnd,
                            gapSize,

                            rayStart2,
                            rayEnd2,
                            gapSize2
                        );
                        for (label j = 0; j < nRays; j++)
                        {
                            cellMap.append(-1); // See above.
                            compactMap.append(i);
                        }
                    }
                }
            }
        }

        Info<< "Shooting " << returnReduce(rayStart.size(), sumOp<label>())
            << " rays from " << returnReduce(testFaces.size(), sumOp<label>())
            << " intersected faces" << endl;

        rayStart.shrink();
        rayEnd.shrink();
        gapSize.shrink();

        rayStart2.shrink();
        rayEnd2.shrink();
        gapSize2.shrink();

        cellMap.shrink();
        compactMap.shrink();

        testFaces.clear();
        ccSurface1.clear();
        ccHit1.clear();
        ccRegion1.clear();
        ccNormal1 = UIndirectList<vector>(ccNormal1, compactMap)();


        // Do intersections in pairs
        labelList surf1;
        List<pointIndexHit> hit1;
        vectorField normal1;
        surfaces_.findNearestIntersection
        (
            rayStart,
            rayEnd,
            surf1,
            hit1,
            normal1
        );

        labelList surf2;
        List<pointIndexHit> hit2;
        vectorField normal2;
        surfaces_.findNearestIntersection
        (
            rayStart2,
            rayEnd2,
            surf2,
            hit2,
            normal2
        );

        forAll(surf1, i)
        {
            // Combine selfProx of shell and surfaces.
            // Ignore regions for now
            const label cellI = cellMap[i];

            const label shelli =
            (
                (cellI != -1 && cellToCompact[cellI] != -1)
              ? gapShell[cellToCompact[cellI]]
              : -1
            );

            bool selfProx = true;
            if (shelli != -1)
            {
                selfProx = shells_.gapSelf()[shelli][0];
            }
            if (surf1[i] != -1 && selfProx)
            {
                const label globalRegioni = surfaces_.globalRegion(surf1[i], 0);
                selfProx = extendedGapSelf[globalRegioni];
            }

            if
            (
                surf1[i] != -1
             && surf2[i] != -1
             && (surf2[i] != surf1[i] || selfProx)
            )
            {
                // Found intersection with surface. Check opposite normal.
                if
                (
                    cellI != -1
                 && (mag(normal1[i]&normal2[i]) > planarCos)
                 && (
                        hit1[i].hitPoint().distSqr(hit2[i].hitPoint())
                      < Foam::sqr(gapSize[i])
                    )
                )
                {
                    if
                    (
                       !markForRefine
                        (
                            surf1[i],
                            nAllowRefine,
                            refineCell[cellI],
                            nRefine
                        )
                    )
                    {
                        break;
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
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


//Foam::meshRefinement::findNearestOppositeOp::findNearestOppositeOp
//(
//    const indexedOctree<treeDataTriSurface>& tree,
//    const point& oppositePoint,
//    const vector& oppositeNormal,
//    const scalar minCos
//)
//:
//    tree_(tree),
//    oppositePoint_(oppositePoint),
//    oppositeNormal_(oppositeNormal),
//    minCos_(minCos)
//{}
//
//
//void Foam::meshRefinement::findNearestOppositeOp::operator()
//(
//    const labelUList& indices,
//    const point& sample,
//    scalar& nearestDistSqr,
//    label& minIndex,
//    point& nearestPoint
//) const
//{
//    const treeDataTriSurface& shape = tree_.shapes();
//    const triSurface& patch = shape.patch();
//    const pointField& points = patch.points();
//
//    forAll(indices, i)
//    {
//        const label index = indices[i];
//        const labelledTri& f = patch[index];
//
//        pointHit nearHit = f.nearestPoint(sample, points);
//        scalar distSqr = sqr(nearHit.distance());
//
//        if (distSqr < nearestDistSqr)
//        {
//            // Nearer. Check if
//            // - a bit way from other hit
//            // - in correct search cone
//            vector d(nearHit.point()-oppositePoint_);
//            scalar normalDist(d&oppositeNormal_);
//
//            if (normalDist > Foam::sqr(SMALL) && normalDist/mag(d) > minCos_)
//            {
//                nearestDistSqr = distSqr;
//                minIndex = index;
//                nearestPoint = nearHit.point();
//            }
//        }
//    }
//}
//
//
//void Foam::meshRefinement::searchCone
//(
//    const label surfI,
//    labelList& nearMap,                 // cells
//    scalarField& nearGap,               // gap size
//    List<pointIndexHit>& nearInfo,      // nearest point on surface
//    List<pointIndexHit>& oppositeInfo   // detected point on gap (or miss)
//) const
//{
//    const labelList& cellLevel = meshCutter_.cellLevel();
//    const pointField& cellCentres = mesh_.cellCentres();
//    const scalar edge0Len = meshCutter_.level0EdgeLength();
//
//    const labelList& surfaceIndices = surfaces_.surfaces();
//    const List<FixedList<label, 3>>& extendedGapLevel =
//        surfaces_.extendedGapLevel();
//    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();
//
//
//    label geomI = surfaceIndices[surfI];
//    const searchableSurface& geom = surfaces_.geometry()[geomI];
//
//    const triSurfaceMesh& s = refCast<const triSurfaceMesh>(geom);
//    const indexedOctree<treeDataTriSurface>& tree = s.tree();
//
//
//    const scalar searchCos = Foam::cos(degToRad(30.0));
//
//    // Normals for ray shooting and inside/outside detection
//    vectorField nearNormal;
//    geom.getNormal(nearInfo, nearNormal);
//    // Regions
//    labelList nearRegion;
//    geom.getRegion(nearInfo, nearRegion);
//
//
//    // Now loop over all near points and search in the half cone
//    labelList map(nearInfo.size());
//    label compactI = 0;
//
//    oppositeInfo.setSize(nearInfo.size());
//
//    forAll(nearInfo, i)
//    {
//        label globalRegionI =
//            surfaces_.globalRegion(surfI, nearRegion[i]);
//
//        // Get updated gap information now we have the region
//        label nGapCells = extendedGapLevel[globalRegionI][0];
//        label minLevel = extendedGapLevel[globalRegionI][1];
//        label maxLevel = extendedGapLevel[globalRegionI][2];
//        volumeType mode = extendedGapMode[globalRegionI];
//
//        label cellI = nearMap[i];
//        label cLevel = cellLevel[cellI];
//
//        if (cLevel >= minLevel && cLevel < maxLevel)
//        {
//            scalar cellSize = edge0Len/pow(2.0, cLevel);
//
//            // Update gap size
//            nearGap[i] = nGapCells*cellSize;
//
//            const point& nearPt = nearInfo[i].hitPoint();
//            vector v(cellCentres[cellI]-nearPt);
//            scalar magV = mag(v);
//
//            // Like with ray shooting we want to
//            // - find triangles up to nearGap away on the wanted side of the
//            //   surface
//            // - find triangles up to 0.5*cellSize away on the unwanted side
//            //   of the surface. This is for cells straddling the surface
//            //   where
//            //   the cell centre might be on the wrong side of the surface
//
//            // Tbd: check that cell centre is inbetween the gap hits
//            // (only if the cell is far enough away)
//
//            scalar posNormalSize = 0.0;
//            scalar negNormalSize = 0.0;
//
//            if (mode == volumeType::OUTSIDE)
//            {
//                posNormalSize = nearGap[i];
//                if (magV < 0.5*cellSize)
//                {
//                    negNormalSize = 0.5*cellSize;
//                }
//            }
//            else if (mode == volumeType::INSIDE)
//            {
//                if (magV < 0.5*cellSize)
//                {
//                    posNormalSize = 0.5*cellSize;
//                }
//                negNormalSize = nearGap[i];
//            }
//            else
//            {
//                posNormalSize = nearGap[i];
//                negNormalSize = nearGap[i];
//            }
//
//            // Test with positive normal
//            oppositeInfo[compactI] = tree.findNearest
//            (
//                nearPt,
//                sqr(posNormalSize),
//                findNearestOppositeOp
//                (
//                    tree,
//                    nearPt,
//                    nearNormal[i],
//                    searchCos
//                )
//            );
//
//            if (oppositeInfo[compactI].hit())
//            {
//                map[compactI++] = i;
//            }
//            else
//            {
//                // Test with negative normal
//                oppositeInfo[compactI] = tree.findNearest
//                (
//                    nearPt,
//                    sqr(negNormalSize),
//                    findNearestOppositeOp
//                    (
//                        tree,
//                        nearPt,
//                        -nearNormal[i],
//                        searchCos
//                    )
//                );
//
//                if (oppositeInfo[compactI].hit())
//                {
//                    map[compactI++] = i;
//                }
//            }
//        }
//    }
//
//    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
//        << " hits on the correct side out of "
//        << returnReduce(map.size(), sumOp<label>()) << endl;
//    map.setSize(compactI);
//    oppositeInfo.setSize(compactI);
//
//    nearMap = labelUIndList(nearMap, map)();
//    nearGap = UIndirectList<scalar>(nearGap, map)();
//    nearInfo = UIndirectList<pointIndexHit>(nearInfo, map)();
//    nearNormal = UIndirectList<vector>(nearNormal, map)();
//
//    // Exclude hits which aren't opposite enough. E.g. you might find
//    // a point on a perpendicular wall - but this does not constitute a gap.
//    vectorField oppositeNormal;
//    geom.getNormal(oppositeInfo, oppositeNormal);
//
//    compactI = 0;
//    forAll(oppositeInfo, i)
//    {
//        if ((nearNormal[i] & oppositeNormal[i]) < -0.707)
//        {
//            map[compactI++] = i;
//        }
//    }
//
//    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
//        << " hits opposite the nearest out of "
//        << returnReduce(map.size(), sumOp<label>()) << endl;
//    map.setSize(compactI);
//
//    nearMap = labelUIndList(nearMap, map)();
//    nearGap = UIndirectList<scalar>(nearGap, map)();
//    nearInfo = UIndirectList<pointIndexHit>(nearInfo, map)();
//    oppositeInfo = UIndirectList<pointIndexHit>(oppositeInfo, map)();
//}


Foam::label Foam::meshRefinement::generateRays
(
    const point& nearPoint,
    const vector& nearNormal,
    const FixedList<label, 3>& gapInfo,
    const volumeType& mode,

    const label cLevel,

    DynamicField<point>& start,
    DynamicField<point>& end
) const
{
    label nOldRays = start.size();

    if (cLevel >= gapInfo[1] && cLevel < gapInfo[2] && gapInfo[0] > 0)
    {
        scalar cellSize = meshCutter_.level0EdgeLength()/pow(2.0, cLevel);

        // Calculate gap size
        scalar nearGap = gapInfo[0]*cellSize;

        const vector& n = nearNormal;

        // Situation 'C' above: cell too close. Use surface
        // -normal and -point to shoot rays

        if (mode == volumeType::OUTSIDE)
        {
            start.append(nearPoint+1e-6*n);
            end.append(nearPoint+nearGap*n);
        }
        else if (mode == volumeType::INSIDE)
        {
            start.append(nearPoint-1e-6*n);
            end.append(nearPoint-nearGap*n);
        }
        else if (mode == volumeType::MIXED)
        {
            start.append(nearPoint+1e-6*n);
            end.append(nearPoint+nearGap*n);

            start.append(nearPoint-1e-6*n);
            end.append(nearPoint-nearGap*n);
        }
    }

    return start.size()-nOldRays;
}


Foam::label Foam::meshRefinement::generateRays
(
    const bool useSurfaceNormal,

    const point& nearPoint,
    const vector& nearNormal,
    const FixedList<label, 3>& gapInfo,
    const volumeType& mode,

    const point& cc,
    const label cLevel,

    DynamicField<point>& start,
    DynamicField<point>& end,
    DynamicField<scalar>& gapSize,

    DynamicField<point>& start2,
    DynamicField<point>& end2,
    DynamicField<scalar>& gapSize2
) const
{
    // We want to handle the following cases:
    // - surface: small gap (marked with 'surface'). gap might be
    //            on inside or outside of surface.
    // - A: cell well inside the gap.
    // - B: cell well outside the gap.
    // - C: cell straddling the gap. cell centre might be inside
    //      or outside
    //
    //       +---+
    //       | B |
    //       +---+
    //
    //            +------+
    //            |      |
    //            |   C  |
    //    --------|------|----surface
    //            +------+
    //
    //        +---+
    //        | A |
    //        +---+
    //
    //
    //    --------------------surface
    //
    // So:
    // - find nearest point on surface
    // - in situation A,B decide if on wanted side of surface
    // - detect if locally a gap (and the cell inside the gap) by
    //   shooting a ray from the point on the surface in the direction
    //   of
    //   - A,B: the cell centre
    //   - C: the surface normal and/or negative surface normal
    //   and see we hit anything
    //
    // Variations of this scheme:
    // - always shoot in the direction of the surface normal. This needs
    //   then an additional check to make sure the cell centre is
    //   somewhere inside the gap
    // - instead of ray shooting use a 'constrained' nearest search
    //   by e.g. looking inside a search cone (implemented in searchCone).
    //   The problem with this constrained nearest is that it still uses
    //   the absolute nearest point on each triangle and only afterwards
    //   checks if it is inside the search cone.


    // Decide which near points are good:
    // - with updated minLevel and maxLevel and nearGap make sure
    //   the cell is still a candidate
    //   NOTE: inside the gap the nearest point on the surface will
    //         be HALF the gap size - otherwise we would have found
    //         a point on the opposite side
    // - if the mode is both sides
    // - or if the hit is inside the current cell (situation 'C',
    //   magV < 0.5cellSize)
    // - or otherwise if on the correct side

    label nOldRays = start.size();

    if (cLevel >= gapInfo[1] && cLevel < gapInfo[2] && gapInfo[0] > 0)
    {
        scalar cellSize = meshCutter_.level0EdgeLength()/pow(2.0, cLevel);

        // Calculate gap size
        scalar nearGap = gapInfo[0]*cellSize;

        // Distance to nearest
        vector v(cc-nearPoint);
        scalar magV = mag(v);

        if (useSurfaceNormal || magV < 0.5*cellSize)
        {
            const vector& n = nearNormal;

            // Situation 'C' above: cell too close. Use surface
            // -normal and -point to shoot rays

            if (mode == volumeType::OUTSIDE)
            {
                start.append(nearPoint+1e-6*n);
                end.append(nearPoint+nearGap*n);
                gapSize.append(nearGap);
                // Second vector so we get pairs of intersections
                start2.append(nearPoint+1e-6*n);
                end2.append(nearPoint-1e-6*n);
                gapSize2.append(gapSize.last());
            }
            else if (mode == volumeType::INSIDE)
            {
                start.append(nearPoint-1e-6*n);
                end.append(nearPoint-nearGap*n);
                gapSize.append(nearGap);
                // Second vector so we get pairs of intersections
                start2.append(nearPoint-1e-6*n);
                end2.append(nearPoint+1e-6*n);
                gapSize2.append(gapSize.last());
            }
            else if (mode == volumeType::MIXED)
            {
                // Do both rays:
                // Outside
                {
                    start.append(nearPoint+1e-6*n);
                    end.append(nearPoint+nearGap*n);
                    gapSize.append(nearGap);
                    // Second vector so we get pairs of intersections
                    start2.append(nearPoint+1e-6*n);
                    end2.append(nearPoint-1e-6*n);
                    gapSize2.append(gapSize.last());
                }
                // Inside
                {
                    start.append(nearPoint-1e-6*n);
                    end.append(nearPoint-nearGap*n);
                    gapSize.append(nearGap);
                    // Second vector so we get pairs of intersections
                    start2.append(nearPoint-1e-6*n);
                    end2.append(nearPoint+1e-6*n);
                    gapSize2.append(gapSize.last());
                }
            }
        }
        else
        {
            // Situation 'A' or 'B' above: cell well away. Test if
            // cell on correct side of surface and shoot ray through
            // cell centre. Note: no need to shoot ray in other
            // direction since we're trying to detect cell inside
            // the gap.

            scalar s = (v&nearNormal);

            if
            (
                (mode == volumeType::MIXED)
             || (mode == volumeType::OUTSIDE && s > SMALL)
             || (mode == volumeType::INSIDE && s < -SMALL)
            )
            {
                //// Use single vector through cell centre
                //vector n(v/(magV+ROOTVSMALL));
                //
                //start.append(cc);
                //end.append(cc+nearGap*n);
                //gapSize.append(nearGap);
                //
                //start2.append(cc);
                //end2.append(cc-nearGap*n);
                //gapSize2.append(nearGap);


                //// Shoot some rays through the cell centre
                //// X-direction:
                //start.append(cc);
                //end.append(cc+nearGap*vector(1, 0, 0));
                //gapSize.append(nearGap);
                //
                //start2.append(cc);
                //end2.append(cc-nearGap*vector(1, 0, 0));
                //gapSize2.append(nearGap);
                //
                //// Y-direction:
                //start.append(cc);
                //end.append(cc+nearGap*vector(0, 1, 0));
                //gapSize.append(nearGap);
                //
                //start2.append(cc);
                //end2.append(cc-nearGap*vector(0, 1, 0));
                //gapSize2.append(nearGap);
                //
                //// Z-direction:
                //start.append(cc);
                //end.append(cc+nearGap*vector(0, 0, 1));
                //gapSize.append(nearGap);
                //
                //start2.append(cc);
                //end2.append(cc-nearGap*vector(0, 0, 1));
                //gapSize2.append(nearGap);


                // 3 axes aligned with normal

                // Use vector through cell centre
                vector n(v/(magV+ROOTVSMALL));

                // Get second vector. Make sure it is sufficiently perpendicular
                vector e2(1, 0, 0);
                scalar s = (e2 & n);
                if (mag(s) < 0.9)
                {
                    e2 -= s*n;
                }
                else
                {
                    e2 = vector(0, 1, 0);
                    e2 -= (e2 & n)*n;
                }
                e2 /= mag(e2);

                // Third vector
                vector e3 = n ^ e2;


                // Rays in first direction
                start.append(cc);
                end.append(cc+nearGap*n);
                gapSize.append(nearGap);

                start2.append(cc);
                end2.append(cc-nearGap*n);
                gapSize2.append(nearGap);

                // Rays in second direction
                start.append(cc);
                end.append(cc+nearGap*e2);
                gapSize.append(nearGap);

                start2.append(cc);
                end2.append(cc-nearGap*e2);
                gapSize2.append(nearGap);

                // Rays in third direction
                start.append(cc);
                end.append(cc+nearGap*e3);
                gapSize.append(nearGap);

                start2.append(cc);
                end2.append(cc-nearGap*e3);
                gapSize2.append(nearGap);
            }
        }
    }

    return start.size()-nOldRays;
}


void Foam::meshRefinement::selectGapCandidates
(
    const labelList& refineCell,
    const label nRefine,

    labelList& cellMap,
    labelList& gapShell,
    List<FixedList<label, 3>>& shellGapInfo,
    List<volumeType>& shellGapMode
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Collect cells to test
    cellMap.setSize(cellLevel.size()-nRefine);
    label compactI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            cellMap[compactI++] = cellI;
        }
    }
    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
        << " unmarked cells out of "
        << mesh_.globalData().nTotalCells() << endl;
    cellMap.setSize(compactI);

    // Do test to see whether cells are inside/outside shell with
    // applicable specification (minLevel <= celllevel < maxLevel)
    shells_.findHigherGapLevel
    (
        pointField(cellCentres, cellMap),
        labelUIndList(cellLevel, cellMap)(),

        gapShell,
        shellGapInfo,
        shellGapMode
    );

    // Compact out hits

    labelList map(shellGapInfo.size());
    compactI = 0;
    forAll(shellGapInfo, i)
    {
        if (shellGapInfo[i][2] > 0)
        {
            map[compactI++] = i;
        }
    }

    Info<< "Selected " << returnReduce(compactI, sumOp<label>())
        << " cells inside gap shells out of "
        << mesh_.globalData().nTotalCells() << endl;

    map.setSize(compactI);
    cellMap = labelUIndList(cellMap, map)();
    gapShell = labelUIndList(gapShell, map)();
    shellGapInfo = UIndirectList<FixedList<label, 3>>(shellGapInfo, map)();
    shellGapMode = UIndirectList<volumeType>(shellGapMode, map)();
}


void Foam::meshRefinement::mergeGapInfo
(
    const FixedList<label, 3>& shellGapInfo,
    const volumeType shellGapMode,
    const FixedList<label, 3>& surfGapInfo,
    const volumeType surfGapMode,

    FixedList<label, 3>& gapInfo,
    volumeType& gapMode
) const
{
    if (surfGapInfo[0] == 0)
    {
        gapInfo = shellGapInfo;
        gapMode = shellGapMode;
    }
    else if (shellGapInfo[0] == 0)
    {
        gapInfo = surfGapInfo;
        gapMode = surfGapMode;
    }
    else
    {
        // Both specify a level. Does surface level win? Or does information
        // need to be merged?

        //gapInfo[0] = max(surfGapInfo[0], shellGapInfo[0]);
        //gapInfo[1] = min(surfGapInfo[1], shellGapInfo[1]);
        //gapInfo[2] = max(surfGapInfo[2], shellGapInfo[2]);
        gapInfo = surfGapInfo;
        gapMode = surfGapMode;
    }
}


Foam::label Foam::meshRefinement::markInternalGapRefinement
(
    const scalar planarCos,
    const bool spreadGapSize,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine,
    labelList& numGapCells,
    scalarField& detectedGapSize
) const
{
    detectedGapSize.setSize(mesh_.nCells());
    detectedGapSize = GREAT;
    numGapCells.setSize(mesh_.nCells());
    numGapCells = -1;

    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();
    const scalar edge0Len = meshCutter_.level0EdgeLength();

    const List<FixedList<label, 3>>& extendedGapLevel =
        surfaces_.extendedGapLevel();
    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();
    const boolList& extendedGapSelf = surfaces_.gapSelf();

    // Get the gap level for the shells
    const labelList maxLevel(shells_.maxGapLevel());

    label oldNRefine = nRefine;

    if (max(maxLevel) > 0)
    {
        // Collect cells to test
        labelList cellMap;
        labelList gapShell;
        List<FixedList<label, 3>> shellGapInfo;
        List<volumeType> shellGapMode;
        selectGapCandidates
        (
            refineCell,
            nRefine,

            cellMap,
            gapShell,
            shellGapInfo,
            shellGapMode
        );

        // Find nearest point and normal on the surfaces
        List<pointIndexHit> nearInfo;
        vectorField nearNormal;
        labelList nearSurface;
        labelList nearRegion;
        {
            // Now we have both the cell-level and the gap size information. Use
            // this to calculate the gap size
            scalarField gapSize(cellMap.size());
            forAll(cellMap, i)
            {
                label cellI = cellMap[i];
                scalar cellSize = edge0Len/pow(2.0, cellLevel[cellI]);
                gapSize[i] = shellGapInfo[i][0]*cellSize;
            }

            surfaces_.findNearestRegion
            (
                identity(surfaces_.surfaces().size()),
                pointField(cellCentres, cellMap),
                sqr(gapSize),
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
                label globalRegionI = surfaces_.globalRegion
                (
                    nearSurface[i],
                    nearRegion[i]
                );

                // Combine info from shell and surface
                FixedList<label, 3> gapInfo;
                volumeType gapMode;
                mergeGapInfo
                (
                    shellGapInfo[i],
                    shellGapMode[i],

                    extendedGapLevel[globalRegionI],
                    extendedGapMode[globalRegionI],

                    gapInfo,
                    gapMode
                );

                // Store wanted number of cells in gap
                label cellI = cellMap[i];
                label cLevel = cellLevel[cellI];
                if (cLevel >= gapInfo[1] && cLevel < gapInfo[2])
                {
                    numGapCells[cellI] = max(numGapCells[cellI], gapInfo[0]);
                }

                // Construct one or more rays to test for oppositeness
                label nRays = generateRays
                (
                    false,
                    nearInfo[i].hitPoint(),
                    nearNormal[i],
                    gapInfo,
                    gapMode,

                    cellCentres[cellI],
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
        shellGapInfo.clear();
        shellGapMode.clear();
        nearInfo.clear();
        nearSurface.clear();
        nearRegion.clear();


        // Do intersections in pairs
        labelList surf1;
        List<pointIndexHit> hit1;
        vectorField normal1;
        surfaces_.findNearestIntersection
        (
            rayStart,
            rayEnd,
            surf1,
            hit1,
            normal1
        );

        labelList surf2;
        List<pointIndexHit> hit2;
        vectorField normal2;
        surfaces_.findNearestIntersection
        (
            rayStart2,
            rayEnd2,
            surf2,
            hit2,
            normal2
        );

        // Extract cell based gap size
        forAll(surf1, i)
        {
            // Combine selfProx of shell and surfaces. Ignore regions for
            // now
            const label shelli = gapShell[map[i]];

            bool selfProx = true;
            if (shelli != -1)
            {
                selfProx = shells_.gapSelf()[shelli][0];
            }
            if (surf1[i] != -1 && selfProx)
            {
                const label globalRegioni = surfaces_.globalRegion(surf1[i], 0);
                selfProx = extendedGapSelf[globalRegioni];
            }

            if
            (
                surf1[i] != -1
             && surf2[i] != -1
             && (surf2[i] != surf1[i] || selfProx)
            )
            {
                // Found intersections with surface. Check for
                // - small gap
                // - coplanar normals

                const label cellI = cellMap[i];

                const scalar d2 =
                    hit1[i].hitPoint().distSqr(hit2[i].hitPoint());

                if
                (
                    cellI != -1
                 && (mag(normal1[i]&normal2[i]) > planarCos)
                 && (d2 < Foam::sqr(gapSize[i]))
                )
                {
                    detectedGapSize[cellI] = min
                    (
                        detectedGapSize[cellI],
                        Foam::sqrt(d2)
                    );
                }
            }
        }

        // Spread it
        if (spreadGapSize)
        {
            scalarField boundaryGapSize;
            syncTools::swapBoundaryCellList
            (
                mesh_,
                detectedGapSize,
                boundaryGapSize
            );

            // Field on cells and faces
            List<transportData> cellData(mesh_.nCells());
            List<transportData> faceData(mesh_.nFaces());

            // Start of walk
            const pointField& faceCentres = mesh_.faceCentres();

            DynamicList<label> frontFaces(mesh_.nFaces());
            DynamicList<transportData> frontData(mesh_.nFaces());
            for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
            {
                label own = mesh_.faceOwner()[faceI];
                label nei = mesh_.faceNeighbour()[faceI];

                scalar minSize = min
                (
                    detectedGapSize[own],
                    detectedGapSize[nei]
                );

                if (minSize < GREAT)
                {
                    frontFaces.append(faceI);
                    frontData.append
                    (
                        transportData
                        (
                            faceCentres[faceI],
                            minSize,
                            0.0
                        )
                    );
                }
            }
            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                label own = mesh_.faceOwner()[faceI];

                scalar minSize = min
                (
                    detectedGapSize[own],
                    boundaryGapSize[faceI-mesh_.nInternalFaces()]
                );

                if (minSize < GREAT)
                {
                    frontFaces.append(faceI);
                    frontData.append
                    (
                        transportData
                        (
                            faceCentres[faceI],
                            minSize,
                            0.0
                        )
                    );
                }
            }

            Info<< "Selected "
                << returnReduce(frontFaces.size(), sumOp<label>())
                << " faces for spreading gap size out of "
                << mesh_.globalData().nTotalFaces() << endl;


            transportData::trackData td(surfaceIndex());

            FaceCellWave<transportData, transportData::trackData> deltaCalc
            (
                mesh_,
                frontFaces,
                frontData,
                faceData,
                cellData,
                mesh_.globalData().nTotalCells()+1,
                td
            );


            forAll(cellMap, i)
            {
                label cellI = cellMap[i];
                if
                (
                    cellI != -1
                 && cellData[cellI].valid(deltaCalc.data())
                 && numGapCells[cellI] != -1
                )
                {
                    // Update transported gap size
                    detectedGapSize[cellI] = min
                    (
                        detectedGapSize[cellI],
                        cellData[cellI].data()
                    );
                }
            }
        }


        // Use it
        forAll(cellMap, i)
        {
            label cellI = cellMap[i];

            if (cellI != -1 && numGapCells[cellI] != -1)
            {
                // Needed gap size
                label cLevel = cellLevel[cellI];
                scalar cellSize =
                    meshCutter_.level0EdgeLength()/pow(2.0, cLevel);
                scalar neededGapSize = numGapCells[cellI]*cellSize;

                if (neededGapSize > detectedGapSize[cellI])
                {
                    if
                    (
                       !markForRefine
                        (
                            123,
                            nAllowRefine,
                            refineCell[cellI],
                            nRefine
                        )
                    )
                    {
                        break;
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
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::label Foam::meshRefinement::markSmallFeatureRefinement
(
    const scalar planarCos,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& surfaceIndices = surfaces_.surfaces();
    const List<FixedList<label, 3>>& extendedGapLevel =
        surfaces_.extendedGapLevel();
    const List<volumeType>& extendedGapMode = surfaces_.extendedGapMode();
    const boolList& extendedGapSelf = surfaces_.gapSelf();

    label oldNRefine = nRefine;

    // Check that we're using any gap refinement
    labelList shellMaxLevel(shells_.maxGapLevel());

    if (max(shellMaxLevel) == 0)
    {
        return 0;
    }

    //- Force calculation of tetBasePt
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();


    forAll(surfaceIndices, surfI)
    {
        label geomI = surfaceIndices[surfI];
        const searchableSurface& geom = surfaces_.geometry()[geomI];


        // Get the element index in a roundabout way. Problem is e.g.
        // distributed surface where local indices differ from global
        // ones (needed for getRegion call)

        pointField ctrs;
        labelList region;
        vectorField normal;
        {
            // Representative local coordinates and bounding sphere
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            List<pointIndexHit> info;
            geom.findNearest(ctrs, radiusSqr, info);

            forAll(info, i)
            {
                if (!info[i].hit())
                {
                    FatalErrorInFunction
                        << "fc:" << ctrs[i]
                        << " radius:" << radiusSqr[i]
                        << exit(FatalError);
                }
            }

            geom.getRegion(info, region);
            geom.getNormal(info, normal);
        }

        // Do test to see whether triangles are inside/outside shell with
        // applicable specification (minLevel <= celllevel < maxLevel)
        List<FixedList<label, 3>> shellGapInfo;
        List<volumeType> shellGapMode;
        labelList gapShell;
        shells_.findHigherGapLevel
        (
            ctrs,
            labelList(ctrs.size(), Zero),

            gapShell,
            shellGapInfo,
            shellGapMode
        );


        DynamicList<label> map(ctrs.size());
        DynamicList<label> cellMap(ctrs.size());

        DynamicField<point> rayStart(ctrs.size());
        DynamicField<point> rayEnd(ctrs.size());
        DynamicField<scalar> gapSize(ctrs.size());

        label nTestCells = 0;

        forAll(ctrs, i)
        {
            if (shellGapInfo[i][2] > 0)
            {
                label globalRegionI = surfaces_.globalRegion(surfI, region[i]);

                // Combine info from shell and surface
                FixedList<label, 3> gapInfo;
                volumeType gapMode;
                mergeGapInfo
                (
                    shellGapInfo[i],
                    shellGapMode[i],

                    extendedGapLevel[globalRegionI],
                    extendedGapMode[globalRegionI],

                    gapInfo,
                    gapMode
                );

                //- Option 1: use octree nearest searching inside polyMesh
                //label cellI = mesh_.findCell(pt, polyMesh::CELL_TETS);

                //- Option 2: use octree 'inside' searching inside polyMesh. Is
                //            much faster.
                label cellI = -1;
                const indexedOctree<treeDataCell>& tree = mesh_.cellTree();
                if (tree.nodes().size() && tree.bb().contains(ctrs[i]))
                {
                    cellI = tree.findInside(ctrs[i]);
                }

                if (cellI != -1 && refineCell[cellI] == -1)
                {
                    // Construct one or two rays to test for oppositeness
                    // Note that we always want to use the surface normal
                    // and not the vector from cell centre to surface point

                    label nRays = generateRays
                    (
                        ctrs[i],
                        normal[i],
                        gapInfo,
                        gapMode,

                        cellLevel[cellI],

                        rayStart,
                        rayEnd
                    );

                    if (nRays > 0)
                    {
                        nTestCells++;
                        for (label j = 0; j < nRays; j++)
                        {
                            cellMap.append(cellI);
                            map.append(i);
                        }
                    }
                }
            }
        }

        Info<< "Selected " << returnReduce(nTestCells, sumOp<label>())
            << " cells containing triangle centres out of "
            << mesh_.globalData().nTotalCells() << endl;
        map.shrink();
        cellMap.shrink();
        rayStart.shrink();
        rayEnd.shrink();

        ctrs.clear();
        region.clear();
        shellGapInfo.clear();
        shellGapMode.clear();
        normal = UIndirectList<vector>(normal, map)();

        // Do intersections.
        labelList surfaceHit;
        vectorField surfaceNormal;
        surfaces_.findNearestIntersection
        (
            rayStart,
            rayEnd,
            surfaceHit,
            surfaceNormal
        );


        label nOldRefine = 0;


        forAll(surfaceHit, i)
        {
            // Combine selfProx of shell and surfaces. Ignore regions for
            // now
            const label shelli = gapShell[map[i]];
            bool selfProx = true;
            if (shelli != -1)
            {
                selfProx = shells_.gapSelf()[shelli][0];
            }
            if (surfI != -1 && selfProx)
            {
                const label globalRegioni = surfaces_.globalRegion(surfI, 0);
                selfProx = extendedGapSelf[globalRegioni];
            }

            if
            (
                surfaceHit[i] != -1
             && (surfaceHit[i] != surfI || selfProx)
            )
            {
                // Found intersection with surface. Check coplanar normals.
                label cellI = cellMap[i];

                if (mag(normal[i]&surfaceNormal[i]) > planarCos)
                {
                    if
                    (
                       !markForRefine
                        (
                            surfaceHit[i],
                            nAllowRefine,
                            refineCell[cellI],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }

        Info<< "For surface " << geom.name() << " found "
            << returnReduce(nRefine-nOldRefine, sumOp<label>())
            << " cells in small gaps" << endl;

        if
        (
            returnReduce(nRefine, sumOp<label>())
          > returnReduce(nAllowRefine, sumOp<label>())
        )
        {
            Info<< "Reached refinement limit." << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


////XXXXXXXX
Foam::label Foam::meshRefinement::markSurfaceFieldRefinement
(
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& surfaceIndices = surfaces_.surfaces();

    label oldNRefine = nRefine;

    //- Force calculation of tetBasePt
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();
    const indexedOctree<treeDataCell>& tree = mesh_.cellTree();


    forAll(surfaceIndices, surfI)
    {
        label geomI = surfaceIndices[surfI];
        const searchableSurface& geom = surfaces_.geometry()[geomI];

        // Get the element index in a roundabout way. Problem is e.g.
        // distributed surface where local indices differ from global
        // ones (needed for getRegion call)

        pointField ctrs;
        labelList region;
        labelList minLevelField;
        {
            // Representative local coordinates and bounding sphere
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            List<pointIndexHit> info;
            geom.findNearest(ctrs, radiusSqr, info);

            forAll(info, i)
            {
                if (!info[i].hit())
                {
                    FatalErrorInFunction
                        << "fc:" << ctrs[i]
                        << " radius:" << radiusSqr[i]
                        << exit(FatalError);
                }
            }

            geom.getRegion(info, region);
            geom.getField(info, minLevelField);
        }

        if (minLevelField.size() != geom.size())
        {
            Pout<< "** no minLevelField" << endl;
            continue;
        }


        label nOldRefine = 0;

        forAll(ctrs, i)
        {
            label cellI = -1;
            if (tree.nodes().size() && tree.bb().contains(ctrs[i]))
            {
                cellI = tree.findInside(ctrs[i]);
            }

            if
            (
                cellI != -1
             && refineCell[cellI] == -1
             && minLevelField[i] > cellLevel[cellI]
            )
            {
                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[cellI],
                        nRefine
                    )
                )
                {
                    break;
                }
            }
        }

        Info<< "For surface " << geom.name() << " found "
            << returnReduce(nRefine-nOldRefine, sumOp<label>())
            << " cells containing cached refinement field" << endl;

        if
        (
            returnReduce(nRefine, sumOp<label>())
          > returnReduce(nAllowRefine, sumOp<label>())
        )
        {
            Info<< "Reached refinement limit." << endl;
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}
////XXXXXXXXX


// ************************************************************************* //
