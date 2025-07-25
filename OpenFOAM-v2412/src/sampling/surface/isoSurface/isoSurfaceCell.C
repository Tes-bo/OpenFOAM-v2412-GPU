/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "isoSurfaceCell.H"
#include "isoSurfacePoint.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "mergePoints.H"
#include "tetMatcher.H"
#include "syncTools.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "Time.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isoSurfaceBaseMethods.H"
defineIsoSurfaceInterpolateMethods(Foam::isoSurfaceCell);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceCell, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::isoSurfaceCell::isoFraction
(
    const scalar s0,
    const scalar s1
) const
{
    const scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        return (iso_-s0)/d;
    }

    return -1.0;
}


Foam::labelPair Foam::isoSurfaceCell::findCommonPoints
(
    const labelledTri& tri0,
    const labelledTri& tri1
)
{
    labelPair common(-1, -1);

    label fp0 = 0;
    label fp1 = tri1.find(tri0[fp0]);

    if (fp1 == -1)
    {
        fp0 = 1;
        fp1 = tri1.find(tri0[fp0]);
    }

    if (fp1 != -1)
    {
        // So tri0[fp0] is tri1[fp1]

        // Find next common point
        label fp0p1 = tri0.fcIndex(fp0);
        label fp1p1 = tri1.fcIndex(fp1);
        label fp1m1 = tri1.rcIndex(fp1);

        if (tri0[fp0p1] == tri1[fp1p1] || tri0[fp0p1] == tri1[fp1m1])
        {
            common[0] = tri0[fp0];
            common[1] = tri0[fp0p1];
        }
    }
    return common;
}


Foam::point Foam::isoSurfaceCell::calcCentre(const triSurface& s)
{
    vector sum = Zero;

    forAll(s, i)
    {
        sum += s[i].centre(s.points());
    }
    return sum/s.size();
}


Foam::pointIndexHit Foam::isoSurfaceCell::collapseSurface
(
    const label celli,
    pointField& localPoints,
    DynamicList<labelledTri, 64>& localTris
) const
{
    pointIndexHit info(false, Zero, localTris.size());

    if (localTris.size() == 1)
    {
        const labelledTri& tri = localTris[0];
        info.hitPoint(tri.centre(localPoints));
    }
    else if (localTris.size() == 2)
    {
        // Check if the two triangles share an edge.
        const labelledTri& tri0 = localTris[0];
        const labelledTri& tri1 = localTris[1];

        labelPair shared = findCommonPoints(tri0, tri1);

        if (shared[0] != -1)
        {
            const vector n0 = tri0.areaNormal(localPoints);
            const vector n1 = tri1.areaNormal(localPoints);

            // Merge any zero-sized triangles,
            // or if they point in the same direction.

            if
            (
                mag(n0) <= ROOTVSMALL
             || mag(n1) <= ROOTVSMALL
             || (n0 & n1) >= 0
            )
            {
                info.hitPoint
                (
                    0.5
                  * (
                        tri0.centre(localPoints)
                      + tri1.centre(localPoints)
                    )
                );
            }
        }
    }
    else if (localTris.size())
    {
        // Check if single region. Rare situation.
        triSurface surf
        (
            localTris,
            geometricSurfacePatchList(),
            localPoints,
            true
        );
        localTris.clearStorage();

        labelList faceZone;
        label nZones = surf.markZones
        (
            boolList(surf.nEdges(), false),
            faceZone
        );

        if (nZones == 1)
        {
            // Check that all normals make a decent angle
            scalar minCos = GREAT;
            const vector& n0 = surf.faceNormals()[0];
            for (label i = 1; i < surf.size(); ++i)
            {
                scalar cosAngle = (n0 & surf.faceNormals()[i]);
                if (cosAngle < minCos)
                {
                    minCos = cosAngle;
                }
            }

            if (minCos > 0)
            {
                info.hitPoint(calcCentre(surf));
            }
        }
    }

    return info;
}


void Foam::isoSurfaceCell::calcSnappedCc
(
    const bitSet& isTet,
    const scalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedCc
) const
{
    const pointField& cc = mesh_.cellCentres();
    const pointField& pts = mesh_.points();

    snappedCc.setSize(mesh_.nCells());
    snappedCc = -1;

    // Work arrays
    DynamicList<point, 64> localPoints(64);
    DynamicList<labelledTri, 64> localTris(64);
    Map<label> pointToLocal(64);

    forAll(mesh_.cells(), celli)
    {
        if (cellCutType_[celli] == cutType::CUT) // No tet cuts
        {
            const scalar cVal = cVals[celli];

            const cell& cFaces = mesh_.cells()[celli];

            localPoints.clear();
            localTris.clear();
            pointToLocal.clear();

            // Create points for all intersections close to cell centre
            // (i.e. from pyramid edges)

            for (const label facei : cFaces)
            {
                const face& f = mesh_.faces()[facei];

                for (const label pointi : f)
                {
                    scalar s = isoFraction(cVal, pVals[pointi]);

                    if (s >= 0.0 && s <= 0.5)
                    {
                        if (pointToLocal.insert(pointi, localPoints.size()))
                        {
                            localPoints.append((1.0-s)*cc[celli]+s*pts[pointi]);
                        }
                    }
                }
            }

            if (localPoints.size() == 1)
            {
                // No need for any analysis.
                snappedCc[celli] = snappedPoints.size();
                snappedPoints.append(localPoints[0]);

                //Pout<< "cell:" << celli
                //    << " at " << mesh_.cellCentres()[celli]
                //    << " collapsing " << localPoints
                //    << " intersections down to "
                //    << snappedPoints[snappedCc[celli]] << endl;
            }
            else if (localPoints.size() == 2)
            {
                //? No need for any analysis.???
                snappedCc[celli] = snappedPoints.size();
                snappedPoints.append(0.5*(localPoints[0]+localPoints[1]));

                //Pout<< "cell:" << celli
                //    << " at " << mesh_.cellCentres()[celli]
                //    << " collapsing " << localPoints
                //    << " intersections down to "
                //    << snappedPoints[snappedCc[celli]] << endl;
            }
            else if (localPoints.size())
            {
                // Need to analyse
                for (const label facei : cFaces)
                {
                    const face& f = mesh_.faces()[facei];

                    // Do a tetrahedralisation. Each face to cc becomes pyr.
                    // Each pyr gets split into tets by diagonalisation
                    // of face.

                    const label fp0 = mesh_.tetBasePtIs()[facei];
                    label fp = f.fcIndex(fp0);
                    for (label i = 2; i < f.size(); ++i)
                    {
                        label nextFp = f.fcIndex(fp);
                        triFace tri(f[fp0], f[fp], f[nextFp]);

                        // Get fractions for the three edges to cell centre
                        FixedList<scalar, 3> s(3);
                        s[0] = isoFraction(cVal, pVals[tri[0]]);
                        s[1] = isoFraction(cVal, pVals[tri[1]]);
                        s[2] = isoFraction(cVal, pVals[tri[2]]);

                        if
                        (
                            (s[0] >= 0.0 && s[0] <= 0.5)
                         && (s[1] >= 0.0 && s[1] <= 0.5)
                         && (s[2] >= 0.0 && s[2] <= 0.5)
                        )
                        {
                            if
                            (
                                (mesh_.faceOwner()[facei] == celli)
                             == (cVal >= pVals[tri[0]])
                            )
                            {
                                localTris.append
                                (
                                    labelledTri
                                    (
                                        pointToLocal[tri[1]],
                                        pointToLocal[tri[0]],
                                        pointToLocal[tri[2]],
                                        0
                                    )
                                );
                            }
                            else
                            {
                                localTris.append
                                (
                                    labelledTri
                                    (
                                        pointToLocal[tri[0]],
                                        pointToLocal[tri[1]],
                                        pointToLocal[tri[2]],
                                        0
                                    )
                                );
                            }
                        }

                        fp = nextFp;
                    }
                }

                pointField surfPoints;
                surfPoints.transfer(localPoints);
                pointIndexHit info = collapseSurface
                (
                    celli,
                    surfPoints,
                    localTris
                );

                if (info.hit())
                {
                    snappedCc[celli] = snappedPoints.size();
                    snappedPoints.append(info.point());

                    //Pout<< "cell:" << celli
                    //    << " at " << mesh_.cellCentres()[celli]
                    //    << " collapsing " << surfPoints
                    //    << " intersections down to "
                    //    << snappedPoints[snappedCc[celli]] << endl;
                }
            }
        }
    }
}


void Foam::isoSurfaceCell::genPointTris
(
    const scalarField& cellValues,
    const scalarField& pointValues,
    const label pointi,
    const label facei,
    const label celli,
    DynamicList<point, 64>& localTriPoints
) const
{
    const pointField& cc = mesh_.cellCentres();
    const pointField& pts = mesh_.points();
    const face& f = mesh_.faces()[facei];

    const label fp0 = mesh_.tetBasePtIs()[facei];
    label fp = f.fcIndex(fp0);
    for (label i = 2; i < f.size(); ++i)
    {
        label nextFp = f.fcIndex(fp);
        triFace tri(f[fp0], f[fp], f[nextFp]);

        label index = tri.find(pointi);

        if (index == -1)
        {
            continue;
        }

        // Tet between index..index-1, index..index+1, index..cc
        label b = tri[tri.fcIndex(index)];
        label c = tri[tri.rcIndex(index)];

        // Get fractions for the three edges emanating from point
        FixedList<scalar, 3> s(3);
        s[0] = isoFraction(pointValues[pointi], pointValues[b]);
        s[1] = isoFraction(pointValues[pointi], pointValues[c]);
        s[2] = isoFraction(pointValues[pointi], cellValues[celli]);

        if
        (
            (s[0] >= 0.0 && s[0] <= 0.5)
         && (s[1] >= 0.0 && s[1] <= 0.5)
         && (s[2] >= 0.0 && s[2] <= 0.5)
        )
        {
            point p0 = (1.0-s[0])*pts[pointi] + s[0]*pts[b];
            point p1 = (1.0-s[1])*pts[pointi] + s[1]*pts[c];
            point p2 = (1.0-s[2])*pts[pointi] + s[2]*cc[celli];

            if
            (
                (mesh_.faceOwner()[facei] == celli)
             == (pointValues[pointi] > cellValues[celli])
            )
            {
                localTriPoints.append(p0);
                localTriPoints.append(p1);
                localTriPoints.append(p2);
            }
            else
            {
                localTriPoints.append(p1);
                localTriPoints.append(p0);
                localTriPoints.append(p2);
            }
        }

        fp = nextFp;
    }
}


void Foam::isoSurfaceCell::genPointTris
(
    const scalarField& pointValues,
    const label pointi,
    const label facei,
    const label celli,
    DynamicList<point, 64>& localTriPoints
) const
{
    const pointField& pts = mesh_.points();
    const cell& cFaces = mesh_.cells()[celli];

    // Make tet from this face to the 4th point (same as cellcentre in
    // non-tet cells)
    const face& f = mesh_.faces()[facei];

    // Find 4th point
    label ccPointi = -1;
    for (const label cfacei : cFaces)
    {
        const face& f1 = mesh_.faces()[cfacei];
        for (const label p1 : f1)
        {
            if (!f.found(p1))
            {
                ccPointi = p1;
                break;
            }
        }
        if (ccPointi != -1)
        {
            break;
        }
    }


    // Tet between index..index-1, index..index+1, index..cc
    label index = f.find(pointi);
    label b = f[f.fcIndex(index)];
    label c = f[f.rcIndex(index)];

    //Pout<< " p0:" << pointi << " b:" << b << " c:" << c
    //<< " d:" << ccPointi << endl;

    // Get fractions for the three edges emanating from point
    FixedList<scalar, 3> s(3);
    s[0] = isoFraction(pointValues[pointi], pointValues[b]);
    s[1] = isoFraction(pointValues[pointi], pointValues[c]);
    s[2] = isoFraction(pointValues[pointi], pointValues[ccPointi]);

    if
    (
        (s[0] >= 0.0 && s[0] <= 0.5)
     && (s[1] >= 0.0 && s[1] <= 0.5)
     && (s[2] >= 0.0 && s[2] <= 0.5)
    )
    {
        point p0 = (1.0-s[0])*pts[pointi] + s[0]*pts[b];
        point p1 = (1.0-s[1])*pts[pointi] + s[1]*pts[c];
        point p2 = (1.0-s[2])*pts[pointi] + s[2]*pts[ccPointi];

        if (mesh_.faceOwner()[facei] != celli)
        {
            localTriPoints.append(p0);
            localTriPoints.append(p1);
            localTriPoints.append(p2);
        }
        else
        {
            localTriPoints.append(p1);
            localTriPoints.append(p0);
            localTriPoints.append(p2);
        }
    }
}


void Foam::isoSurfaceCell::calcSnappedPoint
(
    const bitSet& isTet,
    const scalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedPoint
) const
{
    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();

    // Determine if point is on boundary. Points on boundaries are never
    // snapped. Coupled boundaries are handled explicitly so not marked here.
    bitSet isBoundaryPoint(mesh_.nPoints());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if (!pp.coupled())
        {
            for (const label facei : pp.range())
            {
                const face& f = mesh_.faces()[facei];

                isBoundaryPoint.set(f);
            }
        }
    }


    const point greatPoint(GREAT, GREAT, GREAT);

    pointField collapsedPoint(mesh_.nPoints(), greatPoint);


    // Work arrays
    DynamicList<point, 64> localTriPoints(100);
    labelHashSet localPointCells(100);

    forAll(mesh_.pointFaces(), pointi)
    {
        constexpr uint8_t realCut(cutType::CUT | cutType::TETCUT);

        if (isBoundaryPoint.test(pointi))
        {
            continue;
        }

        const labelList& pFaces = mesh_.pointFaces()[pointi];

        bool anyCut = false;

        for (const label facei : pFaces)
        {
            if
            (
                (cellCutType_[faceOwn[facei]] & realCut) != 0
             || (
                    mesh_.isInternalFace(facei)
                 && (cellCutType_[faceNei[facei]] & realCut) != 0
                )
            )
            {
                anyCut = true;
                break;
            }
        }

        if (!anyCut)
        {
            continue;
        }


        // Do a pointCells walk (by using pointFaces)

        localPointCells.clear();
        localTriPoints.clear();

        for (const label facei : pFaces)
        {
            const label own = faceOwn[facei];

            if (isTet.test(own))
            {
                // Since tets have no cell centre to include make sure
                // we only generate a triangle once per point.
                if (localPointCells.insert(own))
                {
                    genPointTris(pVals, pointi, facei, own, localTriPoints);
                }
            }
            else
            {
                genPointTris
                (
                    cVals,
                    pVals,
                    pointi,
                    facei,
                    own,
                    localTriPoints
                );
            }

            if (mesh_.isInternalFace(facei))
            {
                const label nei = faceNei[facei];

                if (isTet.test(nei))
                {
                    if (localPointCells.insert(nei))
                    {
                        genPointTris(pVals, pointi, facei, nei, localTriPoints);
                    }
                }
                else
                {
                    genPointTris
                    (
                        cVals,
                        pVals,
                        pointi,
                        facei,
                        nei,
                        localTriPoints
                    );
                }
            }
        }

        if (localTriPoints.size() == 3)
        {
            // Single triangle. No need for any analysis. Average points.
            pointField points;
            points.transfer(localTriPoints);
            collapsedPoint[pointi] = sum(points)/points.size();

            //Pout<< "    point:" << pointi
            //    << " replacing coord:" << mesh_.points()[pointi]
            //    << " by average:" << collapsedPoint[pointi] << endl;
        }
        else if (localTriPoints.size())
        {
            // Convert points into triSurface.

            // Merge points and compact out non-valid triangles
            labelList triMap;               // merged to unmerged triangle
            labelList triPointReverseMap;   // unmerged to merged point
            triSurface surf
            (
                stitchTriPoints
                (
                    false,                  // do not check for duplicate tris
                    localTriPoints,
                    triPointReverseMap,
                    triMap
                )
            );

            labelList faceZone;
            label nZones = surf.markZones
            (
                boolList(surf.nEdges(), false),
                faceZone
            );

            if (nZones == 1)
            {
                // Check that all normals make a decent angle
                scalar minCos = GREAT;
                const vector& n0 = surf.faceNormals()[0];
                for (label i = 1; i < surf.size(); ++i)
                {
                    const vector& n = surf.faceNormals()[i];
                    scalar cosAngle = (n0 & n);
                    if (cosAngle < minCos)
                    {
                        minCos = cosAngle;
                    }
                }
                if (minCos > 0)
                {
                    collapsedPoint[pointi] = calcCentre(surf);
                }
            }
        }
    }

    syncTools::syncPointPositions
    (
        mesh_,
        collapsedPoint,
        minMagSqrEqOp<point>(),
        greatPoint
    );

    snappedPoint.setSize(mesh_.nPoints());
    snappedPoint = -1;

    forAll(collapsedPoint, pointi)
    {
        // Cannot do == comparison since might be transformed so have
        // truncation errors.
        if (magSqr(collapsedPoint[pointi]) < 0.5*magSqr(greatPoint))
        {
            snappedPoint[pointi] = snappedPoints.size();
            snappedPoints.append(collapsedPoint[pointi]);
        }
    }
}


Foam::triSurface Foam::isoSurfaceCell::stitchTriPoints
(
    const bool checkDuplicates,
    const List<point>& triPoints,
    labelList& triPointReverseMap,  // unmerged to merged point
    labelList& triMap               // merged to unmerged triangle
) const
{
    label nTris = triPoints.size()/3;

    if ((triPoints.size() % 3) != 0)
    {
        FatalErrorInFunction
            << "Problem: number of points " << triPoints.size()
            << " not a multiple of 3." << abort(FatalError);
    }

    pointField newPoints;
    Foam::mergePoints
    (
        triPoints,
        mergeDistance_,
        false,
        triPointReverseMap,
        newPoints
    );

    // Check that enough merged.
    if (debug)
    {
        Pout<< "isoSurfaceCell : merged from " << triPointReverseMap.size()
            << " points down to " << newPoints.size() << endl;

        labelList oldToNew;
        const label nUnique = Foam::mergePoints
        (
            newPoints,
            mergeDistance_,
            true,
            oldToNew
        );

        if (nUnique != newPoints.size())
        {
            FatalErrorInFunction
                << "Merged points contain duplicates"
                << " when merging with distance " << mergeDistance_ << endl
                << "merged:" << newPoints.size() << " re-merged:"
                << nUnique
                << abort(FatalError);
        }
    }


    List<labelledTri> tris;
    {
        DynamicList<labelledTri> dynTris(nTris);
        label rawPointi = 0;
        DynamicList<label> newToOldTri(nTris);

        for (label oldTriI = 0; oldTriI < nTris; ++oldTriI)
        {
            labelledTri tri
            (
                triPointReverseMap[rawPointi],
                triPointReverseMap[rawPointi+1],
                triPointReverseMap[rawPointi+2],
                0
            );
            if (tri.good())
            {
                newToOldTri.append(oldTriI);
                dynTris.append(tri);
            }

            rawPointi += 3;
        }

        triMap.transfer(newToOldTri);
        tris.transfer(dynTris);
    }


    // Use face centres to determine 'flat hole' situation (see RMT paper).
    // Two unconnected triangles get connected because (some of) the edges
    // separating them get collapsed. Below only checks for duplicate triangles,
    // not non-manifold edge connectivity.
    if (checkDuplicates)
    {
        if (debug)
        {
            Pout<< "isoSurfaceCell : merged from " << nTris
                << " down to " << tris.size() << " triangles." << endl;
        }

        pointField centres(tris.size());
        forAll(tris, triI)
        {
            centres[triI] = tris[triI].centre(newPoints);
        }

        labelList oldToMerged;
        label nUnique = Foam::mergePoints
        (
            centres,
            mergeDistance_,
            false,
            oldToMerged
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : detected "
                << (oldToMerged.size() - nUnique)
                << " duplicate triangles." << endl;
        }

        if (oldToMerged.size() != nUnique)
        {
            // Filter out duplicates.
            label newTriI = 0;
            DynamicList<label> newToOldTri(tris.size());
            labelList newToMaster(nUnique, -1);
            forAll(tris, triI)
            {
                label mergedI = oldToMerged[triI];

                if (newToMaster[mergedI] == -1)
                {
                    newToMaster[mergedI] = triI;
                    newToOldTri.append(triMap[triI]);
                    tris[newTriI++] = tris[triI];
                }
            }

            triMap.transfer(newToOldTri);
            tris.setSize(newTriI);
        }
    }

    return triSurface(tris, geometricSurfacePatchList(), newPoints, true);
}


void Foam::isoSurfaceCell::calcAddressing
(
    const triSurface& surf,
    List<FixedList<label, 3>>& faceEdges,
    labelList& edgeFace0,
    labelList& edgeFace1,
    Map<labelList>& edgeFacesRest
) const
{
    const pointField& points = surf.points();

    pointField edgeCentres(3*surf.size());
    label edgeI = 0;
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];
        edgeCentres[edgeI++] = 0.5*(points[tri[0]]+points[tri[1]]);
        edgeCentres[edgeI++] = 0.5*(points[tri[1]]+points[tri[2]]);
        edgeCentres[edgeI++] = 0.5*(points[tri[2]]+points[tri[0]]);
    }

    labelList oldToMerged;
    const label nUnique = Foam::mergePoints
    (
        edgeCentres,
        mergeDistance_,
        false,
        oldToMerged
    );

    if (debug)
    {
        Pout<< "isoSurfaceCell : detected "
            << nUnique
            << " edges on " << surf.size() << " triangles." << endl;
    }

    if (nUnique == edgeCentres.size())
    {
        // Nothing to do
        return;
    }


    // Determine faceEdges
    faceEdges.setSize(surf.size());
    edgeI = 0;
    forAll(surf, triI)
    {
        faceEdges[triI][0] = oldToMerged[edgeI++];
        faceEdges[triI][1] = oldToMerged[edgeI++];
        faceEdges[triI][2] = oldToMerged[edgeI++];
    }


    // Determine edgeFaces
    edgeFace0.resize_nocopy(nUnique);
    edgeFace0 = -1;
    edgeFace1.resize_nocopy(nUnique);
    edgeFace1 = -1;
    edgeFacesRest.clear();

    forAll(oldToMerged, oldEdgeI)
    {
        label triI = oldEdgeI / 3;
        label edgeI = oldToMerged[oldEdgeI];

        if (edgeFace0[edgeI] == -1)
        {
            edgeFace0[edgeI] = triI;
        }
        else if (edgeFace1[edgeI] == -1)
        {
            edgeFace1[edgeI] = triI;
        }
        else
        {
            //WarningInFunction
            //    << "Edge " << edgeI << " with centre "
            //    << " used by more than two triangles: " << edgeFace0[edgeI]
            //    << ", "
            //    << edgeFace1[edgeI] << " and " << triI << endl;

            edgeFacesRest(edgeI).push_back(triI);
        }
    }
}


bool Foam::isoSurfaceCell::danglingTriangle
(
    const FixedList<label, 3>& fEdges,
    const labelList& edgeFace1
)
{
    label nOpen = 0;
    for (const label edgei : fEdges)
    {
        if (edgeFace1[edgei] == -1)
        {
            ++nOpen;
        }
    }

    return (nOpen == 1 || nOpen == 2 || nOpen == 3);
}


Foam::label Foam::isoSurfaceCell::markDanglingTriangles
(
    const List<FixedList<label, 3>>& faceEdges,
    const labelList& edgeFace0,
    const labelList& edgeFace1,
    const Map<labelList>& edgeFacesRest,
    boolList& keepTriangles
)
{
    keepTriangles.setSize(faceEdges.size());
    keepTriangles = true;

    label nDangling = 0;

    // Remove any dangling triangles
    forAllConstIters(edgeFacesRest, iter)
    {
        // These are all the non-manifold edges. Filter out all triangles
        // with only one connected edge (= this edge)

        const label edgeI = iter.key();
        const labelList& otherEdgeFaces = iter.val();

        // Remove all dangling triangles
        if (danglingTriangle(faceEdges[edgeFace0[edgeI]], edgeFace1))
        {
            keepTriangles[edgeFace0[edgeI]] = false;
            ++nDangling;
        }
        if (danglingTriangle(faceEdges[edgeFace1[edgeI]], edgeFace1))
        {
            keepTriangles[edgeFace1[edgeI]] = false;
            ++nDangling;
        }
        for (const label triI : otherEdgeFaces)
        {
            if (danglingTriangle(faceEdges[triI], edgeFace1))
            {
                keepTriangles[triI] = false;
                ++nDangling;
            }
        }
    }
    return nDangling;
}


Foam::triSurface Foam::isoSurfaceCell::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldFaces,
    labelList& oldToNewPoints,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        ListOps::createWithValue<bool>(s.size(), newToOldFaces, true, false)
    );

    newToOldPoints.setSize(s.points().size());
    oldToNewPoints.setSize(s.points().size());
    oldToNewPoints = -1;
    {
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Renumber labels for face
                for (const label oldPointi : s[oldFacei])
                {
                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }

    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = s.points()[newToOldPoints[i]];
    }
    // Extract faces
    List<labelledTri> newTriangles(newToOldFaces.size());

    forAll(newToOldFaces, i)
    {
        // Get old vertex labels
        const labelledTri& tri = s[newToOldFaces[i]];

        newTriangles[i][0] = oldToNewPoints[tri[0]];
        newTriangles[i][1] = oldToNewPoints[tri[1]];
        newTriangles[i][2] = oldToNewPoints[tri[2]];
        newTriangles[i].region() = tri.region();
    }

    // Reuse storage.
    return triSurface(newTriangles, s.patches(), newPoints, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceCell::isoSurfaceCell
(
    const polyMesh& mesh,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const isoSurfaceParams& params,
    const bitSet& ignoreCells
)
:
    isoSurfaceBase(mesh, cellValues, pointValues, iso, params),
    mergeDistance_(params.mergeTol()*mesh.bounds().mag()),
    cellCutType_(mesh.nCells(), cutType::UNVISITED)
{
    const bool regularise = (params.filter() != filterType::NONE);

    if (debug)
    {
        Pout<< "isoSurfaceCell:" << nl
            << "    cell min/max  : " << minMax(cVals_) << nl
            << "    point min/max : " << minMax(pVals_) << nl
            << "    isoValue      : " << iso << nl
            << "    filter        : " << Switch(regularise) << nl
            << "    mergeTol      : " << params.mergeTol() << nl
            << "    mesh span     : " << mesh.bounds().mag() << nl
            << "    mergeDistance : " << mergeDistance_ << nl
            << "    ignoreCells   : " << ignoreCells.count()
            << " / " << cVals_.size() << nl
            << endl;
    }


    label nBlockedCells = 0;

    // Mark ignoreCells as blocked
    nBlockedCells += blockCells(cellCutType_, ignoreCells);


    // Some processor domains may require tetBasePtIs and others do not.
    // Do now to ensure things stay synchronized.
    (void)mesh_.tetBasePtIs();


    // Calculate a tet/non-tet filter
    bitSet isTet(mesh_.nCells());
    {
        for (label celli = 0; celli < mesh_.nCells(); ++celli)
        {
            if (tetMatcher::test(mesh_, celli))
            {
                isTet.set(celli);
            }
        }
    }

    // Determine cell cuts
    nCutCells_ = calcCellCuts(cellCutType_);

    if (debug)
    {
        Pout<< "isoSurfaceCell : candidate cells cut "
            << nCutCells_
            << " blocked " << nBlockedCells
            << " total " << mesh_.nCells() << endl;
    }

    if (debug && isA<fvMesh>(mesh))
    {
        const auto& fvmesh = refCast<const fvMesh>(mesh);

        volScalarField debugField
        (
            IOobject
            (
                "isoSurfaceCell.cutType",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvmesh,
            dimensionedScalar(dimless, Zero)
        );

        auto& debugFld = debugField.primitiveFieldRef();

        forAll(cellCutType_, celli)
        {
            debugFld[celli] = cellCutType_[celli];
        }

        Pout<< "Writing cut types:"
            << debugField.objectPath() << endl;

        debugField.write();
    }


    DynamicList<point> snappedPoints(nCutCells_);

    // Per cc -1 or a point inside snappedPoints.
    labelList snappedCc;
    if (regularise)
    {
        calcSnappedCc
        (
            isTet,
            cellValues,
            pointValues,
            snappedPoints,
            snappedCc
        );
    }
    else
    {
        snappedCc.setSize(mesh_.nCells());
        snappedCc = -1;
    }

    if (debug)
    {
        Pout<< "isoSurfaceCell : shifted " << snappedPoints.size()
            << " cell centres to intersection." << endl;
    }

    snappedPoints.shrink();
    label nCellSnaps = snappedPoints.size();

    // Per point -1 or a point inside snappedPoints.
    labelList snappedPoint;
    if (regularise)
    {
        calcSnappedPoint
        (
            isTet,
            cellValues,
            pointValues,
            snappedPoints,
            snappedPoint
        );
    }
    else
    {
        snappedPoint.setSize(mesh_.nPoints());
        snappedPoint = -1;
    }

    if (debug)
    {
        Pout<< "isoSurfaceCell : shifted " << snappedPoints.size()-nCellSnaps
            << " vertices to intersection." << endl;
    }


    // Use a triSurface as a temporary for various operations
    triSurface tmpsurf;

    {
        DynamicList<point> triPoints(nCutCells_);
        DynamicList<label> triMeshCells(nCutCells_);

        generateTriPoints
        (
            cellValues,
            pointValues,

            mesh_.cellCentres(),
            mesh_.points(),

            snappedPoints,
            snappedCc,
            snappedPoint,

            triPoints,
            triMeshCells
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : generated " << triMeshCells.size()
                << " unmerged triangles." << endl;
        }


        label nOldPoints = triPoints.size();

        // Trimmed to original triangle
        DynamicList<label> trimTriMap;
        // Trimmed to original point
        labelList trimTriPointMap;
        if (getClipBounds().good())
        {
            isoSurfacePoint::trimToBox
            (
                treeBoundBox(getClipBounds()),
                triPoints,              // new points
                trimTriMap,             // map from (new) triangle to original
                trimTriPointMap,        // map from (new) point to original
                interpolatedPoints_,    // labels of newly introduced points
                interpolatedOldPoints_, // and their interpolation
                interpolationWeights_
            );
            triMeshCells = labelField(triMeshCells, trimTriMap);
        }


        // Merge points and compact out non-valid triangles
        labelList triMap;           // merged to unmerged triangle
        tmpsurf = stitchTriPoints
        (
            regularise,         // check for duplicate tris
            triPoints,
            triPointMergeMap_,  // unmerged to merged point
            triMap              // merged to unmerged triangle
        );

        if (debug)
        {
            Pout<< "isoSurfaceCell : generated " << triMap.size()
                << " merged triangles." << endl;
        }

        if (getClipBounds().good())
        {
            // Adjust interpolatedPoints_
            inplaceRenumber(triPointMergeMap_, interpolatedPoints_);

            // Adjust triPointMergeMap_
            labelList newTriPointMergeMap(nOldPoints, -1);
            forAll(trimTriPointMap, trimPointI)
            {
                label oldPointI = trimTriPointMap[trimPointI];
                if (oldPointI >= 0)
                {
                    label pointI = triPointMergeMap_[trimPointI];
                    if (pointI >= 0)
                    {
                        newTriPointMergeMap[oldPointI] = pointI;
                    }
                }
            }
            triPointMergeMap_.transfer(newTriPointMergeMap);
        }

        meshCells_.setSize(triMap.size());
        forAll(triMap, i)
        {
            meshCells_[i] = triMeshCells[triMap[i]];
        }
    }


    if (debug)
    {
        Pout<< "isoSurfaceCell : checking " << tmpsurf.size()
            << " triangles for validity." << endl;

        forAll(tmpsurf, triI)
        {
            triSurfaceTools::validTri(tmpsurf, triI);
        }
    }


    if (regularise)
    {
        List<FixedList<label, 3>> faceEdges;
        labelList edgeFace0, edgeFace1;
        Map<labelList> edgeFacesRest;


        while (true)
        {
            // Calculate addressing
            calcAddressing
            (
                tmpsurf,
                faceEdges,
                edgeFace0,
                edgeFace1,
                edgeFacesRest
            );

            // See if any dangling triangles
            boolList keepTriangles;
            label nDangling = markDanglingTriangles
            (
                faceEdges,
                edgeFace0,
                edgeFace1,
                edgeFacesRest,
                keepTriangles
            );

            if (debug)
            {
                Pout<< "isoSurfaceCell : detected " << nDangling
                    << " dangling triangles." << endl;
            }

            if (nDangling == 0)
            {
                break;
            }

            // Create face map (new to old)
            labelList subsetTriMap(findIndices(keepTriangles, true));

            labelList subsetPointMap;
            labelList reversePointMap;
            tmpsurf = subsetMesh
            (
                tmpsurf,
                subsetTriMap,
                reversePointMap,
                subsetPointMap
            );
            meshCells_ = labelField(meshCells_, subsetTriMap);
            inplaceRenumber(reversePointMap, triPointMergeMap_);
        }
    }


    // Transfer to mesh storage. Note, an iso-surface has no zones
    {
        // Recover the pointField
        pointField pts;
        tmpsurf.swapPoints(pts);

        // Transcribe from triFace to face
        faceList faces;
        tmpsurf.triFaceFaces(faces);

        tmpsurf.clearOut();

        Mesh updated(std::move(pts), std::move(faces), surfZoneList());

        this->Mesh::transfer(updated);
    }
}


// ************************************************************************* //
