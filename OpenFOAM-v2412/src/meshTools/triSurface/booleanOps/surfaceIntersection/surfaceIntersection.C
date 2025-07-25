/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceIntersection.H"
#include "triSurfaceSearch.H"
#include "OBJstream.H"
#include "labelPairHashes.H"
#include "bitSet.H"
#include "triSurface.H"
#include "pointIndexHit.H"
#include "mergePoints.H"
#include "plane.H"
#include "edgeIntersections.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceIntersection, 0);
}

const Foam::Enum
<
    Foam::surfaceIntersection::intersectionType
>
Foam::surfaceIntersection::selfIntersectionNames
({
    { intersectionType::SELF, "self" },
    { intersectionType::SELF_REGION, "region" },
    { intersectionType::NONE, "none" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceIntersection::setOptions(const dictionary& dict)
{
    dict.readIfPresent("tolerance",       tolerance_);
    dict.readIfPresent("allowEdgeHits",   allowEdgeHits_);
    dict.readIfPresent("snap",            snapToEnd_);
    dict.readIfPresent("warnDegenerate",  warnDegenerate_);
}


void Foam::surfaceIntersection::storeIntersection
(
    const enum intersectionType cutFrom,
    const labelList& facesA,
    const label faceB,
    const UList<point>& allCutPoints,
    const label cutPointId,
    DynamicList<edge>& allCutEdges
)
{
    // Our lookup for two faces - populate with faceB (invariant)
    // Normally always have face from the first surface as first element
    labelPair twoFaces(faceB, faceB);

    forAll(facesA, facesAI)
    {
        const label faceA = facesA[facesAI];

        switch (cutFrom)
        {
            case surfaceIntersection::FIRST:
            {
                // faceA from 1st, faceB from 2nd
                twoFaces.first() = faceA;
                break;
            }
            case surfaceIntersection::SECOND:
            {
                // faceA from 2nd, faceB from 1st
                twoFaces.second() = faceA;
                break;
            }
            case surfaceIntersection::SELF:
            case surfaceIntersection::SELF_REGION:
            {
                // Lookup should be commutativity - use sorted order
                if (faceA < faceB)
                {
                    twoFaces.first()  = faceA;
                    twoFaces.second() = faceB;
                }
                else
                {
                    twoFaces.first()  = faceB;
                    twoFaces.second() = faceA;
                }
                break;
            }

            case surfaceIntersection::NONE:
                return;
                break;
        }


        // Get existing edge, or create a null edge (with -1)
        edge& thisEdge = facePairToEdge_(twoFaces);
        const label pointCount = thisEdge.count();

        if (pointCount == 0)
        {
            // First intersection of the faces - record it.
            thisEdge.insert(cutPointId);

            if (debug & 4)
            {
                Pout<< "intersect faces " << twoFaces
                    << " point-1: " << cutPointId << " = "
                    << allCutPoints[cutPointId] << endl;
            }
            continue;
        }
        else if (pointCount == 2)
        {
            // This occurs for ugly surfaces with shards that result in multiple
            // cuts very near a snapped end point.
            if (debug & 4)
            {
                Pout<< "suppressed double intersection " << twoFaces
                    << thisEdge << endl;
            }
            continue;
        }

        if (thisEdge.insert(cutPointId))
        {
            // Second intersection of the faces - this is an edge,
            // with special treatment:
            // - avoid duplicate points: addressed by the insert() above
            // - avoid degenerate lengths
            // - avoid duplicate edges - can occur with really dirty geometry

            if (edgeToId_.found(thisEdge))
            {
                // Already have this edgeId, but not for this intersection.
                thisEdge.sort();
                if (facePairToEdge_.insert(twoFaces, thisEdge))
                {
                    if (debug & 4)
                    {
                        Pout<< "reuse edge - faces " << twoFaces << " edge#"
                            << edgeToId_[thisEdge] << " edge " << thisEdge
                            << " = " << thisEdge.line(allCutPoints)
                            << endl;
                    }
                }
            }
            else if (thisEdge.mag(allCutPoints) < SMALL)
            {
                // Degenerate length
                // - eg, end snapping was disabled or somehow failed.

                // Don't normally emit warnings, since these also arise for
                // manifold connections. For example,
                //
                //   e1|  /e2
                //     | /
                //     |/
                // ----.---- plane
                //
                // The plane is correctly pierced at the '.' by both edge-1
                // and edge-2, which belong to the same originating face.

                // Filter/merge away the extraneous points later.
                if (warnDegenerate_ > 0)
                {
                    --warnDegenerate_;
                    WarningInFunction
                        << "Degenerate edge between faces " << twoFaces
                        << " on 1st/2nd surface with points "
                        << thisEdge.line(allCutPoints)
                        << endl;
                }
                else if (debug & 4)
                {
                    Pout<< "degenerate edge face-pair " << twoFaces << " "
                        << thisEdge[0] << " point " << allCutPoints[thisEdge[0]]
                        << endl;
                }

                // This is a failed edge - undo this second interaction
                thisEdge.erase(cutPointId);
            }
            else
            {
                // This is a new edge.
                const label edgeId = allCutEdges.size();

                if (facePairToEdgeId_.insert(twoFaces, edgeId))
                {
                    // Record complete (line) intersection of two faces
                    thisEdge.sort();
                    edgeToId_.insert(thisEdge, edgeId);
                    allCutEdges.append(thisEdge);

                    if (debug & 4)
                    {
                        Pout<< "create edge - faces " << twoFaces << " edge#"
                            << edgeId << " edge " << thisEdge
                            << " = " << thisEdge.line(allCutPoints)
                            << endl;
                    }
                }
                else
                {
                    // Faces already had an intersection
                    // This should not fail, but for safety.
                    Info<<"WARN " << twoFaces
                        << " already intersected= " << thisEdge << endl;
                    thisEdge.erase(cutPointId);
                }
            }
        }
        else
        {
            // Duplicate point - usually zero-length edge from snapping
            // - can discard this face/face interaction entirely
            facePairToEdge_.erase(twoFaces);
        }
    }
}


// Classify cut of edge of surface1 with surface2:
// 1- point of edge hits point on surface2
// 2- edge pierces point on surface2
// 3- point of edge hits edge on surface2
// 4- edge pierces edge on surface2
// 5- point of edge hits face on surface2
// 6- edge pierces face on surface2
//
// Note that handling of 2 and 4 should be the same but with surface1 and
// surface2 reversed.
void Foam::surfaceIntersection::classifyHit
(
    const triSurface& surf1,
    const scalarField& surf1PointTol,
    const triSurface& surf2,
    const enum intersectionType cutFrom,
    const label edgeI,
    const pointIndexHit& pHit,

    DynamicList<point>& allCutPoints,
    DynamicList<edge>& allCutEdges,
    List<DynamicList<label>>& surfEdgeCuts
)
{
    const edge& e1 = surf1.edges()[edgeI];

    const labelList& facesA = surf1.edgeFaces()[edgeI];

    // Label of face on surface2 edgeI intersected
    label surf2Facei = pHit.index();

    // Classify point on surface2

    const triSurface::face_type& f2 = surf2.localFaces()[surf2Facei];
    const pointField& surf1Pts = surf1.localPoints();
    const pointField& surf2Pts = surf2.localPoints();

    label nearType, nearLabel;

    f2.nearestPointClassify(pHit.hitPoint(), surf2Pts, nearType, nearLabel);

    // Classify points on edge of surface1
    const label edgeEnd =
        classify
        (
            surf1PointTol[e1.start()],
            surf1PointTol[e1.end()],
            pHit.hitPoint(),
            e1,
            surf1Pts
        );

    if (nearType == triPointRef::POINT)
    {
        if (edgeEnd >= 0)
        {
            // 1. Point hits point. Do nothing.
            if (debug & 2)
            {
                Pout<< "hit-type[1] " << pHit.point() << " is surf1:"
                    << " end point of edge[" << edgeI << "] " << e1
                    << "==" << e1.line(surf1Pts)
                    << " surf2: vertex " << f2[nearLabel]
                    << " coord:" << surf2Pts[f2[nearLabel]]
                    << " - suppressed" << endl;
            }
        }
        else
        {
            // 2. Edge hits point. Cut edge with new point.
            label cutPointId = -1;
            const label nearVert = f2[nearLabel];

            // For self-intersection, we have tolerances for each point
            // (surf2 is actually surf1) so we shift the hit to coincide
            // identically.
            if
            (
                cutFrom == surfaceIntersection::SELF
             || cutFrom == surfaceIntersection::SELF_REGION
            )
            {
                const point& nearPt = surf1Pts[nearVert];

                if
                (
                    pHit.hitPoint().dist(nearPt)
                  < surf1PointTol[nearVert]
                )
                {
                    cutPointId = allCutPoints.size();

                    if (snapToEnd_)
                    {
                        if (snappedEnds_.insert(nearVert, cutPointId))
                        {
                            // Initial snap
                            allCutPoints.append(nearPt);
                        }
                        else
                        {
                            // Already snapped this point.
                            cutPointId = snappedEnds_[nearVert];
                        }
                    }
                    else
                    {
                        allCutPoints.append(nearPt);
                    }
                }
            }

            if (debug & 2)
            {
                Pout<< "hit-type[2] " << pHit.point() << " is surf1:"
                    << " from edge[" << edgeI << "] " << e1
                    << " surf2: vertex " << f2[nearLabel]
                    << " coord:" << surf2Pts[f2[nearLabel]]
                    << " - "
                    << (cutPointId >= 0 ? "snapped" : "stored") << endl;
            }

            if (cutPointId == -1)
            {
                cutPointId = allCutPoints.size();
                allCutPoints.append(pHit.hitPoint());
            }
            surfEdgeCuts[edgeI].append(cutPointId);

            const labelList& facesB = surf2.pointFaces()[f2[nearLabel]];
            forAll(facesB, faceBI)
            {
                storeIntersection
                (
                    cutFrom,
                    facesA,
                    facesB[faceBI],
                    allCutPoints,
                    cutPointId,
                    allCutEdges
                );
            }
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        if (edgeEnd >= 0)
        {
            // 3. Point hits edge.
            // Normally do nothing on this side since the reverse
            // (edge hits point) is handled by 2.
            // However, if the surfaces are separated by a minor gap,
            // the end-point of a tolerance-extended edge can intersect another
            // edge without itself being intersected by an edge.

            const label edge2I = getEdge(surf2, surf2Facei, nearLabel);
            const edge& e2 = surf2.edges()[edge2I];
            const label nearVert = e1[edgeEnd];

            label cutPointId = -1;

            // Storage treatment
            // =0: nothing/ignore
            // >0: store point/edge-cut. Attempt to create new edge.
            // <0: store point/edge-cut only
            int handling = (allowEdgeHits_ ? 1 : 0);
            if
            (
                allowEdgeHits_
             &&
                (
                    cutFrom == surfaceIntersection::SELF
                 || cutFrom == surfaceIntersection::SELF_REGION
                )
            )
            {
                // The edge-edge intersection is hashed as an 'edge' to
                // exploit the commutative lookup.
                // Ie, only do the cut once
                const edge intersect(edgeI, edge2I);

                if (e2.found(nearVert))
                {
                    // Actually the same as #1 above, but missed due to
                    // tolerancing
                    handling = 0; // suppress
                }
                else if (edgeEdgeIntersection_.insert(intersect))
                {
                    const point& nearPt = surf1Pts[nearVert];

                    if
                    (
                        pHit.hitPoint().dist(nearPt)
                      < surf1PointTol[nearVert]
                    )
                    {
                        cutPointId = allCutPoints.size();

                        if (snapToEnd_)
                        {
                            if (snappedEnds_.insert(nearVert, cutPointId))
                            {
                                // Initial snap
                                allCutPoints.append(nearPt);
                            }
                            else
                            {
                                // Already snapped this point.
                                cutPointId = snappedEnds_[nearVert];
                                handling = 2;  // cached
                            }
                        }
                        else
                        {
                            allCutPoints.append(nearPt);
                        }
                    }
                }
                else
                {
                    handling = 0; // ignore - already did this interaction
                }
            }

            if (debug & 2)
            {
                Pout<< "hit-type[3] " << pHit.point() << " is surf1:"
                    << " end point of edge[" << edgeI << "] " << e1
                    << "==" << e1.line(surf1Pts)
                    << " surf2: edge[" << edge2I << "] " << e2
                    << " coords:" << e2.line(surf2Pts)
                    << " - "
                    << (
                           handling > 1
                         ? "cached" : handling
                         ? "stored" : "suppressed"
                       ) << endl;
            }

            if (handling)
            {
                if (cutPointId == -1)
                {
                    cutPointId = allCutPoints.size();
                    allCutPoints.append(pHit.hitPoint());
                }
                surfEdgeCuts[edgeI].append(cutPointId);
            }

            if (handling > 0)
            {
                const labelList& facesB = surf2.edgeFaces()[edge2I];
                forAll(facesB, faceBI)
                {
                    storeIntersection
                    (
                        cutFrom,
                        facesA,
                        facesB[faceBI],
                        allCutPoints,
                        cutPointId,
                        allCutEdges
                    );
                }
            }
        }
        else
        {
            // 4. Edge hits edge.

            // Cut edge with new point (creates duplicates when
            // doing the surf2 with surf1 intersection but these
            // are merged later on)

            // edge hits all faces on surf2 connected to the edge
            //
            // The edge-edge intersection is symmetric, store only once.
            // - When intersecting two surfaces, note which edges are cut each
            //   time, but only create an edge from the first pass.
            // - For self-intersection, it is slightly trickier if we don't
            //   want too many duplicate points.

            const label edge2I = getEdge(surf2, surf2Facei, nearLabel);
            const edge& e2 = surf2.edges()[edge2I];
            label cutPointId = -1;

            // Storage treatment
            // =0: nothing/ignore
            // >0: store point/edge-cut. Attempt to create new edge.
            // <0: store point/edge-cut only
            int handling = 0;
            switch (cutFrom)
            {
                case surfaceIntersection::FIRST:
                {
                    handling = 1;
                    break;
                }
                case surfaceIntersection::SECOND:
                {
                    handling = -1;
                    break;
                }
                case surfaceIntersection::SELF:
                case surfaceIntersection::SELF_REGION:
                {
                    // The edge-edge intersection is hashed as an 'edge' to
                    // exploit the commutative lookup.
                    // Ie, only do the cut once
                    const edge intersect(edgeI, edge2I);

                    if (edgeEdgeIntersection_.insert(intersect))
                    {
                        handling = 1;
                        forAll(e1, edgepti)
                        {
                            const label endId = e1[edgepti];
                            const point& nearPt = surf1Pts[endId];

                            if
                            (
                                pHit.hitPoint().dist(nearPt)
                              < surf1PointTol[endId]
                            )
                            {
                                cutPointId = allCutPoints.size();

                                if (snapToEnd_)
                                {
                                    if (snappedEnds_.insert(endId, cutPointId))
                                    {
                                        // First time with this end-point
                                        allCutPoints.append(nearPt);
                                    }
                                    else
                                    {
                                        // Already seen this end point
                                        cutPointId = snappedEnds_[endId];
                                        handling = 2;  // cached
                                    }
                                }
                                else
                                {
                                    allCutPoints.append(nearPt);
                                }

                                break;
                            }
                        }
                    }

                    break;
                }

                case surfaceIntersection::NONE:
                    return;
                    break;
            }

            if (debug & 2)
            {
                Pout<< "hit-type[4] " << pHit.point() << " is surf1:"
                    << " from edge[" << edgeI << "] " << e1
                    << "==" << e1.line(surf1Pts)
                    << " surf2: edge[" << edge2I << "] " << e2
                    << " coords:" << e2.line(surf2Pts)
                    << " - "
                    << (
                         handling < 0
                       ? "cut-point" : handling
                       ? "stored" : "suppressed"
                       )
                    << endl;
            }

            if (handling)
            {
                if (cutPointId == -1)
                {
                    cutPointId = allCutPoints.size();
                    allCutPoints.append(pHit.hitPoint());
                }
                surfEdgeCuts[edgeI].append(cutPointId);
            }

            if (handling)
            {
                const vector eVec = e1.unitVec(surf1Pts);

                const labelList& facesB = surf2.edgeFaces()[edge2I];
                forAll(facesB, faceBI)
                {
                    // Intersecting edge should be non-coplanar with face
                    if
                    (
                        mag((surf2.faceNormals()[facesB[faceBI]] & eVec))
                      > 0.01
                    )
                    {
                        storeIntersection
                        (
                            cutFrom,
                            facesA,
                            facesB[faceBI],
                            allCutPoints,
                            cutPointId,
                            allCutEdges
                        );
                    }
                }
            }
        }
    }
    else
    {
        if (edgeEnd >= 0)
        {
            // 5. Point hits face. Do what? Introduce
            // point & triangulation in face?

            // Look exactly at what side (of surf2) edge is. Leave out ones on
            // inside of surf2 (i.e. on opposite side of normal)

            // Vertex on/near surf2; vertex away from surf2
            // otherVert on outside of surf2
            const label nearVert  = (edgeEnd == 0 ? e1.start() : e1.end());
            const label otherVert = (edgeEnd == 0 ? e1.end() : e1.start());

            const point& nearPt  = surf1Pts[nearVert];
            const point& otherPt = surf1Pts[otherVert];

            const vector eVec = otherPt - nearPt;

            if ((surf2.faceNormals()[surf2Facei] & eVec) > 0)
            {
                // map to nearVert
                // Reclassify as normal edge-face pierce (see below)
                bool cached = false;

                label cutPointId = allCutPoints.size();
                if (snapToEnd_)
                {
                    if (snappedEnds_.insert(nearVert, cutPointId))
                    {
                        // First time with this end-point
                        allCutPoints.append(nearPt);
                    }
                    else
                    {
                        // Already seen this end point
                        cutPointId = snappedEnds_[nearVert];
                        cached = true;
                    }
                }
                else
                {
                    allCutPoints.append(nearPt);
                }

                surfEdgeCuts[edgeI].append(cutPointId);

                if (debug & 2)
                {
                    Pout<< "hit-type[5] " << pHit.point()
                        << " shifted to " << nearPt
                        << " from edge[" << edgeI << "] " << e1
                        << "==" << e1.line(surf1Pts)
                        << " hits surf2 face[" << surf2Facei << "]"
                        << " - "
                        << (cached ? "cached" : "stored") << endl;
                }

                // edge hits single face only
                storeIntersection
                (
                    cutFrom,
                    facesA,
                    surf2Facei,
                    allCutPoints,
                    cutPointId,
                    allCutEdges
                );
            }
            else
            {
                if (debug & 2)
                {
                    Pout<< "hit-type[5] " << pHit.point()
                        << " from edge[" << edgeI << "] " << e1
                        << " hits inside of surf2 face[" << surf2Facei << "]"
                        << " - discarded" << endl;
                }
            }
        }
        else
        {
            // 6. Edge pierces face. 'Normal' situation.
            if (debug & 2)
            {
                Pout<< "hit-type[6] " << pHit.point()
                    << " from edge[" << edgeI << "] " << e1
                    << "==" << e1.line(surf1Pts)
                    << " hits surf2 face[" << surf2Facei << "]"
                    << " - stored" << endl;
            }

            const label cutPointId = allCutPoints.size();
            allCutPoints.append(pHit.hitPoint());
            surfEdgeCuts[edgeI].append(cutPointId);

            // edge hits single face only
            storeIntersection
            (
                cutFrom,
                facesA,
                surf2Facei,
                allCutPoints,
                cutPointId,
                allCutEdges
            );
        }
    }
}


// Cut all edges of surf1 with surf2. Sets
// - cutPoints          : coordinates of cutPoints
// - cutEdges           : newly created edges between cutPoints
// - facePairToVertex   : hash from face1I and face2I to edge
// - facePairToEdgeId   : hash from face1I and face2I to index in cutEdge
// - surfEdgeCuts       : gives for each edge the cutPoints
//                        (in order from start to end)
//
void Foam::surfaceIntersection::doCutEdges
(
    const triSurface& surf1,
    const triSurfaceSearch& querySurf2,
    const enum intersectionType cutFrom,

    DynamicList<point>& allCutPoints,
    DynamicList<edge>& allCutEdges,
    List<DynamicList<label>>& surfEdgeCuts
)
{
    const scalar oldTol = intersection::setPlanarTol(tolerance_);

    const pointField& surf1Pts = surf1.localPoints();

    // Calculate local (to point) tolerance based on min edge length.
    scalarField surf1PointTol(surf1Pts.size());

    forAll(surf1PointTol, pointi)
    {
        surf1PointTol[pointi] = tolerance_ * minEdgeLen(surf1, pointi);
    }

    const indexedOctree<treeDataPrimitivePatch<triSurface>>& searchTree
        = querySurf2.tree();

    if
    (
        cutFrom == surfaceIntersection::SELF
     || cutFrom == surfaceIntersection::SELF_REGION
    )
    {
        // An edge may intersect multiple faces
        // - mask out faces that have already been hit before trying again
        // - never intersect with faces attached to the edge itself
        DynamicList<label> maskFaces(32);

        // Optionally prevent intersection within a single region.
        // Like self-intersect, but only if regions are different
        bitSet maskRegions(32);

        treeDataTriSurface::findAllIntersectOp
            allIntersectOp(searchTree, maskFaces);

        forAll(surf1.edges(), edgeI)
        {
            const edge& e = surf1.edges()[edgeI];
            const vector edgeVec = e.vec(surf1Pts);

            // Extend start/end by 1/2 tolerance - ensures cleaner cutting
            const point ptStart =
                surf1Pts[e.start()] - 0.5*surf1PointTol[e.start()]*edgeVec;
            const point ptEnd =
                surf1Pts[e.end()]   + 0.5*surf1PointTol[e.end()]*edgeVec;

            maskRegions.clear();
            if (cutFrom == surfaceIntersection::SELF_REGION)
            {
                for (auto& facei : surf1.edgeFaces()[edgeI])
                {
                    maskRegions.set(surf1[facei].region());
                }
            }

            // Never intersect with faces attached directly to the edge itself,
            // nor with faces attached to its end points. This mask contains
            // some duplicates, but filtering them out is less efficient.
            maskFaces = surf1.pointFaces()[e.start()];
            maskFaces.append(surf1.pointFaces()[e.end()]);

            while (true)
            {
                pointIndexHit pHit = searchTree.findLine
                (
                    ptStart,
                    ptEnd,
                    allIntersectOp
                );

                if (!pHit.hit())
                {
                    break;
                }

                maskFaces.append(pHit.index());

                if (maskRegions.test(surf1[pHit.index()].region()))
                {
                    continue;
                }

                classifyHit
                (
                    surf1,
                    surf1PointTol,
                    surf1,
                    cutFrom,
                    edgeI,
                    pHit,

                    allCutPoints,
                    allCutEdges,
                    surfEdgeCuts
                );
            }
        }
    }
    else
    {
        const triSurface& surf2 = querySurf2.surface();

        forAll(surf1.edges(), edgeI)
        {
            const edge& e = surf1.edges()[edgeI];
            const vector edgeVec = e.vec(surf1Pts);

            const point tolVec = intersection::planarTol()*(edgeVec);
            const scalar tolDim = mag(tolVec);

            // Extend start/end by 1/2 tolerance - ensures cleaner cutting
            point ptStart =
                surf1Pts[e.start()] - 0.5*surf1PointTol[e.start()]*edgeVec;
            const point ptEnd =
                surf1Pts[e.end()]   + 0.5*surf1PointTol[e.end()]*edgeVec;

            bool doTrack = false;
            do
            {
                pointIndexHit pHit = searchTree.findLine(ptStart, ptEnd);

                if (!pHit.hit())
                {
                    break;
                }

                classifyHit
                (
                    surf1,
                    surf1PointTol,
                    surf2,
                    cutFrom,
                    edgeI,
                    pHit,

                    allCutPoints,
                    allCutEdges,
                    surfEdgeCuts
                );

                if (tolerance_ > 0)
                {
                    if (pHit.hitPoint().dist(ptEnd) < tolDim)
                    {
                        // Near the end => done
                        doTrack = false;
                    }
                    else
                    {
                        // Continue tracking a bit further on
                        ptStart = pHit.hitPoint() + tolVec;
                        doTrack = true;
                    }
                }
            }
            while (doTrack);  // execute at least once
        }
    }
    if (debug & 2)
    {
        Pout<< endl;
    }

    // These temporaries are now unneeded:
    edgeEdgeIntersection_.clear();
    snappedEnds_.clear();

    intersection::setPlanarTol(oldTol);
}


void Foam::surfaceIntersection::joinDisconnected
(
    DynamicList<edge>& allCutEdges
)
{
    // This simple heuristic seems to work just as well (or better) than
    // more complicated schemes
    //
    // For any face/face intersection that only appears once,
    // consider which other faces/points are involved and connect between
    // those points.
    // Just do a simple connect-the-dots?

    Pair<Map<labelPairHashSet>> missedFacePoint;

    // Stage 1:
    // - Extract "faceId -> (faceId, pointId)"
    //   for all face/face pairs that only have one interaction
    forAllConstIters(facePairToEdge_, iter)
    {
        const labelPair& twoFaces = iter.key();
        const edge& e = iter.val();

        if (e.count() == 1)
        {
            // minVertex = -1 (unused), maxVertex = pointId
            const label pointId = e.maxVertex();

            missedFacePoint[0](twoFaces[0]).insert
            (
                labelPair(twoFaces[1], pointId)
            );

            missedFacePoint[1](twoFaces[1]).insert
            (
                labelPair(twoFaces[0], pointId)
            );
        }
    }


    // Stage 2:
    // - anything with two cross-interactions could cause a new edge:

    edgeHashSet newEdges;
    forAll(missedFacePoint, sidei)
    {
        const auto& mapping = missedFacePoint[sidei];

        forAllConstIters(mapping, iter)
        {
            const auto& connect = iter.val();

            if (connect.size() == 2)
            {
                // exactly two face/face cross-interactions

                edge e;
                for (const auto& facePoint : connect)
                {
                    e.insert(facePoint.second());
                }
                e.sort();

                // Only consider edges with two unique ends,
                // and do not introduce duplicates
                if (e.count() == 2 && !edgeToId_.found(e))
                {
                    newEdges.insert(e);
                }
            }
        }
    }

    label edgeId = allCutEdges.size();
    edgeList newEdgesLst = newEdges.sortedToc();
    for (const auto& e : newEdgesLst)
    {
        // Record complete (line) intersection of two faces
        allCutEdges.append(e);
        edgeToId_.insert(e, edgeId);
        ++edgeId;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceIntersection::surfaceIntersection()
:
    tolerance_(1e-3),
    allowEdgeHits_(true),
    snapToEnd_(true),
    warnDegenerate_(0),
    cutPoints_(0),
    cutEdges_(0),
    facePairToEdge_(0),
    facePairToEdgeId_(0),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurfaceSearch& query1,
    const triSurfaceSearch& query2,
    const dictionary& dict
)
:
    tolerance_(1e-3),
    allowEdgeHits_(true),
    snapToEnd_(true),
    warnDegenerate_(0),
    cutPoints_(0),
    cutEdges_(0),
    facePairToEdge_(2*max(query1.surface().size(), query2.surface().size())),
    facePairToEdgeId_(2*max(query1.surface().size(), query2.surface().size())),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{
    setOptions(dict);

    const triSurface& surf1 = query1.surface();
    const triSurface& surf2 = query2.surface();

    //
    // Cut all edges of surf1 with surf2.
    //
    if (debug)
    {
        Pout<< "Cutting surf1 edges" << endl;
    }


    DynamicList<edge>  allCutEdges(surf1.nEdges()/20);
    DynamicList<point> allCutPoints(surf1.nPoints()/20);


    // From edge to cut index on surface1
    List<DynamicList<label>> edgeCuts1(query1.surface().nEdges());

    // 1st surf (for labelPair order)
    doCutEdges
    (
        surf1,
        query2,
        surfaceIntersection::FIRST,
        allCutPoints,
        allCutEdges,
        edgeCuts1
    );
    // Transfer to straight labelListList
    transfer(edgeCuts1, surf1EdgeCuts_);


    //
    // Cut all edges of surf2 with surf1.
    //
    if (debug)
    {
        Pout<< "Cutting surf2 edges" << endl;
    }

    // From edge to cut index
    List<DynamicList<label>> edgeCuts2(query2.surface().nEdges());

    // 2nd surf (for labelPair order)
    doCutEdges
    (
        surf2,
        query1,
        surfaceIntersection::SECOND,
        allCutPoints,
        allCutEdges,
        edgeCuts2
    );

    // join disconnected intersection points
    joinDisconnected(allCutEdges);

    // Transfer to straight label(List)List
    transfer(edgeCuts2, surf2EdgeCuts_);
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);

    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;

        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OBJstream("intEdges.obj").write(cutEdges_, cutPoints_);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);

        Pout<< "Dumping cut edges of surface2 to surf2EdgeCuts.obj" << endl;
        OFstream edge2Stream("surf2EdgeCuts.obj");
        writeIntersectedEdges(surf2, surf2EdgeCuts_, edge2Stream);
    }

    // Temporaries
    facePairToEdge_.clear();

    // Cleanup any duplicate cuts?
    // mergeEdges();
}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurfaceSearch& query1,
    const dictionary& dict
)
:
    tolerance_(1e-3),
    allowEdgeHits_(true),
    snapToEnd_(true),
    warnDegenerate_(0),
    cutPoints_(0),
    cutEdges_(0),
    facePairToEdge_(2*query1.surface().size()),
    facePairToEdgeId_(2*query1.surface().size()),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{
    setOptions(dict);

    const intersectionType cutFrom = selfIntersectionNames.getOrDefault
    (
        "intersectionMethod",
        dict,
        intersectionType::SELF
    );

    if (cutFrom == intersectionType::NONE)
    {
        if (debug)
        {
            Pout<< "Skipping self-intersection (selected: none)" << endl;
        }

        // Temporaries
        facePairToEdge_.clear();
        facePairToEdgeId_.clear();

        return;
    }

    const triSurface& surf1 = query1.surface();

    //
    // Cut all edges of surf1 with surf1 itself.
    //
    if (debug)
    {
        Pout<< "Cutting surf1 edges" << endl;
    }

    DynamicList<edge>  allCutEdges;
    DynamicList<point> allCutPoints;

    // From edge to cut index on surface1
    List<DynamicList<label>> edgeCuts1(query1.surface().nEdges());

    // self-intersection
    doCutEdges
    (
        surf1,
        query1,
        cutFrom,
        allCutPoints,
        allCutEdges,
        edgeCuts1
    );

    // join disconnected intersection points
    joinDisconnected(allCutEdges);

    // Transfer to straight label(List)List
    transfer(edgeCuts1, surf1EdgeCuts_);
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);

    // Short-circuit
    if (cutPoints_.empty() && cutEdges_.empty())
    {
        if (debug)
        {
            Pout<< "Empty intersection" << endl;
        }
        return;
    }

    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated and compressed:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;


        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OBJstream("intEdges.obj").write(cutEdges_, cutPoints_);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);
    }

    // Temporaries
    facePairToEdge_.clear();

    // // Cleanup any duplicate cuts?
    // mergeEdges();
}


Foam::surfaceIntersection::surfaceIntersection
(
    const triSurface& surf1,
    const edgeIntersections& intersections1,
    const triSurface& surf2,
    const edgeIntersections& intersections2
)
:
    tolerance_(1e-3),
    allowEdgeHits_(true),
    snapToEnd_(true),
    warnDegenerate_(0),
    cutPoints_(0),
    cutEdges_(0),
    facePairToEdge_(2*max(surf1.size(), surf2.size())),
    facePairToEdgeId_(2*max(surf1.size(), surf2.size())),
    surf1EdgeCuts_(0),
    surf2EdgeCuts_(0)
{

    // All intersection Pout (so for both surfaces)
    DynamicList<edge>  allCutEdges((surf1.nEdges() + surf2.nEdges())/20);
    DynamicList<point> allCutPoints((surf1.nPoints() + surf2.nPoints())/20);


    // Cut all edges of surf1 with surf2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "Storing surf1 intersections" << endl;
    }

    {
        // From edge to cut index on surface1
        List<DynamicList<label>> edgeCuts1(surf1.nEdges());

        forAll(intersections1, edgeI)
        {
            const List<pointIndexHit>& intersections = intersections1[edgeI];

            forAll(intersections, i)
            {
                // edgeI intersects surf2. Store point.
                const pointIndexHit& pHit = intersections[i];
                const label cutPointId = allCutPoints.size();

                allCutPoints.append(pHit.hitPoint());
                edgeCuts1[edgeI].append(cutPointId);

                storeIntersection
                (
                    surfaceIntersection::FIRST,
                    surf1.edgeFaces()[edgeI],
                    pHit.index(),
                    allCutPoints,
                    cutPointId,
                    allCutEdges
                );
            }
        }

        // Transfer to straight labelListList
        transfer(edgeCuts1, surf1EdgeCuts_);
    }


    // Cut all edges of surf2 with surf1
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "Storing surf2 intersections" << endl;
    }

    {
        // From edge to cut index on surface2
        List<DynamicList<label>> edgeCuts2(surf2.nEdges());

        forAll(intersections2, edgeI)
        {
            const List<pointIndexHit>& intersections = intersections2[edgeI];

            forAll(intersections, i)
            {
                // edgeI intersects surf1. Store point.
                const pointIndexHit& pHit = intersections[i];
                const label cutPointId = allCutPoints.size();

                allCutPoints.append(pHit.hitPoint());
                edgeCuts2[edgeI].append(cutPointId);

                storeIntersection
                (
                    surfaceIntersection::SECOND,
                    surf2.edgeFaces()[edgeI],
                    pHit.index(),
                    allCutPoints,
                    cutPointId,
                    allCutEdges
                );
            }
        }

        // Transfer to surf2EdgeCuts_ (straight labelListList)
        transfer(edgeCuts2, surf2EdgeCuts_);
    }


    // Transfer to straight label(List)List
    cutEdges_.transfer(allCutEdges);
    cutPoints_.transfer(allCutPoints);


    if (debug)
    {
        Pout<< "surfaceIntersection : Intersection generated:"
            << endl
            << "    points:" << cutPoints_.size() << endl
            << "    edges :" << cutEdges_.size() << endl;

        Pout<< "surfaceIntersection : Writing intersection to intEdges.obj"
            << endl;

        OBJstream("intEdges.obj").write(cutEdges_, cutPoints_);

        // Dump all cut edges to files
        Pout<< "Dumping cut edges of surface1 to surf1EdgeCuts.obj" << endl;
        OFstream edge1Stream("surf1EdgeCuts.obj");
        writeIntersectedEdges(surf1, surf1EdgeCuts_, edge1Stream);

        Pout<< "Dumping cut edges of surface2 to surf2EdgeCuts.obj" << endl;
        OFstream edge2Stream("surf2EdgeCuts.obj");
        writeIntersectedEdges(surf2, surf2EdgeCuts_, edge2Stream);
    }

    // Temporaries
    facePairToEdge_.clear();

    // // Cleanup any duplicate cuts?
    // mergeEdges();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::surfaceIntersection::cutPoints() const
{
    return cutPoints_;
}


const Foam::edgeList& Foam::surfaceIntersection::cutEdges() const
{
    return cutEdges_;
}


const Foam::labelPairLookup& Foam::surfaceIntersection::facePairToEdgeId() const
{
    return facePairToEdgeId_;
}


const Foam::labelListList& Foam::surfaceIntersection::edgeCuts
(
    const bool isFirstSurf
) const
{
    if (isFirstSurf)
    {
        return surf1EdgeCuts_;
    }
    else
    {
        return surf2EdgeCuts_;
    }
}


const Foam::labelListList& Foam::surfaceIntersection::surf1EdgeCuts() const
{
    return surf1EdgeCuts_;
}


const Foam::labelListList& Foam::surfaceIntersection::surf2EdgeCuts() const
{
    return surf2EdgeCuts_;
}


void Foam::surfaceIntersection::mergePoints(const scalar mergeDist)
{
    labelList pointMap;

    label nChanged = Foam::inplaceMergePoints
    (
        cutPoints_,
        mergeDist,
        false,
        pointMap
    );

    if (nChanged)
    {
        forAll(cutEdges_, edgei)
        {
            edge& e = cutEdges_[edgei];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }

        forAll(surf1EdgeCuts_, edgei)
        {
            inplaceRenumber(pointMap, surf1EdgeCuts_[edgei]);
            inplaceUniqueSort(surf1EdgeCuts_[edgei]);
        }
        forAll(surf2EdgeCuts_, edgei)
        {
            inplaceRenumber(pointMap, surf2EdgeCuts_[edgei]);
            inplaceUniqueSort(surf2EdgeCuts_[edgei]);
        }
    }

    this->mergeEdges();
}


void Foam::surfaceIntersection::mergeEdges()
{
    edgeHashSet uniqEdges(2*cutEdges_.size());

    label nUniqEdges = 0;
    labelList edgeNumbering(cutEdges_.size(), -1);

    forAll(cutEdges_, edgeI)
    {
        const edge& e = cutEdges_[edgeI];

        // Remove degenerate and repeated edges
        // - reordering (e[0] < e[1]) is not really necessary
        if (e[0] != e[1] && uniqEdges.insert(e))
        {
            edgeNumbering[edgeI] = nUniqEdges;
            if (nUniqEdges != edgeI)
            {
                cutEdges_[nUniqEdges] = e;
            }
            cutEdges_[nUniqEdges].sort();
            ++nUniqEdges;
        }
    }

    // if (nUniqEdges < cutEdges_.size())
    // {
    //     // Additional safety, in case the edge was replaced?
    //     forAllIters(facePairToEdge_, iter)
    //     {
    //         iter.val() = edgeNumbering[iter.val()];
    //     }
    // }

    cutEdges_.setSize(nUniqEdges);  // truncate
}


// ************************************************************************* //
