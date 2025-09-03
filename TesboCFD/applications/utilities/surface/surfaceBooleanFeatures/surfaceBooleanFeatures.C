/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Application
    surfaceBooleanFeatures

Group
    grpSurfaceUtilities

Description
    Generates the extendedFeatureEdgeMesh for the interface between a boolean
    operation on two surfaces.

    Assumes that the orientation of the surfaces is correct:
    - if the operation is union or intersection, that both surface's normals
      (n) have the same orientation with respect to a point, i.e. surfaces
      A and B are orientated the same with respect to point x:

    \verbatim
       _______
      |       |--> n
      |    ___|___             x
      |A  |   |   |--> n
      |___|___|  B|
          |       |
          |_______|

    \endverbatim

    - if the operation is a subtraction, the surfaces should be oppositely
    oriented with respect to a point, i.e. for (A - B), then B's orientation
    should be such that x is "inside", and A's orientation such that x is
    "outside"

    \verbatim
       _______
      |       |--> n
      |    ___|___             x
      |A  |   |   |
      |___|___|  B|
          |  n <--|
          |_______|

    \endverbatim

    When the operation is performed - for union, all of the edges generated
    where one surfaces cuts another are all "internal" for union,
    and "external" for intersection, (B - A) and (A - B).
    This has been assumed, formal (dis)proof is invited.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "Time.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "triSurfaceSearch.H"
#include "triSurfaceMesh.H"
#include "OFstream.H"
#include "booleanSurface.H"
#include "edgeIntersections.H"
#include "meshTools.H"
#include "DynamicField.H"
#include "Enum.H"

#ifndef NO_CGAL

// Silence boost bind deprecation warnings (before CGAL-5.2.1)
#include "CGAL/version.h"
#if defined(CGAL_VERSION_NR) && (CGAL_VERSION_NR < 1050211000)
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#endif
#pragma clang diagnostic ignored "-Wbitwise-instead-of-logical"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include "CGALIndexedPolyhedron.H"
#include "PolyhedronReader.H"
typedef CGAL::AABB_face_graph_triangle_primitive
<
    Polyhedron, CGAL::Default, CGAL::Tag_false
> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

// Used boost::optional prior to CGAL-6.0
#if (CGAL_VERSION_NR < 1060000910)
typedef boost::optional
#else
typedef std::optional
#endif
<
    Tree::Intersection_and_primitive_id<Segment>::Type
> Segment_intersection;

#endif // NO_CGAL


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool intersectSurfaces
(
    triSurface& surf,
    edgeIntersections& edgeCuts
)
{
    bool hasMoved = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surface edges with itself" << endl;

        // Determine surface edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters = 0;
        {
            triSurfaceSearch querySurf(surf);

            scalarField surfPointTol
            (
                max(1e-3*edgeIntersections::minEdgeLength(surf), SMALL)
            );

            // Determine raw intersections
            edgeCuts = edgeIntersections
            (
                surf,
                querySurf,
                surfPointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points(surf.points());

                nIters =
                    edgeCuts.removeDegenerates
                    (
                        5,              // max iterations
                        surf,
                        querySurf,
                        surfPointTol,
                        points         // work array
                    );

                if (nIters != 0)
                {
                    // Update geometric quantities
                    surf.movePoints(points);
                    hasMoved = true;
                }
            }
        }
    }

    if (hasMoved)
    {
        fileName newFile("surf.obj");
        Info<< "Surface has been moved. Writing to " << newFile << endl;
        surf.write(newFile);
    }

    return hasMoved;
}


// Keep on shuffling surface points until no more degenerate intersections.
// Moves both surfaces and updates set of edge cuts.
bool intersectSurfaces
(
    triSurface& surf1,
    edgeIntersections& edgeCuts1,
    triSurface& surf2,
    edgeIntersections& edgeCuts2
)
{
    bool hasMoved1 = false;
    bool hasMoved2 = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surf1 edges with surf2"
            << " faces" << endl;

        // Determine surface1 edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters1 = 0;
        {
            triSurfaceSearch querySurf2(surf2);

            scalarField surf1PointTol
            (
                max(1e-3*edgeIntersections::minEdgeLength(surf1), SMALL)
            );

            // Determine raw intersections
            edgeCuts1 = edgeIntersections
            (
                surf1,
                querySurf2,
                surf1PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points1(surf1.points());

                nIters1 =
                    edgeCuts1.removeDegenerates
                    (
                        5,              // max iterations
                        surf1,
                        querySurf2,
                        surf1PointTol,
                        points1         // work array
                    );

                if (nIters1 != 0)
                {
                    // Update geometric quantities
                    surf1.movePoints(points1);
                    hasMoved1 = true;
                }
            }
        }

        Info<< "Determining intersections of surf2 edges with surf1"
            << " faces" << endl;

        label nIters2 = 0;
        {
            triSurfaceSearch querySurf1(surf1);

            scalarField surf2PointTol
            (
                max(1e-3*edgeIntersections::minEdgeLength(surf2), SMALL)
            );

            // Determine raw intersections
            edgeCuts2 = edgeIntersections
            (
                surf2,
                querySurf1,
                surf2PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points2(surf2.points());

                nIters2 =
                    edgeCuts2.removeDegenerates
                    (
                        5,              // max iterations
                        surf2,
                        querySurf1,
                        surf2PointTol,
                        points2         // work array
                    );

                if (nIters2 != 0)
                {
                    // Update geometric quantities
                    surf2.movePoints(points2);
                    hasMoved2 = true;
                }
            }
        }

        if (nIters1 == 0 && nIters2 == 0)
        {
            Info<< "** Resolved all intersections to be proper edge-face pierce"
                << endl;
            break;
        }
    }

    if (hasMoved1)
    {
        fileName newFile("surf1.obj");
        Info<< "Surface 1 has been moved. Writing to " << newFile
            << endl;
        surf1.write(newFile);
    }

    if (hasMoved2)
    {
        fileName newFile("surf2.obj");
        Info<< "Surface 2 has been moved. Writing to " << newFile
            << endl;
        surf2.write(newFile);
    }

    return hasMoved1 || hasMoved2;
}


label calcNormalDirection
(
    const vector& normal,
    const vector& otherNormal,
    const vector& edgeDir,
    const vector& faceCentre,
    const vector& pointOnEdge
)
{
    const vector cross = normalised(normal ^ edgeDir);

    const vector fC0tofE0 = normalised(faceCentre - pointOnEdge);

    label nDir = ((cross & fC0tofE0) > 0.0 ? 1 : -1);

    nDir *= ((otherNormal & fC0tofE0) > 0.0 ? -1 : 1);

    return nDir;
}


void calcEdgeCuts
(
    triSurface& surf1,
    triSurface& surf2,
    const bool perturb,
    edgeIntersections& edgeCuts1,
    edgeIntersections& edgeCuts2
)
{
    if (perturb)
    {
        intersectSurfaces
        (
            surf1,
            edgeCuts1,
            surf2,
            edgeCuts2
        );
    }
    else
    {
        triSurfaceSearch querySurf2(surf2);

        Info<< "Determining intersections of surf1 edges with surf2 faces"
            << endl;

        edgeCuts1 = edgeIntersections
        (
            surf1,
            querySurf2,
            max(1e-3*edgeIntersections::minEdgeLength(surf1), SMALL)
        );

        triSurfaceSearch querySurf1(surf1);

        Info<< "Determining intersections of surf2 edges with surf1 faces"
            << endl;

        edgeCuts2 = edgeIntersections
        (
            surf2,
            querySurf1,
            max(1e-3*edgeIntersections::minEdgeLength(surf2), SMALL)
        );
    }
}


// CGAL variants

#ifndef NO_CGAL

void visitPointRegion
(
    const triSurface& s,
    const label zoneI,
    const label pointI,
    const label startEdgeI,
    const label startFaceI,
    labelList& pFacesZone
)
{
    const labelList& eFaces = s.edgeFaces()[startEdgeI];

    if (eFaces.size() == 2)
    {
        label nextFaceI;
        if (eFaces[0] == startFaceI)
        {
            nextFaceI = eFaces[1];
        }
        else if (eFaces[1] == startFaceI)
        {
            nextFaceI = eFaces[0];
        }
        else
        {
            FatalErrorInFunction
                << "problem" << exit(FatalError);

            nextFaceI = -1;
        }



        label index = s.pointFaces()[pointI].find(nextFaceI);

        if (pFacesZone[index] == -1)
        {
            // Mark face as been visited.
            pFacesZone[index] = zoneI;

            // Step to next edge on face which is still using pointI
            const labelList& fEdges = s.faceEdges()[nextFaceI];

            label nextEdgeI = -1;

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];
                const edge& e = s.edges()[edgeI];

                if (edgeI != startEdgeI && (e[0] == pointI || e[1] == pointI))
                {
                    nextEdgeI = edgeI;

                    break;
                }
            }

            if (nextEdgeI == -1)
            {
                FatalErrorInFunction
                    << "Problem: cannot find edge out of " << fEdges
                    << "on face " << nextFaceI << " that uses point " << pointI
                    << " and is not edge " << startEdgeI << abort(FatalError);
            }


            visitPointRegion
            (
                s,
                zoneI,
                pointI,
                nextEdgeI,
                nextFaceI,
                pFacesZone
            );
        }
    }
}


label dupNonManifoldPoints(triSurface& s, labelList& pointMap)
{
    const labelListList& pf = s.pointFaces();
    const labelListList& fe = s.faceEdges();
    const edgeList& edges = s.edges();


    DynamicField<point> newPoints(s.points());
    // From dupSurf back to s.pointa
    DynamicList<label> newPointMap(identity(newPoints.size()));
    List<labelledTri> newFaces(s);
    label nNonManifold = 0;

    forAll(pf, pointI)
    {
        const labelList& pFaces = pf[pointI];

        // Visited faces (as indices into pFaces)
        labelList pFacesZone(pFaces.size(), -1);

        label nZones = 0;
        label index = 0;
        for (; index < pFacesZone.size(); index++)
        {
            if (pFacesZone[index] == -1)
            {
                label zoneI = nZones++;
                pFacesZone[index] = zoneI;

                label faceI = pFaces[index];
                const labelList& fEdges = fe[faceI];

                // Find starting edge
                forAll(fEdges, fEdgeI)
                {
                    label edgeI = fEdges[fEdgeI];
                    const edge& e = edges[edgeI];

                    if (e[0] == pointI || e[1] == pointI)
                    {
                        visitPointRegion
                        (
                            s,
                            zoneI,
                            pointI,
                            edgeI,
                            faceI,
                            pFacesZone
                        );
                    }
                }
            }
        }


        // Subset
        if (nZones > 1)
        {
            for (label zoneI = 1; zoneI < nZones; zoneI++)
            {
                label newPointI = newPoints.size();
                newPointMap.append(s.meshPoints()[pointI]);
                newPoints.append(s.points()[s.meshPoints()[pointI]]);

                forAll(pFacesZone, index)
                {
                    if (pFacesZone[index] == zoneI)
                    {
                        label faceI = pFaces[index];
                        const labelledTri& localF = s.localFaces()[faceI];
                        forAll(localF, fp)
                        {
                            if (localF[fp] == pointI)
                            {
                                newFaces[faceI][fp] = newPointI;
                            }
                        }
                    }
                }
            }
            nNonManifold++;
        }
    }


    Info<< "Duplicating " << nNonManifold << " points out of " << s.nPoints()
        << endl;
    if (nNonManifold > 0)
    {
        triSurface dupSurf(newFaces, s.patches(), newPoints, true);

        // Create map from dupSurf localPoints to s.localPoints
        const Map<label> mpm = s.meshPointMap();

        const labelList& dupMp = dupSurf.meshPoints();

        labelList dupPointMap(dupMp.size());
        forAll(dupMp, pointI)
        {
            label dupMeshPointI = dupMp[pointI];
            label meshPointI = newPointMap[dupMeshPointI];
            dupPointMap[pointI] = mpm[meshPointI];
        }


        forAll(dupPointMap, pointI)
        {
            const point& dupPt = dupSurf.points()[dupMp[pointI]];
            label sLocalPointI = dupPointMap[pointI];
            label sMeshPointI = s.meshPoints()[sLocalPointI];
            const point& sPt = s.points()[sMeshPointI];

            if (mag(dupPt-sPt) > SMALL)
            {
                FatalErrorInFunction
                    << "dupPt:" << dupPt
                    << " sPt:" << sPt
                    << exit(FatalError);
            }
        }


        //s.transfer(dupSurf);
        s = dupSurf;
        pointMap = labelUIndList(pointMap, dupPointMap)();
    }

    return nNonManifold;
}


// Find intersections of surf1 by edges of surf2. Return number of degenerate
// cuts (cuts through faces/edges/points)
labelPair edgeIntersectionsCGAL
(
    const Tree& tree,
    const triSurface& surf,
    const pointField& points,
    edgeIntersections& edgeCuts
)
{
    const edgeList& edges = surf.edges();
    const labelList& meshPoints = surf.meshPoints();

    //Info<< "Intersecting CGAL surface ..." << endl;
    List<List<pointIndexHit>> intersections(edges.size());
    labelListList classifications(edges.size());

    label nPoints = 0;
    label nSegments = 0;

    std::vector<Segment_intersection> segments;
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        const point& a = points[meshPoints[e[0]]];
        const point& b = points[meshPoints[e[1]]];

        K::Segment_3 segment_query
        (
            Point(a.x(), a.y(), a.z()),
            Point(b.x(), b.y(), b.z())
        );

        segments.clear();
        tree.all_intersections(segment_query, std::back_inserter(segments));


        for (const Segment_intersection& intersect : segments)
        {
            // Get intersection object
            if
            (
                // Used boost::variant prior to CGAL-6.0
                #if (CGAL_VERSION_NR < 1060000910)
                const Point* ptPtr = boost::get<Point>(&(intersect->first))
                #else
                const Point* ptPtr = std::get_if<Point>(&(intersect->first))
                #endif
            )
            {
                point pt
                (
                    CGAL::to_double(ptPtr->x()),
                    CGAL::to_double(ptPtr->y()),
                    CGAL::to_double(ptPtr->z())
                );

                #if defined (CGAL_VERSION_NR) && (CGAL_VERSION_NR < 1041400000)
                Polyhedron::Face_handle f = (intersect->second);
                #else
                // 1.14 and later
                Polyhedron::Face_handle f = (intersect->second).first;
                #endif

                intersections[edgeI].append
                (
                    pointIndexHit
                    (
                        true,
                        pt,
                        f->index
                    )
                );
                // Intersection on edge interior
                classifications[edgeI].append(-1);
                ++nPoints;
            }
            else if
            (
                // Used boost::variant prior to CGAL-6.0
                #if (CGAL_VERSION_NR < 1060000910)
                const Segment* sPtr = boost::get<Segment>(&(intersect->first))
                #else
                const Segment* sPtr = std::get_if<Segment>(&(intersect->first))
                #endif
            )
            {
                #if defined (CGAL_VERSION_NR) && (CGAL_VERSION_NR < 1041400000)
                Polyhedron::Face_handle f = (intersect->second);
                #else
                // 1.14 and later
                Polyhedron::Face_handle f = (intersect->second).first;
                #endif

                //std::cout
                //    << "intersection object is a segment:" << sPtr->source()
                //    << " " << sPtr->target() << std::endl;

                //std::cout<< "triangle:" << f->index
                //    << " region:" << f->region << std::endl;

                const point source
                (
                    CGAL::to_double(sPtr->source().x()),
                    CGAL::to_double(sPtr->source().y()),
                    CGAL::to_double(sPtr->source().z())
                );

                const point target
                (
                    CGAL::to_double(sPtr->target().x()),
                    CGAL::to_double(sPtr->target().y()),
                    CGAL::to_double(sPtr->target().z())
                );

                // Make up some intersection point
                intersections[edgeI].append
                (
                    pointIndexHit
                    (
                        true,
                        0.5*(source+target),
                        f->index
                    )
                );
                // Intersection aligned with face. Tbd: enums
                classifications[edgeI].append(2);
                ++nSegments;
            }
        }
    }

    edgeCuts = edgeIntersections(intersections, classifications);

    return labelPair(nPoints, nSegments);
}


// Intersect edges of surf1 until no more degenerate intersections. Return
// number of degenerates
labelPair edgeIntersectionsAndShuffleCGAL
(
    Random& rndGen,
    const triSurface& surf2,
    const scalarField& surf1PointTol,
    triSurface& surf1,
    edgeIntersections& edgeCuts1
)
{
    //Info<< "Constructing CGAL surface ..." << endl;
    Polyhedron p;
    PolyhedronReader(surf2, p);


    //Info<< "Constructing CGAL tree ..." << endl;
    const Tree tree(p.facets_begin(), p.facets_end(), p);


    labelPair cuts(0, 0);

    // Update surface1 points until no longer intersecting
    pointField surf1Points(surf1.points());

    for (label iter = 0; iter < 5; iter++)
    {
        // See which edges of 1 intersect 2
        Info<< "Determining intersections of " << surf1.nEdges()
            << " edges with surface of " << label(tree.size()) << " triangles"
            << endl;
        cuts = edgeIntersectionsCGAL
        (
            tree,
            surf1,
            surf1Points,
            edgeCuts1
        );
        Info<< "Determined intersections:" << nl
            << "    edges               : " << surf1.nEdges() << nl
            << "    non-degenerate cuts : " << cuts.first() << nl
            << "    degenerate cuts     : " << cuts.second() << nl
            << endl;

        if (cuts.second() == 0)
        {
            break;
        }

        Info<< "Shuffling conflicting points" << endl;

        const labelListList& edgeStat = edgeCuts1.classification();
        const edgeList& edges = surf1.edges();
        const labelList& mp = surf1.meshPoints();
        const point p05(0.5, 0.5, 0.5);

        forAll(edgeStat, edgeI)
        {
            const labelList& stat = edgeStat[edgeI];
            forAll(stat, i)
            {
                if (stat[i] == 2)
                {
                    const edge& e = edges[edgeI];
                    forAll(e, eI)
                    {
                        vector d = rndGen.sample01<vector>() - p05;
                        surf1Points[mp[e[eI]]] += surf1PointTol[e[eI]]*d;
                    }
                }
            }
        }
    }
    surf1.movePoints(surf1Points);

    return cuts;
}


// Return map from subSurf.edges() back to surf.edges()
labelList matchEdges
(
    const triSurface& surf,
    const triSurface& subSurf,
    const labelList& pointMap
)
{
    if (pointMap.size() != subSurf.nPoints())
    {
        FatalErrorInFunction
            << "problem" << exit(FatalError);
    }

    labelList edgeMap(subSurf.nEdges(), -1);

    const edgeList& edges = surf.edges();
    const labelListList& pointEdges = surf.pointEdges();

    const edgeList& subEdges = subSurf.edges();


    forAll(subEdges, subEdgeI)
    {
        const edge& subE = subEdges[subEdgeI];

        // Match points on edge to those on surf
        label start = pointMap[subE[0]];
        label end = pointMap[subE[1]];
        const labelList& pEdges = pointEdges[start];
        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];
            const edge& e = edges[edgeI];

            if (e.otherVertex(start) == end)
            {
                if (edgeMap[subEdgeI] == -1)
                {
                    edgeMap[subEdgeI] = edgeI;
                }
                else if (edgeMap[subEdgeI] != edgeI)
                {
                    FatalErrorInFunction
                        << subE << " points:"
                        << subE.line(subSurf.localPoints())
                        << " matches to " << edgeI
                        << " points:" << e.line(surf.localPoints())
                        << " but also " << edgeMap[subEdgeI]
                        << " points:"
                        << edges[edgeMap[subEdgeI]].line(surf.localPoints())
                        << exit(FatalError);
                }
                break;
            }
        }

        if (edgeMap[subEdgeI] == -1)
        {
            FatalErrorInFunction
                << subE << " at:" << subSurf.localPoints()[subE[0]]
                << subSurf.localPoints()[subE[1]]
                << exit(FatalError);
        }
    }

    return edgeMap;
}


void calcEdgeCutsCGAL
(
    triSurface& surf1,
    triSurface& surf2,
    const bool perturb,
    edgeIntersections& edgeCuts1,
    edgeIntersections& edgeCuts2
)
{
    if (!perturb)
    {
        // See which edges of 1 intersect 2
        {
            Info<< "Intersect surface 1 edges with surface 2:" << nl;
            Info<< "    constructing CGAL surface ..." << endl;
            Polyhedron p;
            PolyhedronReader(surf2, p);

            Info<< "    constructing CGAL tree ..." << endl;
            const Tree tree(p.facets_begin(), p.facets_end(), p);

            edgeIntersectionsCGAL
            (
                tree,
                surf1,
                surf1.points(),
                edgeCuts1
            );
        }
        // See which edges of 2 intersect 1
        {
            Info<< "Intersect surface 2 edges with surface 1:" << nl;
            Info<< "    constructing CGAL surface ..." << endl;
            Polyhedron p;
            PolyhedronReader(surf1, p);

            Info<< "    constructing CGAL tree ..." << endl;
            const Tree tree(p.facets_begin(), p.facets_end(), p);

            edgeIntersectionsCGAL
            (
                tree,
                surf2,
                surf2.points(),
                edgeCuts2
            );
        }
        Info<< endl;
    }
    else
    {
        const scalarField surf1PointTol
        (
            max(1e-8*edgeIntersections::minEdgeLength(surf1), SMALL)
        );
        const scalarField surf2PointTol
        (
            max(1e-8*edgeIntersections::minEdgeLength(surf2), SMALL)
        );


        Random rndGen(0);

        labelPair cuts1;
        labelPair cuts2;

        for (label iter = 0; iter < 10; iter++)
        {
            // Find intersections of surf1 edges with surf2 triangles
            cuts1 = edgeIntersectionsAndShuffleCGAL
            (
                rndGen,
                surf2,
                surf1PointTol,
                surf1,
                edgeCuts1
            );

            // Find intersections of surf2 edges with surf1 triangles
            cuts2 = edgeIntersectionsAndShuffleCGAL
            (
                rndGen,
                surf1,
                surf2PointTol,
                surf2,
                edgeCuts2
            );

            if (cuts1.second() + cuts2.second() == 0)
            {
                break;
            }
        }
    }
}


void calcEdgeCutsBitsCGAL
(
    triSurface& surf1,
    triSurface& surf2,
    const bool perturb,
    edgeIntersections& edgeCuts1,
    edgeIntersections& edgeCuts2
)
{
    label nZones1 = 1;
    labelList zone1;
    {
        labelHashSet orientationEdge(surf1.size()/1000);
        PatchTools::checkOrientation(surf1, false, &orientationEdge);
        nZones1 = PatchTools::markZones(surf1, orientationEdge, zone1);

        Info<< "Surface triangles " << surf1.size()
            << " number of zones (orientation or non-manifold):"
            << nZones1 << endl;
    }

    label nZones2 = 1;
    labelList zone2;
    {
        labelHashSet orientationEdge(surf2.size()/1000);
        PatchTools::checkOrientation(surf2, false, &orientationEdge);
        nZones2 = PatchTools::markZones(surf2, orientationEdge, zone2);

        Info<< "Surface triangles " << surf2.size()
            << " number of zones (orientation or non-manifold):"
            << nZones2 << endl;
    }


    if (nZones1 == 1 && nZones2 == 1)
    {
        calcEdgeCutsCGAL
        (
            surf1,
            surf2,
            perturb,
            edgeCuts1,
            edgeCuts2
        );
    }
    else
    {
        edgeCuts1 = edgeIntersections
        (
            List<List<pointIndexHit>>(surf1.nEdges()),
            labelListList(surf1.nEdges())
        );
        edgeCuts2 = edgeIntersections
        (
            List<List<pointIndexHit>>(surf2.nEdges()),
            labelListList(surf2.nEdges())
        );


        for (label zone1I = 0; zone1I < nZones1; zone1I++)
        {
            // Generate sub surface for zone1I
            boolList includeMap1(surf1.size(), false);

            forAll(zone1, faceI)
            {
                if (zone1[faceI] == zone1I)
                {
                    includeMap1[faceI] = true;
                }
            }

            // Subset. Map from local points on subset to local points on
            // original
            labelList pointMap1;
            labelList faceMap1;
            triSurface subSurf1
            (
                surf1.subsetMesh
                (
                    includeMap1,
                    pointMap1,
                    faceMap1
                )
            );

            // Split non-manifold points; update pointMap
            dupNonManifoldPoints(subSurf1, pointMap1);

            const boundBox subBb1
            (
                subSurf1.points(),
                subSurf1.meshPoints(),
                false
            );

            const labelList edgeMap1
            (
                matchEdges
                (
                    surf1,
                    subSurf1,
                    pointMap1
                )
            );


            for (label zone2I = 0; zone2I < nZones2; zone2I++)
            {
                // Generate sub surface for zone2I
                boolList includeMap2(surf2.size(), false);

                forAll(zone2, faceI)
                {
                    if (zone2[faceI] == zone2I)
                    {
                        includeMap2[faceI] = true;
                    }
                }

                labelList pointMap2;
                labelList faceMap2;
                triSurface subSurf2
                (
                    surf2.subsetMesh
                    (
                        includeMap2,
                        pointMap2,
                        faceMap2
                    )
                );


                const boundBox subBb2
                (
                    subSurf2.points(),
                    subSurf2.meshPoints(),
                    false
                );

                // Short-circuit expensive calculations
                if (!subBb2.overlaps(subBb1))
                {
                    continue;
                }


                // Split non-manifold points; update pointMap
                dupNonManifoldPoints(subSurf2, pointMap2);

                const labelList edgeMap2
                (
                    matchEdges
                    (
                        surf2,
                        subSurf2,
                        pointMap2
                    )
                );


                // Do cuts
                edgeIntersections subEdgeCuts1;
                edgeIntersections subEdgeCuts2;
                calcEdgeCutsCGAL
                (
                    subSurf1,
                    subSurf2,
                    perturb,
                    subEdgeCuts1,
                    subEdgeCuts2
                );

                // Move original surface
                {
                    pointField points2(surf2.points());
                    forAll(pointMap2, i)
                    {
                        label subMeshPointI = subSurf2.meshPoints()[i];
                        label meshPointI = surf2.meshPoints()[pointMap2[i]];
                        points2[meshPointI] = subSurf2.points()[subMeshPointI];
                    }
                    surf2.movePoints(points2);
                }

                // Insert into main structure
                edgeCuts1.merge(subEdgeCuts1, edgeMap1, faceMap2);
                edgeCuts2.merge(subEdgeCuts2, edgeMap2, faceMap1);
            }


            // Move original surface
            {
                pointField points1(surf1.points());
                forAll(pointMap1, i)
                {
                    label subMeshPointI = subSurf1.meshPoints()[i];
                    label meshPointI = surf1.meshPoints()[pointMap1[i]];
                    points1[meshPointI] = subSurf1.points()[subMeshPointI];
                }
                surf1.movePoints(points1);
            }
        }
    }
}

#endif // NO_CGAL


//void calcFeaturePoints(const pointField& points, const edgeList& edges)
//{
//    edgeMesh eMesh(points, edges);
//
//    const labelListList& pointEdges = eMesh.pointEdges();
//
//
//    // Get total number of feature points
//    label nFeaturePoints = 0;
//    forAll(pointEdges, pI)
//    {
//        const labelList& pEdges = pointEdges[pI];
//
//        if (pEdges.size() == 1)
//        {
//            nFeaturePoints++;
//        }
//    }
//
//
//    // Calculate addressing from feature point to cut point and cut edge
//    labelList featurePointToCutPoint(nFeaturePoints);
//    labelList featurePointToCutEdge(nFeaturePoints);
//
//    label nFeatPts = 0;
//    forAll(pointEdges, pI)
//    {
//        const labelList& pEdges = pointEdges[pI];
//
//        if (pEdges.size() == 1)
//        {
//            featurePointToCutPoint[nFeatPts] = pI;
//            featurePointToCutEdge[nFeatPts] = pEdges[0];
//            nFeatPts++;
//        }
//    }
//
//
//
//    label concaveStart = 0;
//    label mixedStart = 0;
//    label nonFeatureStart = nFeaturePoints;
//
//
//    labelListList featurePointNormals(nFeaturePoints);
//    labelListList featurePointEdges(nFeaturePoints);
//    labelList regionEdges;
//}


autoPtr<extendedFeatureEdgeMesh> createEdgeMesh
(
    const IOobject& io,
    const booleanSurface::booleanOpType action,
    const bool surf1Baffle,
    const bool surf2Baffle,
    const bool invertedSpace,
    const triSurface& surf1,
    const edgeIntersections& edgeCuts1,
    const triSurface& surf2,
    const edgeIntersections& edgeCuts2
)
{
    // Determine intersection edges from the edge cuts
    surfaceIntersection inter
    (
        surf1,
        edgeCuts1,
        surf2,
        edgeCuts2
    );

    label nFeatEds = inter.cutEdges().size();

    DynamicList<vector> normals(2*nFeatEds);
    vectorField edgeDirections(nFeatEds, Zero);
    DynamicList<extendedFeatureEdgeMesh::sideVolumeType> normalVolumeTypes
    (
        2*nFeatEds
    );
    List<DynamicList<label>> edgeNormals(nFeatEds);
    List<DynamicList<label>> normalDirections(nFeatEds);


    const triSurface& s1 = surf1;
    const triSurface& s2 = surf2;

    forAllConstIters(inter.facePairToEdgeId(), iter)
    {
        const labelPair& facePair = iter.key();
        const label cutEdgeI = iter.val();

        const edge& fE = inter.cutEdges()[cutEdgeI];

        const vector& norm1 = s1.faceNormals()[facePair.first()];
        const vector& norm2 = s2.faceNormals()[facePair.second()];

        DynamicList<label>& eNormals = edgeNormals[cutEdgeI];
        DynamicList<label>& nDirections = normalDirections[cutEdgeI];

        edgeDirections[cutEdgeI] = fE.vec(inter.cutPoints());

        normals.append(norm1);
        eNormals.append(normals.size() - 1);

        if (surf1Baffle)
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

            nDirections.append(1);
        }
        else
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);
            nDirections.append
            (
                calcNormalDirection
                (
                    norm1,
                    norm2,
                    edgeDirections[cutEdgeI],
                    s1[facePair.first()].centre(s1.points()),
                    inter.cutPoints()[fE.start()]
                )
            );
        }

        normals.append(norm2);
        eNormals.append(normals.size() - 1);

        if (surf2Baffle)
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

            nDirections.append(1);
        }
        else
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

            nDirections.append
            (
                calcNormalDirection
                (
                    norm2,
                    norm1,
                    edgeDirections[cutEdgeI],
                    s2[facePair.second()].centre(s2.points()),
                    inter.cutPoints()[fE.start()]
                )
            );
        }


        if (surf1Baffle)
        {
            normals.append(norm2);

            if (surf2Baffle)
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

                nDirections.append(1);
            }
            else
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

                nDirections.append
                (
                    calcNormalDirection
                    (
                        norm2,
                        norm1,
                        edgeDirections[cutEdgeI],
                        s2[facePair.second()].centre(s2.points()),
                        inter.cutPoints()[fE.start()]
                    )
                );
            }

            eNormals.append(normals.size() - 1);
        }

        if (surf2Baffle)
        {
            normals.append(norm1);

            if (surf1Baffle)
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

                nDirections.append(1);
            }
            else
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

                nDirections.append
                (
                    calcNormalDirection
                    (
                        norm1,
                        norm2,
                        edgeDirections[cutEdgeI],
                        s1[facePair.first()].centre(s1.points()),
                        inter.cutPoints()[fE.start()]
                    )
                );
            }

            eNormals.append(normals.size() - 1);
        }
    }


    label internalStart = -1;
    label nIntOrExt = 0;
    label nFlat = 0;
    label nOpen = 0;
    // label nMultiple = 0;

    forAll(edgeNormals, eI)
    {
        label nEdNorms = edgeNormals[eI].size();

        if (nEdNorms == 1)
        {
            nOpen++;
        }
        else if (nEdNorms == 2)
        {
            const vector& n0(normals[edgeNormals[eI][0]]);
            const vector& n1(normals[edgeNormals[eI][1]]);

            if ((n0 & n1) > extendedFeatureEdgeMesh::cosNormalAngleTol_)
            {
                nFlat++;
            }
            else
            {
                nIntOrExt++;
            }
        }
        // else if (nEdNorms > 2)
        // {
        //     ++nMultiple;
        // }
    }

    if (action == booleanSurface::UNION)
    {
        if (!invertedSpace)
        {
            // All edges are internal
            internalStart = 0;
        }
        else
        {
            // All edges are external
            internalStart = nIntOrExt;
        }
    }
    else if (action == booleanSurface::INTERSECTION)
    {
        if (!invertedSpace)
        {
            // All edges are external
            internalStart = nIntOrExt;
        }
        else
        {
            // All edges are internal
            internalStart = 0;
        }
    }
    else if (action == booleanSurface::DIFFERENCE)
    {
        // All edges are external
        internalStart = nIntOrExt;
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported booleanSurface:booleanOpType and space "
            << action << " " << invertedSpace
            << abort(FatalError);
    }

    // There are no feature points supported by surfaceIntersection
    // Flat, open or multiple edges are assumed to be impossible
    // Region edges are not explicitly supported by surfaceIntersection

    vectorField normalsTmp(normals);
    List<extendedFeatureEdgeMesh::sideVolumeType> normalVolumeTypesTmp
    (
        normalVolumeTypes
    );
    labelListList edgeNormalsTmp(edgeNormals.size());
    forAll(edgeNormalsTmp, i)
    {
        edgeNormalsTmp[i] = edgeNormals[i];
    }
    labelListList normalDirectionsTmp(normalDirections.size());
    forAll(normalDirectionsTmp, i)
    {
        normalDirectionsTmp[i] = normalDirections[i];
    }

    //calcFeaturePoints(inter.cutPoints(), inter.cutEdges());

    return autoPtr<extendedFeatureEdgeMesh>::New
    (
        io,
        inter.cutPoints(),
        inter.cutEdges(),

        0,                  // concaveStart,
        0,                  // mixedStart,
        0,                  // nonFeatureStart,

        internalStart,      // internalStart,
        nIntOrExt,          // flatStart,
        nIntOrExt + nFlat,  // openStart,
        nIntOrExt + nFlat + nOpen,   // multipleStart,

        normalsTmp,
        normalVolumeTypesTmp,
        edgeDirections,
        normalDirectionsTmp,
        edgeNormalsTmp,

        labelListList(),    // featurePointNormals,
        labelListList(),    // featurePointEdges,
        labelList()         // regionEdges
    );
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generates a extendedFeatureEdgeMesh for the interface created by"
        " a boolean operation on two surfaces."
        #ifndef NO_CGAL
        " [Compiled with CGAL]"
        #else
        " [Compiled without CGAL]"
        #endif
    );

    argList::noParallel();
    argList::addArgument
    (
        "action",
        "One of (intersection | union | difference)"
    );
    argList::addArgument("surface1", "The input surface file 1");
    argList::addArgument("surface2", "The input surface file 2");

    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor (both surfaces)"
    );
    argList::addBoolOption
    (
        "surf1Baffle",
        "Mark surface 1 as a baffle"
    );

    argList::addBoolOption
    (
        "surf2Baffle",
        "Mark surface 2 as a baffle"
    );

    argList::addBoolOption
    (
        "perturb",
        "Perturb surface points to escape degenerate intersections"
    );

    argList::addBoolOption
    (
        "no-cgal",
        #ifndef NO_CGAL
        "Do not use CGAL algorithms"
        #else
        "Ignored, compiled without CGAL"
        #endif
    );

    argList::addBoolOption
    (
        "invertedSpace",
        "Do the surfaces have inverted space orientation, "
        "i.e. a point at infinity is considered inside. "
        "This is only sensible for union and intersection."
    );

    argList::addOption
    (
        "trim",
        "((surface1 volumeType) .. (surfaceN volumeType))",
        "Trim resulting intersection with additional surfaces;"
        " volumeType is 'inside' (keep (parts of) edges that are inside)"
        ", 'outside' (keep (parts of) edges that are outside) or"
        " 'mixed' (keep all)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const word action(args[1]);

    const Enum<booleanSurface::booleanOpType> validActions
    {
        { booleanSurface::INTERSECTION, "intersection" },
        { booleanSurface::UNION, "union" },
        { booleanSurface::DIFFERENCE, "difference" }
    };

    if (!validActions.found(action))
    {
        FatalErrorInFunction
            << "Unsupported action " << action << endl
            << "Supported actions:" << validActions << nl
            << abort(FatalError);
    }


    List<Pair<word>> surfaceAndSide;
    if (args.readIfPresent("trim", surfaceAndSide))
    {
        Info<< "Trimming edges with " << surfaceAndSide << endl;
    }


    // Scale factor for both surfaces:
    const scalar scaleFactor = args.getOrDefault<scalar>("scale", -1);

    const word surf1Name(args[2]);
    Info<< "Reading surface " << surf1Name << endl;
    triSurfaceMesh surf1
    (
        IOobject
        (
            surf1Name,
            runTime.constant(),
            triSurfaceMesh::meshSubDir,
            runTime
        )
    );
    if (scaleFactor > 0)
    {
        Info<< "Scaling  : " << scaleFactor << nl;
        surf1.scalePoints(scaleFactor);
    }

    Info<< surf1Name << " statistics:" << endl;
    surf1.writeStats(Info);
    Info<< endl;

    const word surf2Name(args[3]);
    Info<< "Reading surface " << surf2Name << endl;
    triSurfaceMesh surf2
    (
        IOobject
        (
            surf2Name,
            runTime.constant(),
            triSurfaceMesh::meshSubDir,
            runTime
        )
    );
    if (scaleFactor > 0)
    {
        Info<< "Scaling  : " << scaleFactor << nl;
        surf2.scalePoints(scaleFactor);
    }

    Info<< surf2Name << " statistics:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    const bool surf1Baffle = args.found("surf1Baffle");
    const bool surf2Baffle = args.found("surf2Baffle");

    edgeIntersections edgeCuts1;
    edgeIntersections edgeCuts2;

    const bool invertedSpace = args.found("invertedSpace");

    if (invertedSpace && validActions[action] == booleanSurface::DIFFERENCE)
    {
        FatalErrorInFunction
            << "Inverted space only makes sense for union or intersection."
            << exit(FatalError);
    }


    // Calculate where edges are cut by the other surface
#ifndef NO_CGAL
    if (!args.found("no-cgal"))
    {
        calcEdgeCutsBitsCGAL
        (
            surf1,
            surf2,
            args.found("perturb"),
            edgeCuts1,
            edgeCuts2
        );
    }
    else
#endif // NO_CGAL
    {
        calcEdgeCuts
        (
            surf1,
            surf2,
            args.found("perturb"),
            edgeCuts1,
            edgeCuts2
        );
    }

    const fileName sFeatFileName
    (
        fileName::stem(surf1Name)
      + "_"
      + fileName::stem(surf2Name)
      + "_"
      + action
    );

    autoPtr<extendedFeatureEdgeMesh> feMeshPtr
    (
        createEdgeMesh
        (
            IOobject
            (
                sFeatFileName + ".extendedFeatureEdgeMesh",
                runTime.constant(),
                "extendedFeatureEdgeMesh",
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            booleanSurface::booleanOpTypeNames[action],
            surf1Baffle,
            surf2Baffle,
            invertedSpace,
            surf1,
            edgeCuts1,
            surf2,
            edgeCuts2
        )
    );


    // Trim intersections
    forAll(surfaceAndSide, trimI)
    {
        const word& trimName = surfaceAndSide[trimI].first();
        const volumeType trimType
        (
            volumeType::names[surfaceAndSide[trimI].second()]
        );

        Info<< "Reading trim surface " << trimName << endl;
        const triSurfaceMesh trimSurf
        (
            IOobject
            (
                trimName,
                runTime.constant(),
                triSurfaceMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        Info<< trimName << " statistics:" << endl;
        trimSurf.writeStats(Info);
        Info<< endl;

        labelList pointMap;
        labelList edgeMap;
        feMeshPtr().trim
        (
            trimSurf,
            trimType,
            pointMap,
            edgeMap
        );
    }


    const extendedFeatureEdgeMesh& feMesh = feMeshPtr();

    feMesh.writeStats(Info);
    Info<< nl << "Writing extendedFeatureEdgeMesh to "
        << feMesh.objectRelPath() << nl
        << endl;
    feMesh.write();
    feMesh.writeObj(feMesh.path()/sFeatFileName);

    {
        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                sFeatFileName + ".eMesh",   // name
                runTime.constant(),         // instance
                triSurfaceMesh::meshSubDir,
                runTime,                    // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectRelPath() << endl;

        bfeMesh.regIOobject::write();
    }

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
