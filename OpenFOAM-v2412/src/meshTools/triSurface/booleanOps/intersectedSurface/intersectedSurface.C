/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "intersectedSurface.H"
#include "surfaceIntersection.H"
#include "faceList.H"
#include "faceTriangulation.H"
#include "treeBoundBox.H"
#include "OFstream.H"
#include "error.H"
#include "meshTools.H"
#include "edgeSurface.H"
#include "DynamicList.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(intersectedSurface, 0);
}


const Foam::label Foam::intersectedSurface::UNVISITED = 0;
const Foam::label Foam::intersectedSurface::STARTTOEND = 1;
const Foam::label Foam::intersectedSurface::ENDTOSTART = 2;
const Foam::label Foam::intersectedSurface::BOTH = STARTTOEND | ENDTOSTART;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::intersectedSurface::writeOBJ
(
    const pointField& points,
    const edgeList& edges,
    Ostream& os
)
{
    for (const point& pt : points)
    {
        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    for (const edge& e : edges)
    {
        os << "l " << e.start()+1 << ' ' << e.end()+1 << nl;
    }
}


void Foam::intersectedSurface::writeOBJ
(
    const pointField& points,
    const edgeList& edges,
    const labelList& faceEdges,
    Ostream& os
)
{
    for (const point& pt : points)
    {
        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }
    for (const label edgei : faceEdges)
    {
        const edge& e = edges[edgei];

        os << "l " << e.start()+1 << ' ' << e.end()+1 << nl;
    }
}


void Foam::intersectedSurface::writeLocalOBJ
(
    const pointField& points,
    const edgeList& edges,
    const labelList& faceEdges,
    const fileName& fName
)
{
    OFstream os(fName);

    labelList pointMap(points.size(), -1);

    label maxVerti = 0;

    for (const label edgei : faceEdges)
    {
        const edge& e = edges[edgei];

        forAll(e, i)
        {
            const label pointi = e[i];

            if (pointMap[pointi] == -1)
            {
                const point& pt = points[pointi];

                os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;

                pointMap[pointi] = maxVerti++;
            }
        }
    }

    for (const label edgei : faceEdges)
    {
        const edge& e = edges[edgei];

        os  << "l " << pointMap[e.start()]+1 << ' ' << pointMap[e.end()]+1
            << nl;
    }
}


void Foam::intersectedSurface::writeOBJ
(
    const pointField& points,
    const face& f,
    Ostream& os
)
{
    for (const point& pt : points)
    {
        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    os << 'f';

    for (const label pointi : f)
    {
        os << ' ' << pointi+1;
    }
    os << nl;
}


void Foam::intersectedSurface::printVisit
(
    const edgeList& edges,
    const labelList& edgeLabels,
    const Map<label>& visited
)
{
    Pout<< "Visited:" << nl;
    for (const label edgei : edgeLabels)
    {
        const edge& e = edges[edgei];

        label stat = visited[edgei];

        if (stat == UNVISITED)
        {
            Pout<< "    edge:" << edgei << "  verts:" << e
                << "  unvisited" << nl;
        }
        else if (stat == STARTTOEND)
        {
            Pout<< "    edge:" << edgei << "  from " << e[0]
                << " to " << e[1] << nl;
        }
        else if (stat == ENDTOSTART)
        {
            Pout<< "    edge:" << edgei << "  from " << e[1]
                << " to " << e[0] << nl;
        }
        else
        {
            Pout<< "    edge:" << edgei << "  both " << e
                << nl;
        }
    }
    Pout<< endl;
}


bool Foam::intersectedSurface::sameEdgeOrder
(
    const labelledTri& fA,
    const labelledTri& fB
)
{
    forAll(fA, fpA)
    {
        label fpB = fB.find(fA[fpA]);

        if (fpB != -1)
        {
            // Get prev/next vertex on fA
            label vA1 = fA[fA.fcIndex(fpA)];
            label vAMin1 = fA[fA.rcIndex(fpA)];

            // Get prev/next vertex on fB
            label vB1 = fB[fB.fcIndex(fpB)];
            label vBMin1 = fB[fB.rcIndex(fpB)];

            if (vA1 == vB1 || vAMin1 == vBMin1)
            {
                return true;
            }
            else if (vA1 == vBMin1 || vAMin1 == vB1)
            {
                // shared vertices in opposite order.
                return false;
            }
            else
            {
                FatalErrorInFunction
                    << "Triangle:" << fA << " and triangle:" << fB
                    << " share a point but not an edge"
                    << abort(FatalError);
            }
        }
    }

    FatalErrorInFunction
        << "Triangle:" << fA << " and triangle:" << fB
        << " do not share an edge"
        << abort(FatalError);

    return false;
}


void Foam::intersectedSurface::incCount
(
    Map<label>& visited,
    const label key,
    const label offset
)
{
    visited(key, 0) += offset;
}


Foam::Map<Foam::DynamicList<Foam::label>>
Foam::intersectedSurface::calcPointEdgeAddressing
(
    const edgeSurface& eSurf,
    const label facei
)
{
    const pointField& points = eSurf.points();
    const edgeList& edges = eSurf.edges();

    const labelList& fEdges = eSurf.faceEdges()[facei];

    Map<DynamicList<label>> facePointEdges(4*fEdges.size());

    for (const label edgei : fEdges)
    {
        const edge& e = edges[edgei];

        // Add e.start to point-edges
        facePointEdges(e.start()).append(edgei);

        // Add e.end to point-edges
        facePointEdges(e.end()).append(edgei);
    }

    // Shrink it
    forAllIters(facePointEdges, iter)
    {
        iter.val().shrink();

        // Check on dangling points.
        if (iter.val().empty())
        {
            FatalErrorInFunction
                << "Point:" << iter.key() << " used by too few edges:"
                << iter.val() << abort(FatalError);
        }
    }

    if (debug & 2)
    {
        // Print facePointEdges
        Pout<< "calcPointEdgeAddressing: face consisting of edges:" << nl;
        for (const label edgei : fEdges)
        {
            const edge& e = edges[edgei];
            Pout<< "    " << edgei << ' ' << e
                << points[e.start()]
                << points[e.end()] << nl;
        }

        Pout<< "    Constructed point-edge addressing:" << nl;
        forAllConstIters(facePointEdges, iter)
        {
            Pout<< "    vertex " << iter.key() << " is connected to edges "
                << iter.val() << nl;
        }
        Pout<< endl;
    }

    return facePointEdges;
}


Foam::label Foam::intersectedSurface::nextEdge
(
    const edgeSurface& eSurf,
    const Map<label>& visited,
    const label facei,
    const vector& n,
    const Map<DynamicList<label>>& facePointEdges,
    const label prevEdgei,
    const label prevVerti
)
{
    const pointField& points = eSurf.points();
    const edgeList& edges = eSurf.edges();
    const labelList& fEdges = eSurf.faceEdges()[facei];


    // Edges connected to prevVerti
    const DynamicList<label>& connectedEdges = facePointEdges[prevVerti];

    if (connectedEdges.size() <= 1)
    {
        // Problem. Point not connected.
        {
            Pout<< "Writing face:" << facei << " to face.obj" << endl;
            OFstream str("face.obj");
            writeOBJ(points, edges, fEdges, str);

            Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
            writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
        }

        FatalErrorInFunction
            << "Problem: prevVertI:" << prevVerti << " on edge " << prevEdgei
            << " has less than 2 connected edges."
            << " connectedEdges:" << connectedEdges << abort(FatalError);

        return -1;
    }

    if (connectedEdges.size() == 2)
    {
        // Simple case. Take other edge
        if (connectedEdges[0] == prevEdgei)
        {
            if (debug & 2)
            {
                Pout<< "Stepped from edge:" << edges[prevEdgei]
                    << " to edge:" << edges[connectedEdges[1]]
                    << " chosen from candidates:" << connectedEdges << endl;
            }
            return connectedEdges[1];
        }
        else
        {
            if (debug & 2)
            {
               Pout<< "Stepped from edge:" << edges[prevEdgei]
                   << " to edge:" << edges[connectedEdges[0]]
                   << " chosen from candidates:" << connectedEdges << endl;
            }
            return connectedEdges[0];
        }
    }


    // Multiple choices. Look at angle between edges.

    const edge& prevE = edges[prevEdgei];

    // x-axis of coordinate system
    const vector e0 =
        normalised
        (
            n ^ (points[prevE.otherVertex(prevVerti)] - points[prevVerti])
        );

    // y-axis of coordinate system
    const vector e1 = n ^ e0;

    if (mag(mag(e1) - 1) > SMALL)
    {
        {
            Pout<< "Writing face:" << facei << " to face.obj" << endl;
            OFstream str("face.obj");
            writeOBJ(points, edges, fEdges, str);

            Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
            writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
        }

        FatalErrorInFunction
            << "Unnormalized normal e1:" << e1
            << " formed from cross product of e0:" << e0 << " n:" << n
            << abort(FatalError);
    }


    //
    // Determine maximum angle over all connected (unvisited) edges.
    //

    scalar maxAngle = -GREAT;
    label maxEdgei = -1;

    for (const label edgei : connectedEdges)
    {
        if (edgei != prevEdgei)
        {
            label stat = visited[edgei];

            const edge& e = edges[edgei];

            // Find out whether walk of edge from prevVert would be acceptable.
            if
            (
                stat == UNVISITED
             || (
                    stat == STARTTOEND
                 && prevVerti == e[1]
                )
             || (
                    stat == ENDTOSTART
                 && prevVerti == e[0]
                )
            )
            {
                // Calculate angle of edge with respect to base e0, e1
                const vector vec =
                    normalised
                    (
                        n
                      ^ (points[e.otherVertex(prevVerti)] - points[prevVerti])
                    );

                scalar angle = pseudoAngle(e0, e1, vec);

                if (angle > maxAngle)
                {
                    maxAngle = angle;
                    maxEdgei = edgei;
                }
            }
        }
    }


    if (maxEdgei == -1)
    {
        // No unvisited edge found
        {
            Pout<< "Writing face:" << facei << " to face.obj" << endl;
            OFstream str("face.obj");
            writeOBJ(points, edges, fEdges, str);

            Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
            writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
        }

        FatalErrorInFunction
            << "Trying to step from edge " << edges[prevEdgei]
            << ", vertex " << prevVerti
            << " but cannot find 'unvisited' edges among candidates:"
            << connectedEdges
            << abort(FatalError);
    }

    if (debug & 2)
    {
        Pout<< "Stepped from edge:" << edges[prevEdgei]
            << " to edge:" << maxEdgei << " angle:" << edges[maxEdgei]
            << " with angle:" << maxAngle
            << endl;
    }

    return maxEdgei;
}


Foam::face Foam::intersectedSurface::walkFace
(
    const edgeSurface& eSurf,
    const label facei,
    const vector& n,
    const Map<DynamicList<label>>& facePointEdges,

    const label startEdgei,
    const label startVerti,

    Map<label>& visited
)
{
    const pointField& points = eSurf.points();
    const edgeList& edges = eSurf.edges();

    // Overestimate size of face
    face f(eSurf.faceEdges()[facei].size());

    label fp = 0;

    label verti = startVerti;
    label edgei = startEdgei;

    while (true)
    {
        const edge& e = edges[edgei];

        if (debug & 2)
        {
            Pout<< "Now at:" << endl
                << "    edge:" << edgei << " vertices:" << e
                << " positions:" << points[e.start()] << ' ' << points[e.end()]
                << "    vertex:" << verti << endl;
        }

        // Mark edge as visited
        if (e[0] == verti)
        {
            visited[edgei] |= STARTTOEND;
        }
        else
        {
            visited[edgei] |= ENDTOSTART;
        }


        // Store face vertex
        f[fp++] = verti;

        // Step to other vertex
        verti = e.otherVertex(verti);

        if (verti == startVerti)
        {
            break;
        }

        // Step from vertex to next edge
        edgei = nextEdge
        (
            eSurf,
            visited,
            facei,
            n,
            facePointEdges,
            edgei,
            verti
        );
    }

    f.setSize(fp);

    return f;
}


void Foam::intersectedSurface::findNearestVisited
(
    const edgeSurface& eSurf,
    const label facei,
    const Map<DynamicList<label>>& facePointEdges,
    const Map<label>& pointVisited,
    const point& pt,
    const label excludePointi,

    label& minVerti,
    scalar& minDist
)
{
    minVerti = -1;
    minDist = GREAT;

    forAllConstIters(pointVisited, iter)
    {
        const label pointi = iter.key();
        const label nVisits = iter.val();

        if (pointi != excludePointi)
        {
            if (nVisits == 2*facePointEdges[pointi].size())
            {
                // Fully visited (i.e. both sides of all edges)
                scalar dist = mag(eSurf.points()[pointi] - pt);

                if (dist < minDist)
                {
                    minDist = dist;
                    minVerti = pointi;
                }
            }
        }
    }

    if (minVerti == -1)
    {
        const labelList& fEdges = eSurf.faceEdges()[facei];

        SeriousErrorInFunction
            << "Dumping face edges to faceEdges.obj" << endl;

        writeLocalOBJ(eSurf.points(), eSurf.edges(), fEdges, "faceEdges.obj");

        FatalErrorInFunction
            << "No fully visited edge found for pt " << pt
            << abort(FatalError);
    }
}


Foam::faceList Foam::intersectedSurface::resplitFace
(
    const triSurface& surf,
    const label facei,
    const Map<DynamicList<label>>& facePointEdges,
    const Map<label>& visited,
    edgeSurface& eSurf
)
{
    // Count the number of times point has been visited so we
    // can compare number to facePointEdges.
    Map<label> pointVisited(2*facePointEdges.size());

    forAllConstIters(visited, iter)
    {
        const label edgei = iter.key();
        const label stat = iter.val();

        const edge& e = eSurf.edges()[edgei];

        if (stat == STARTTOEND || stat == ENDTOSTART)
        {
            incCount(pointVisited, e[0], 1);
            incCount(pointVisited, e[1], 1);
        }
        else if (stat == BOTH)
        {
            incCount(pointVisited, e[0], 2);
            incCount(pointVisited, e[1], 2);
        }
        else if (stat == UNVISITED)
        {
            incCount(pointVisited, e[0], 0);
            incCount(pointVisited, e[1], 0);
        }
    }


    if (debug)
    {
        forAllConstIters(pointVisited, iter)
        {
            const label pointi = iter.key();
            const label nVisits = iter.val();

            Pout<< "point:" << pointi
                << "  nVisited:" << nVisits
                << "  pointEdges:" << facePointEdges[pointi].size() << endl;
        }
    }


    // Find nearest point pair where one is not fully visited and
    // the other is.

    label visitedVert0 = -1;
    label unvisitedVert0 = -1;

    {
        scalar minDist = GREAT;

        forAllConstIters(facePointEdges, iter)
        {
            const label pointi = iter.key();

            const DynamicList<label>& pEdges = iter.val();

            const label nVisits = pointVisited[pointi];

            if (nVisits < 2*pEdges.size())
            {
                // Not fully visited. Find nearest fully visited.

                scalar nearDist;
                label nearVerti;

                findNearestVisited
                (
                    eSurf,
                    facei,
                    facePointEdges,
                    pointVisited,
                    eSurf.points()[pointi],
                    -1,                         // Do not exclude vertex
                    nearVerti,
                    nearDist
                );


                if (nearDist < minDist)
                {
                    minDist = nearDist;
                    visitedVert0 = nearVerti;
                    unvisitedVert0 = pointi;
                }
            }
        }
    }


    // Find second intersection.
    {
        scalar minDist = GREAT;

        forAllConstIters(facePointEdges, iter)
        {
            const label pointi = iter.key();

            const DynamicList<label>& pEdges = iter.val();

            if (pointi != unvisitedVert0)
            {
                const label nVisits = pointVisited[pointi];

                if (nVisits < 2*pEdges.size())
                {
                    // Not fully visited. Find nearest fully visited.

                    scalar nearDist;
                    label nearVerti;

                    findNearestVisited
                    (
                        eSurf,
                        facei,
                        facePointEdges,
                        pointVisited,
                        eSurf.points()[pointi],
                        visitedVert0,           // vertex to exclude
                        nearVerti,
                        nearDist
                    );


                    if (nearDist < minDist)
                    {
                        minDist = nearDist;
                    }
                }
            }
        }
    }


    // Add the new intersection edges to the edgeSurface
    edgeList additionalEdges(1);
    additionalEdges[0] = edge(visitedVert0, unvisitedVert0);

    eSurf.addIntersectionEdges(facei, additionalEdges);

    if (debug)
    {
        fileName newFName("face_" + Foam::name(facei) + "_newEdges.obj");
        Pout<< "Dumping face:" << facei << " to " << newFName << endl;
        writeLocalOBJ
        (
            eSurf.points(),
            eSurf.edges(),
            eSurf.faceEdges()[facei],
            newFName
        );
    }

    // Retry splitFace. Use recursion since is rare situation.
    return splitFace(surf, facei, eSurf);
}


Foam::faceList Foam::intersectedSurface::splitFace
(
    const triSurface& surf,
    const label facei,
    edgeSurface& eSurf
)
{
    // Alias
    const pointField& points = eSurf.points();
    const edgeList& edges = eSurf.edges();
    const labelList& fEdges = eSurf.faceEdges()[facei];

    // Create local (for the face only) point-edge connectivity.
    Map<DynamicList<label>> facePointEdges
    (
        calcPointEdgeAddressing
        (
            eSurf,
            facei
        )
    );

    // Order in which edges have been walked. Initialize outside edges.
    Map<label> visited(fEdges.size()*2);

    forAll(fEdges, i)
    {
        label edgei = fEdges[i];

        if (eSurf.isSurfaceEdge(edgei))
        {
            // Edge is edge from original surface so an outside edge for
            // the current face.
            label surfEdgei = eSurf.parentEdge(edgei);

            label owner = surf.edgeOwner()[surfEdgei];

            if
            (
                owner == facei
             || sameEdgeOrder
                (
                    surf.localFaces()[owner],
                    surf.localFaces()[facei]
                )
            )
            {
                // Edge is in same order as current face.
                // Mark off the opposite order.
                visited.insert(edgei, ENDTOSTART);
            }
            else
            {
                // Edge is in same order as current face.
                // Mark off the opposite order.
                visited.insert(edgei, STARTTOEND);
            }
        }
        else
        {
            visited.insert(edgei, UNVISITED);
        }
    }



    // Storage for faces
    DynamicList<face> faces(fEdges.size());

    while (true)
    {
        // Find starting edge:
        // - unvisited triangle edge
        // - once visited intersection edge
        // Give priority to triangle edges.
        label startEdgei = -1;
        label startVerti = -1;

        forAll(fEdges, i)
        {
            label edgei = fEdges[i];

            const edge& e = edges[edgei];

            label stat = visited[edgei];

            if (stat == STARTTOEND)
            {
                startEdgei = edgei;
                startVerti = e[1];

                if (eSurf.isSurfaceEdge(edgei))
                {
                    break;
                }
            }
            else if (stat == ENDTOSTART)
            {
                startEdgei = edgei;
                startVerti = e[0];

                if (eSurf.isSurfaceEdge(edgei))
                {
                    break;
                }
            }
        }

        if (startEdgei == -1)
        {
            break;
        }

        //Pout<< "splitFace: starting walk from edge:" << startEdgei
        //    << ' ' << edges[startEdgei] << " vertex:" << startVerti << endl;

        //// Print current visited state.
        //printVisit(eSurf.edges(), fEdges, visited);

        //{
        //    Pout<< "Writing face:" << faceI << " to face.obj" << endl;
        //    OFstream str("face.obj");
        //    writeOBJ(eSurf.points(), eSurf.edges(), fEdges, str);
        //}

        faces.append
        (
            walkFace
            (
                eSurf,
                facei,
                surf.faceNormals()[facei],
                facePointEdges,

                startEdgei,
                startVerti,

                visited
            )
        );
    }


    // Check if any unvisited edges left.
    forAll(fEdges, i)
    {
        label edgei = fEdges[i];

        label stat = visited[edgei];

        if (eSurf.isSurfaceEdge(edgei) && stat != BOTH)
        {
            SeriousErrorInFunction
                << "Dumping face edges to faceEdges.obj" << endl;

            writeLocalOBJ(points, edges, fEdges, "faceEdges.obj");

            FatalErrorInFunction
               << "Problem: edge " << edgei << " vertices "
                << edges[edgei] << " on face " << facei
                << " has visited status " << stat << " from a "
                 << "right-handed walk along all"
                << " of the triangle edges. Are the original surfaces"
                << " closed and non-intersecting?"
                << abort(FatalError);
        }
        else if (stat != BOTH)
        {
            // Redo face after introducing extra edge. Edge introduced
            // should be one nearest to any fully visited edge.
            return resplitFace
            (
                surf,
                facei,
                facePointEdges,
                visited,
                eSurf
            );
        }
    }


    // See if normal needs flipping.
    faces.shrink();

    const vector n = faces[0].areaNormal(eSurf.points());

    if ((n & surf.faceNormals()[facei]) < 0)
    {
        for (face& f : faces)
        {
            reverse(f);
        }
    }

    return faces;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::intersectedSurface::intersectedSurface()
:
    triSurface(),
    intersectionEdges_(0),
    faceMap_(0),
    nSurfacePoints_(0)
{}


Foam::intersectedSurface::intersectedSurface(const triSurface& surf)
:
    triSurface(surf),
    intersectionEdges_(0),
    faceMap_(0),
    nSurfacePoints_(0)
{}


Foam::intersectedSurface::intersectedSurface
(
    const triSurface& surf,
    const bool isFirstSurface,
    const surfaceIntersection& inter
)
:
    triSurface(),
    intersectionEdges_(0),
    faceMap_(0),
    nSurfacePoints_(surf.nPoints())
{
    if (inter.cutPoints().empty() && inter.cutEdges().empty())
    {
        // No intersection. Make straight copy.
        triSurface::operator=(surf);

        // Identity for face map
        faceMap_.setSize(size());

        forAll(faceMap_, facei)
        {
            faceMap_[facei] = facei;
        }
        return;
    }


    // Calculate face-edge addressing for surface + intersection.
    edgeSurface eSurf(surf, isFirstSurface, inter);

    // Update point information for any removed points. (when are they removed?
    // - but better make sure)
    nSurfacePoints_ = eSurf.nSurfacePoints();

    // Now we have full points, edges and edge addressing for surf. Start
    // extracting faces and triangulate them.

    // Storage for new triangles of surface.
    DynamicList<labelledTri> newTris(eSurf.edges().size()/2);

    // Start in newTris for decomposed face.
    labelList startTrii(surf.size(), Zero);

    forAll(surf, facei)
    {
        startTrii[facei] = newTris.size();

        if (eSurf.faceEdges()[facei].size() != surf.faceEdges()[facei].size())
        {
            // Face has been cut by intersection.
            // Cut face into multiple subfaces. Use faceEdge information
            // from edgeSurface only. (original triSurface 'surf' is used for
            // faceNormals and region number only)
            faceList newFaces
            (
                splitFace
                (
                    surf,
                    facei,              // current triangle
                    eSurf               // face-edge description of surface
                                        // + intersection
                )
            );
            forAll(newFaces, newFacei)
            {
                const face& newF = newFaces[newFacei];
                const vector& n = surf.faceNormals()[facei];
                const label region = surf[facei].region();

                faceTriangulation tris(eSurf.points(), newF, n);

                forAll(tris, trii)
                {
                    const triFace& t = tris[trii];

                    forAll(t, i)
                    {
                        if (t[i] < 0 || t[i] >= eSurf.points().size())
                        {
                            FatalErrorInFunction
                                << "Face triangulation of face " << facei
                                << " uses points outside range 0.."
                                << eSurf.points().size()-1 << endl
                                << "Triangulation:"
                                << tris << abort(FatalError);
                        }
                    }

                    newTris.append(labelledTri(t[0], t[1], t[2], region));
                }
            }
        }
        else
        {
            // Face has not been cut at all. No need to renumber vertices since
            // eSurf keeps surface vertices first.
            newTris.append(surf.localFaces()[facei]);
        }
    }

    newTris.shrink();


    // Construct triSurface. Note that addressing will be compact since
    // edgeSurface is compact.
    triSurface::operator=
    (
        triSurface
        (
            newTris,
            surf.patches(),
            eSurf.points()
        )
    );

    // Construct mapping back into original surface
    faceMap_.setSize(size());

    for (label facei = 0; facei < surf.size()-1; facei++)
    {
        for (label trii = startTrii[facei]; trii < startTrii[facei+1]; trii++)
        {
            faceMap_[trii] = facei;
        }
    }
    for (label trii = startTrii[surf.size()-1]; trii < size(); trii++)
    {
        faceMap_[trii] = surf.size()-1;
    }


    // Find edges on *this which originate from 'cuts'. (i.e. newEdgei >=
    // nSurfaceEdges) Renumber edges into local triSurface numbering.

    intersectionEdges_.setSize(eSurf.edges().size() - eSurf.nSurfaceEdges());

    label intersectionEdgei = 0;

    for
    (
        label edgei = eSurf.nSurfaceEdges();
        edgei < eSurf.edges().size();
        edgei++
    )
    {
        // Get edge vertices in triSurface local numbering
        const edge& e = eSurf.edges()[edgei];
        label surfStarti = meshPointMap()[e.start()];
        label surfEndi = meshPointMap()[e.end()];

        // Find edge connected to surfStarti which also uses surfEndi.
        label surfEdgei = -1;

        const labelList& pEdges = pointEdges()[surfStarti];

        forAll(pEdges, i)
        {
            const edge& surfE = edges()[pEdges[i]];

            // Edge already connected to surfStart for sure. See if also
            // connects to surfEnd
            if (surfE.found(surfEndi))
            {
                surfEdgei = pEdges[i];

                break;
            }
        }

        if (surfEdgei != -1)
        {
            intersectionEdges_[intersectionEdgei++] = surfEdgei;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find edge among candidates " << pEdges
                << " which uses points " << surfStarti
                << " and " << surfEndi
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
