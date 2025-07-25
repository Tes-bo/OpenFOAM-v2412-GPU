/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::triSurfaceTools

Description
    A collection of tools for triSurface.

Note
    The curvature calculation is an implementation of the algorithm from:
    \verbatim
        "Estimating Curvatures and their Derivatives on Triangle Meshes"
        by S. Rusinkiewicz
        3DPVT'04 Proceedings of the 3D Data Processing,
        Visualization, and Transmission, 2nd International Symposium
        Pages 486-493
        http://gfx.cs.princeton.edu/pubs/_2004_ECA/curvpaper.pdf
    \endverbatim

SourceFiles
    triSurfaceTools.C
    triSurfaceCloseness.C
    triSurfaceCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceTools_H
#define triSurfaceTools_H

#include "boolList.H"
#include "pointField.H"
#include "vectorField.H"
#include "triadFieldFwd.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "Map.H"
#include "FixedList.H"
#include "Pair.H"
#include "vector2D.H"
#include "triangle.H"
#include "surfaceLocation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBox;
class edge;
class labelledTri;
class polyBoundaryMesh;
class plane;
class triSurface;
class face;
class Time;
template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                           Class triSurfaceTools Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceTools
{
    // Private Member Functions

        // Refinement

            enum refineType
            {
                NONE,
                RED,
                GREEN
            };
            static void calcRefineStatus
            (
                const triSurface& surf,
                const label facei,
                List<refineType>& refine
            );
            static void greenRefine
            (
                const triSurface& surf,
                const label facei,
                const label edgeI,
                const label newPointi,
                DynamicList<labelledTri>& newFaces
            );
            static triSurface doRefine
            (
                const triSurface& surf,
                const List<refineType>& refineStatus
            );


        // Coarsening

            static scalar faceCosAngle
            (
                const point& pStart,
                const point& pEnd,
                const point& pLeft,
                const point& pRight
            );

            static void protectNeighbours
            (
                const triSurface& surf,
                const label vertI,
                labelList& faceStatus
            );

            //- Faces to collapse because of edge collapse
            static labelHashSet getCollapsedFaces
            (
                const triSurface& surf,
                label edgeI
            );

            // Return value of faceUsed for faces using vertI (local numbering).
            // Used internally.
            static label vertexUsesFace
            (
                const triSurface& surf,
                const labelHashSet& faceUsed,
                const label vertI
            );

            // Get new connections between faces (because of edge collapse)
            // in form of tables:
            //  - given edge get other edge
            //  - given edge get other face
            // A face using point v1 on edge will get connected to a face using
            //  point v2 if they share a common vertex
            //  (but not a common edge since then the triangles collapse to
            //  nothing)
            static void getMergedEdges
            (
                const triSurface& surf,
                const label edgeI,
                const labelHashSet& collapsedFaces,
                Map<label>& edgeToEdge,
                Map<label>& edgeToFace
            );

            //- Calculates (cos of) angle across edgeI of facei,
            //  taking into account updated addressing (resulting from edge
            //  collapse)
            static scalar edgeCosAngle
            (
                const triSurface& surf,
                const label v1,
                const point& pt,
                const labelHashSet& collapsedFaces,
                const Map<label>& edgeToEdge,
                const Map<label>& edgeToFace,
                const label facei,
                const label edgeI
            );

            //- Calculate minimum (cos of) edge angle using addressing from
            //  collapsing
            //  edge to v1 at pt. Returns 1 if v1 is on edge without neighbours
            //  (and hence no edge angle can be defined)
            static scalar collapseMinCosAngle
            (
                const triSurface& surf,
                const label v1,
                const point& pt,
                const labelHashSet& collapsedFaces,
                const Map<label>& edgeToEdge,
                const Map<label>& edgeToFace
            );

            //- Like collapseMinCosAngle but return true for value < minCos
            bool collapseCreatesFold
            (
                const triSurface& surf,
                const label v1,
                const point& pt,
                const labelHashSet& collapsedFaces,
                const Map<label>& edgeToEdge,
                const Map<label>& edgeToFace,
                const scalar minCos
            );

            ////- Checks if edge collapse creates triangles on top of each
            ////  other
            //static bool collapseCreatesDuplicates
            //(
            //    const triSurface& surf,
            //    const label edgeI,
            //    const Map<bool>& collapsedFaces
            //);

        // Tracking

            //- Finds the triangle edge/point cut by the plane between
            //  a point inside/on edge of a triangle and a point outside.
            //  Returns
            //  - location on edge/point and hit()
            //  - or miss() if no intersection found
            static surfaceLocation cutEdge
            (
                const triSurface& s,
                const label triI,
                const label excludeEdgeI,
                const label excludePointi,
                const point& triPoint,
                const plane& cutPlane,
                const point& toPoint
            );

            //- Checks if current is on the same triangle as the endpoint
            //  and shifts it there. If so updates current and sets a hit.
            static void snapToEnd
            (
                const triSurface& s,
                const surfaceLocation& endInfo,
                surfaceLocation& current
            );

            //- Visits faces eFaces around start. Does not visit triangle
            //  start.triangle() nor edge excludeEdgeI.
            //  Returns edge, triangle (if more than one choice) which gets
            //  us nearer endpoint.
            //  Returns
            //  - hit() if triangle contains endpoint
            //  - triangle()=-1 if no triangle found
            //  - nearest triangle/edge otherwise
            static surfaceLocation visitFaces
            (
                const triSurface& s,
                const labelList& eFaces,
                const surfaceLocation& start,
                const label excludeEdgeI,
                const label excludePointi,
                const surfaceLocation& end,
                const plane& cutPlane
            );


public:

    // OBJ writing

        //- Write pointField to OBJ format file
        static void writeOBJ
        (
            const fileName& fName,
            const pointField& pts
        );

        //- Write vertex subset to OBJ format file
        static void writeOBJ
        (
            const triSurface& surf,
            const fileName& fName,
            const boolList& markedVerts
        );


    // Additional addressing

        //- Get all triangles using edge endpoint
        static void getVertexTriangles
        (
            const triSurface& surf,
            const label edgeI,
            labelList& edgeTris
        );

        //- Get all vertices (local numbering) connected to vertices of edge
        static labelList getVertexVertices
        (
            const triSurface& surf,
            const edge& e
        );

        //- Get face connected to edge not facei
        static label otherFace
        (
            const triSurface& surf,
            const label facei,
            const label edgeI
        );

        //- Get the two edges on facei counterclockwise after edgeI
        static void otherEdges
        (
            const triSurface& surf,
            const label facei,
            const label edgeI,
            label& e1,
            label& e2
        );

        //- Get the two vertices (local numbering) on facei counterclockwise
        //  vertI
        static void otherVertices
        (
            const triSurface& surf,
            const label facei,
            const label vertI,
            label& vert1I,
            label& vert2I
        );

        //- Get edge opposite vertex (local numbering)
        static label oppositeEdge
        (
            const triSurface& surf,
            const label facei,
            const label vertI
        );

        //- Get vertex (local numbering) opposite edge
        static label oppositeVertex
        (
            const triSurface& surf,
            const label facei,
            const label edgeI
        );

        //- Returns edge label connecting v1, v2 (local numbering)
        static label getEdge
        (
            const triSurface& surf,
            const label vert1I,
            const label vert2I
        );

        //- Return index of triangle (or -1) using all three edges
        static label getTriangle
        (
            const triSurface& surf,
            const label e0I,
            const label e1I,
            const label e2I
        );

    // Coarsening

        //- Create new triSurface by collapsing edges to edge mids.
        static triSurface collapseEdges
        (
            const triSurface& surf,
            const labelList& collapsableEdges
        );


        //- Face collapse status.
        //  anyEdge: any edge can be collapsed
        //  noEdge: no edge can be collapsed
        //  collapsed: already collapsed
        //  >0: edge label that can be collapsed
        static const label ANYEDGE;
        static const label NOEDGE;
        static const label COLLAPSED;

        //- Create new triSurface by collapsing edges to specified
        //  positions. faceStatus allows
        //  explicit control over which faces need to be protected (see above).
        //  faceStatus gets updated to protect collapsing already collapsed
        //  faces.
        static triSurface collapseEdges
        (
            const triSurface& surf,
            const labelList& collapsableEdges,
            const pointField& edgeMids,
            labelList& faceStatus
        );


    // Refinement

        //- Refine edges by splitting to opposite vertex
        static triSurface greenRefine
        (
            const triSurface& surf,
            const labelList& refineEdges
        );

        //- Refine face by splitting all edges. Neighbouring face is
        //  greenRefine'd.
        static triSurface redGreenRefine
        (
            const triSurface& surf,
            const labelList& refineFaces
        );


    // Geometric

        //- Returns element in edgeIndices with minimum length
        static label minEdge
        (
            const triSurface& surf,
            const labelList& edgeIndices
        );

        //- Returns element in edgeIndices with minimum length
        static label maxEdge
        (
            const triSurface& surf,
            const labelList& edgeIndices
        );

        //- Merge points within distance
        static triSurface mergePoints
        (
            const triSurface& surf,
            const scalar mergeTol
        );

        //- Triangle (unit) normal. If nearest point to triangle on edge use
        //  edge normal (calculated on the fly); if on vertex use vertex normal.
        //  Uses planarTol.
        static vector surfaceNormal
        (
            const triSurface& surf,
            const label nearestFacei,
            const point& nearestPt
        );

        //- On which side of surface
        enum sideType
        {
            UNKNOWN,    // cannot be determined (e.g. non-manifold)
            INSIDE,     // inside
            OUTSIDE     // outside
        };

        //- If nearest point is on edgeI, determine on which side of surface
        //  sample is.
        static sideType edgeSide
        (
            const triSurface& surf,
            const point& sample,
            const point& nearestPoint,
            const label edgeI
        );

        //- Given nearest point (to sample) on surface determines which side
        //  sample is. Uses either face normal, edge normal or point normal
        //  (non-trivial). Uses triangle::classify.
        static sideType surfaceSide
        (
            const triSurface& surf,
            const point& sample,
            const label nearestFacei
        );

    // Triangulation of faces

        //- Simple triangulation of (selected patches of) boundaryMesh. Needs
        //  polyMesh (or polyBoundaryMesh) since only at this level are the
        //  triangles on neighbouring patches connected.
        //  Return faceMap from triI to faceI
        static triSurface triangulate
        (
            const polyBoundaryMesh& mBesh,
            const labelHashSet& includePatches,
            labelList& faceMap,
            const bool verbose = false
        );


        static triSurface triangulate
        (
            const polyBoundaryMesh& bMesh,
            const labelHashSet& includePatches,
            const boundBox& bBox,
            const bool verbose = false
        );


        //- Face-centre triangulation of (selected patches of) boundaryMesh.
        //  Needs
        //  polyMesh (or polyBoundaryMesh) since only at this level are the
        //  triangles on neighbouring patches connected.
        triSurface triangulateFaceCentre
        (
            const polyBoundaryMesh& mBesh,
            const labelHashSet& includePatches,
            const bool verbose = false
        );


    // Triangulation and interpolation

        //- Do unconstrained Delaunay of points. Returns triSurface with 3D
        //  points with z=0. All triangles in region 0.
        static triSurface delaunay2D(const List<vector2D>&);

        //- Calculate linear interpolation weights for point (guaranteed to be
        //  inside triangle)
        static void calcInterpolationWeights
        (
            const triPointRef& tri,
            const point& p,
            FixedList<scalar, 3>& weights
        );

        // Calculate weighting factors from samplePts to triangle it is in.
        // Uses linear search to find triangle.
        // Vertices are:
        //   (a b c)  : vertices of the triangle abc the point is in
        // or if the point is outside all triangles:
        //   (a b -1) : the edge ab the point is nearest to.
        //   (a -1 -1) : the vertex a the point is nearest to
        static void calcInterpolationWeights
        (
            const triSurface& s,
            const pointField& samplePts,
            List<FixedList<label, 3>>& allVerts,
            List<FixedList<scalar, 3>>& allWeights
        );


    // Curvature

        //- Weighting for normals of faces attached to vertex
        static scalar vertexNormalWeight
        (
            const triFace& f,
            const label pI,
            const vector& fN,
            const UList<point>& points
        );

        //- Weighted average of normals of attached faces
        static tmp<vectorField> vertexNormals(const triSurface& surf);

        //- Local coordinate-system for each point normal
        static tmp<triadField> vertexTriads
        (
            const triSurface& surf,
            const vectorField& pointNormals
        );

        //- Surface curvatures at the vertex points
        static tmp<scalarField> curvatures
        (
            const triSurface& surf,
            const vectorField& pointNormals,
            const triadField& pointTriads
        );

        //- Surface curvatures at the vertex points
        static tmp<scalarField> curvatures
        (
            const triSurface& surf
        );

        //- Write surface curvature at the vertex points and return the field
        static tmp<scalarField> writeCurvature
        (
            const Time& runTime,
            const word& basename,
            const triSurface& surf
        );


    // Closeness

        //- Check and write internal/external closeness fields
        static Pair<tmp<scalarField>> writeCloseness
        (
            const Time& runTime,
            const word& basename,
            const triSurface& surf,
            const scalar internalAngleTolerance = 45,
            const scalar externalAngleTolerance = 10
        );


    // Feature Proximity

        //- Calculate feature proximity
        scalarField featureProximity
        (
            const triSurface& surf,
            const scalar searchDistance
        );

        //- Check and write internal/external closeness fields
        static void writeFeatureProximity
        (
            const Time& runTime,
            const word& basename,
            const triSurface& surf,
            const bool writeVTK,
            const scalar searchDistance
        );


    // Surface checking functionality

        //- Check single triangle for (topological) validity
        static bool validTri
        (
            const triSurface&,
            const label facei,
            const bool verbose = true
        );

        //- Check single triangle for (topological) validity
        static bool validTri
        (
            const MeshedSurface<face>&,
            const label facei,
            const bool verbose = true
        );


    // Tracking

        //- Test point on plane of triangle to see if on edge or point or inside
        static surfaceLocation classify
        (
            const triSurface&,
            const label triI,
            const point& trianglePoint
        );

        //- Track on surface to get closer to point.
        //  Possible situations:
        //  - 1. reached endpoint
        //  - 2. reached edge (normal situation)
        //  - 3. reached end of surface (edge on single face)
        //  Input:
        //  - starting position+triangle/edge/point (so has to be on surface!)
        //  - (optional) previous position+triangle/edge/point to prevent
        //    going back. Set index (of triangle/edge/point) to -1 if not
        //    used.
        //  - end position+triangle/edge/point (so has to be on surface!)
        //  - plane to follow. Has to go through end point!
        //  Returns:
        //  - true if end point reached (situation 1)
        //  - new position+triangle/edge/point
        //  Caller has to check for situation 3 by checking that triangle()
        //  is not set.
        static surfaceLocation trackToEdge
        (
            const triSurface&,
            const surfaceLocation& start,
            const surfaceLocation& end,
            const plane& cutPlane
        );

        //- Track from edge to edge across surface. Uses trackToEdge.
        //  Not really useful by itself, more example of how to use trackToEdge.
        //  endInfo should be location on surface.
        //  hitInfo should be initialised to starting location (on surface as
        //  well). Upon return is set to end location.
        static void track
        (
            const triSurface&,
            const surfaceLocation& endInfo,
            const plane& cutPlane,
            surfaceLocation& hitInfo
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
