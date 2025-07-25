/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::surfaceIntersection

Description
    Basic surface-surface intersection description. Constructed from two
    surfaces it creates a description of the intersection.

    The intersection information consists of the intersection line(s)
    with new points, new edges between points (note that these edges and
    points are on both surfaces) and various addressing from original
    surface faces/edges to intersection and vice versa.

    Gets either precalculated intersection information or calculates it
    itself.
    Algorithm works by intersecting all edges of one surface with the other
    surface and storing a reference from both faces (one on surface1, one on
    surface 2) to the vertex. If the reference re-occurs we have the second
    hit of both faces and an edge is created between the retrieved vertex and
    the new one.

    Note: when doing intersecting itself uses 'tolerance' as a fraction of
    current edge length to determine if intersection is a point-touching one
    instead of an edge-piercing action.

    Some constructors allow a dictionary of intersection controls:
    \table
        Property       | Description                    | Type   | Default value
        tolerance      | Edge-length tolerance          | scalar | 1e-3
        allowEdgeHits  | Edge-end cuts another edge     | bool   | true
        snap           | Snap near end-points           | bool   | true
        warnDegenerate | Number of warnings about degenerate edges | label | 0
        intersectionMethod | Control for self intersection | (self,region,none)
    \endtable

SourceFiles
    surfaceIntersection.C
    surfaceIntersectionFuncs.C
    surfaceIntersectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIntersection_H
#define surfaceIntersection_H

#include "DynamicList.H"
#include "point.H"
#include "edgeHashes.H"
#include "edgeList.H"
#include "labelPairHashes.H"
#include "pointIndexHit.H"
#include "typeInfo.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurfaceSearch;
class triSurface;
class edgeIntersections;


//- Key is non-commutative pair of labels. Value is commutative pair of labels
typedef LabelPairMap<edge> labelPairEdgeLookup;

//- Map from edge back to all parents (pairs of faces)
typedef EdgeMap<labelPairHashSet> edgelabelPairHashLookup;


/*---------------------------------------------------------------------------*\
                     Class surfaceIntersection Declaration
\*---------------------------------------------------------------------------*/

class surfaceIntersection
{
public:

    //- Surface intersection types for classify, doCutEdges
    enum intersectionType
    {
        FIRST,          //!< First surface
        SECOND,         //!< Second surface
        SELF,           //!< Self-intersection
        SELF_REGION,    //!< Self-intersection, region-wise only
        NONE            //!< None = invalid (for input only)
    };

    //- The user-selectable self-intersection enumeration names
    static const Enum<intersectionType> selfIntersectionNames;


private:

    // Private data

        //- Tolerance for intersections
        scalar tolerance_;

        //- Allow edge-ends to cut another edge.
        bool allowEdgeHits_;

        //- Snap cut points near edge ends (default: true)
        bool snapToEnd_;

        //- Maximum number of warnings about degenerate edges
        label warnDegenerate_;

        //- Newly introduced points.
        pointField cutPoints_;

        //- Newly introduced edges (are on both surfaces).
        //  Reference into cutPoints.
        edgeList cutEdges_;

        //- From face on surf1 and face on surf2 to intersection edge
        labelPairEdgeLookup facePairToEdge_;

        //- From face on surf1 and face on surf2 to intersection edgeId
        //  (label in cutEdges)
        labelPairLookup facePairToEdgeId_;

        //- Edges on surf1 that are cut.
        //  From edgeId on surf1 to location in cutPoint
        labelListList surf1EdgeCuts_;

        //- Edges on surf2 that are cut.
        //  From edgeId on surf2 to location in cutPoint
        labelListList surf2EdgeCuts_;

        //- Temporary storage to manage edge-edge self-intersections.
        edgeHashSet edgeEdgeIntersection_;

        //- Temporary storage to manage cuts/intersections from the edge ends
        Map<label> snappedEnds_;

        //- Temporary storage
        EdgeMap<label> edgeToId_;


    // Private Member Functions

        //- Adjust intersection options according to the dictionary entries
        void setOptions(const dictionary& dict);

        //- Transfer contents of List<DynamicList<..>> to List<List<..>>
        template<class T>
        static void transfer(List<DynamicList<T>>&,  List<List<T>>&);

        //- Get minimum length of all edges connected to point
        static scalar minEdgeLen(const triSurface& surf, const label pointi);

        //- Get edge label of edge between face vertices fp and fp+1
        static label getEdge
        (
            const triSurface& surf,
            const label facei,
            const label fp
        );

        //- Remove duplicates from ordered dynamic list. Returns map from old
        //  to new (-1 if element removed)
        static void removeDuplicates(const labelList& map, labelList& labels);

        // Remove all duplicate and degenerate elements. Return unique elements
        // and map from old to new.
        static edgeList filterEdges(const edgeList&, labelList& map);

        //- Remove all duplicate elements.
        static labelList filterLabels(const labelList& elems, labelList& map);

        //- Debugging: Dump intersected edges to stream
        void writeIntersectedEdges
        (
            const triSurface& surf,
            const labelListList& edgeCutVerts,
            Ostream& os
        ) const;

        //- Detect if point close to edge of end. Returns -1: not close.
        //  0:close (within startTol) to start, 1:close (within endTol) to end
        static label classify
        (
            const scalar startTol,
            const scalar endTol,
            const point& p,
            const edge& e,
            const UList<point>& points
        );

        //- Update reference between faceA and faceB.
        //  Updates facePairToEdge_ and facePairToEdgeId_ (on the second hit)
        void storeIntersection
        (
            const enum intersectionType cutFrom,
            const labelList& facesA,
            const label faceB,
            const UList<point>& allCutPoints,
            const label cutPointId,
            DynamicList<edge>& allCutEdges
        );

        //- Investigate pHit to whether is case of point hits point,
        //  point hits edge, point hits face or edge hits face.
        void classifyHit
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
        );

        //- Cut edges of surf1 with surface 2.
        void doCutEdges
        (
            const triSurface& surf1,
            const triSurfaceSearch& querySurf2,
            const enum intersectionType cutFrom,

            DynamicList<point>& allCutPoints,
            DynamicList<edge>& allCutEdges,
            List<DynamicList<label>>& surfEdgeCuts
        );

        //- Join disconnected intersection points
        void joinDisconnected(DynamicList<edge>& allCutEdges);


public:

    // Public Data, Declarations

    ClassName("surfaceIntersection");


    // Constructors

        //- Construct null
        surfaceIntersection();

        //- Construct from two surfaces.
        //  Does its own cutting. Has problems with degenerate cuts
        surfaceIntersection
        (
            const triSurfaceSearch& querySurf1,
            const triSurfaceSearch& querySurf2,
            const dictionary& dict = dictionary::null
        );

        //- Construct from self-intersections.
        //  Does its own cutting, but has problems with degenerate cuts
        surfaceIntersection
        (
            const triSurfaceSearch& querySurf1,
            const dictionary& dict = dictionary::null
        );

        //- Construct from precalculated intersection information.
        //  Advantage: intersection information is guaranteed to have no
        //  degenerate cuts.
        surfaceIntersection
        (
            const triSurface& surf1,
            const edgeIntersections& intersections1,
            const triSurface& surf2,
            const edgeIntersections& intersections2
        );


    // Member Functions

        //- The list of cut points
        const pointField& cutPoints() const;

        //- The list of created edges
        const edgeList& cutEdges() const;

        //- Lookup of pairs of faces to created edges
        const labelPairLookup& facePairToEdgeId() const;

        //- Access either surf1EdgeCuts (isFirstSurface = true) or
        //  surf2EdgeCuts
        const labelListList& edgeCuts(const bool isFirstSurf) const;

        //- List of cut points on edges of surface1
        const labelListList& surf1EdgeCuts() const;

        //- List of cut points on edges of surface2
        const labelListList& surf2EdgeCuts() const;


        //- Geometric merge points (points within mergeDist) prior to
        //  automatically calling mergeEdges().
        void mergePoints(const scalar mergeDist);

        //- Merge duplicate edges
        void mergeEdges();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceIntersectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
