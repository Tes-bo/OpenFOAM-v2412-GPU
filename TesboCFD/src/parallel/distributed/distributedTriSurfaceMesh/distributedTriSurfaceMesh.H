/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::distributedTriSurfaceMesh

Description
    IOoject and searching on distributed triSurface. All processor hold
    (possibly overlapping) part of the overall surface. All queries are
    distributed to the processor that can answer it and the result sent back.

    Can work in three modes:
    - follow : makes sure each processor has all the triangles inside the
    externally provided bounding box (usually the mesh bounding box).
    Guarantees minimum amount of communication since mesh-local queries
    should be answerable without any comms.
    - independent : surface is decomposed according to the triangle centres
    so the decomposition might be radically different from the mesh
    decomposition. Guarantees best memory balance but at the expense of
    more communication.
    - frozen : no change

    Note: addressing used:
    distributedTriSurfaceMesh: none

    triSurfaceMesh:
    - surf.pointFaces()     : edge addressing (for volume tests only)
    - surf.edges()          : edgeTree
    - surf.faceFaces()      : only if minQuality > 0


SourceFiles
    distributedTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_distributedTriSurfaceMesh_H
#define Foam_distributedTriSurfaceMesh_H

#include "triSurfaceMesh.H"
#include "localIOdictionary.H"
#include "IOdictionary.H"
#include "Pair.H"
#include "globalIndex.H"
#include "DynamicField.H"
#include "triangleFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;
class decompositionMethod;

// Typedefs
typedef Pair<point> segment;


/*---------------------------------------------------------------------------*\
                   Class distributedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class distributedTriSurfaceMesh
:
    public triSurfaceMesh
{
public:

    // Static data

        enum distributionType
        {
            FOLLOW = 0,
            INDEPENDENT = 1,
            DISTRIBUTED = 2,
            FROZEN = 3
        };

        static const Enum<distributionType> distributionTypeNames_;

private:

    // Private member data

        //- Merging distance
        scalar mergeDist_;

        autoPtr<IOdictionary> decomposeParDict_;

        //- Decomposition used when independently decomposing surface.
        autoPtr<decompositionMethod> decomposer_;

        //- Bounding box settings
        localIOdictionary dict_;

        //- Bounding boxes of all processors
        List<List<treeBoundBox>> procBb_;

        //- Global triangle numbering
        mutable autoPtr<globalIndex> globalTris_;

        //- The (wanted) distribution type.
        distributionType distType_;

        //- The (current) distribution type. Used to trigger re-distribution
        //  when starting from undecomposed surface.
        distributionType currentDistType_;


    // Private Member Functions

        // Read

            //- Search for io.local directory (=triSurface) either in case
            //  directory or in parent directory
            static word findLocalInstance(const IOobject& io);

            //- Read my additional data
            bool read();


        // Line intersection

            static bool isLocal
            (
                const List<treeBoundBox>& myBbs,
                const point& start,
                const point& end
            );

            //- Split segment into subsegments overlapping the processor
            //  bounding box.
            //void Foam::distributedTriSurfaceMesh::splitSegment
            //(
            //    const label segmentI,
            //    const point& start,
            //    const point& end,
            //    const treeBoundBox& bb,
            //
            //    DynamicList<segment>& allSegments,
            //    DynamicList<label>& allSegmentMap,
            //    DynamicList<label> sendMap
            //) const

            //- Distribute segments into overlapping processor
            //  bounding boxes. Sort per processor.
            void distributeSegment
            (
                const label,
                const point& start,
                const point& end,

                DynamicList<segment>&,
                DynamicList<label>&,
                List<DynamicList<label>>&
            ) const;

            //- Divide edges into local and remote segments. Construct map to
            //  distribute and collect data.
            autoPtr<mapDistribute> distributeSegments
            (
                const pointField& start,
                const pointField& end,

                List<segment>& allSegments,
                List<label>& allSegmentMap
            ) const;

            //- Split edges, distribute, test and collect.
            void findLine
            (
                const bool nearestIntersection,
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>& info
            ) const;


        // Triangle index


            //- Helper: convert local triangle indices to global ones
            void convertTriIndices(List<pointIndexHit>& info) const;

            //- Obtains global indices from pointIndexHit and swaps them back
            //  to their original processor. Used to calculate local region
            //  and normal.
            autoPtr<mapDistribute> calcLocalQueries
            (
                const List<pointIndexHit>&,
                labelList& triangleIndex
            ) const;


        // Nearest

            //- Is location inside any of the bounding boxes
            bool contains
            (
                const List<treeBoundBox>& bbs,
                const point& sample
            ) const;

            //- Is location inside any of the processors bb or if not
            //  does it overlap
            Tuple2<label, scalar> findBestProcs
            (
                const point& centre,
                const scalar radiusSqr,
                boolList& procContains,
                boolList& procOverlaps,
                label& minProci
            ) const;

            label calcOverlappingProcs
            (
                const point& centre,
                const scalar radiusSqr,
                boolList& overlaps
            ) const;

            //- Calculate map to send centres+radius to processors
            autoPtr<mapDistribute> calcLocalQueries
            (
                const bool includeLocalProcessor,
                const pointField& centres,
                const scalarField& radiusSqr,

                pointField& allCentres,
                scalarField& allRadiusSqr,
                labelList& allSegmentMap
            ) const;


        // Side

            //- Side of nearest point w.r.t. edge between face0 and face1
            volumeType edgeSide
            (
                const point& sample,
                const point& nearestPoint,
                const label face0,
                const label face1
            ) const;

            //- Find edge-connected face
            label findOtherFace
            (
                const labelListList& pointFaces,
                const label nearFacei,
                const label nearLabel
            ) const;

            //- Side of nearest point on faces w.r.t. samples. Handles nearest
            //  on edge/point
            void surfaceSide
            (
                const pointField& samples,
                const List<pointIndexHit>& nearestInfo,
                List<volumeType>& region
            ) const;


            // Caching of volume type (based on indexedOctree)

            //- Collect mid points of tree boxes
            void collectLeafMids
            (
                const label nodeI,
                DynamicField<point>& midPoints
            ) const;

            //- Find volume type of tree boxes
            volumeType calcVolumeType
            (
                const List<volumeType>& midPointTypes,
                label& index,
                PackedList<2>& nodeTypes,
                const label nodeI
            ) const;

            //- Look up any cached data. Return unknown if cannot be determined.
            volumeType cachedVolumeType
            (
                const label nodeI,
                const point& sample
            ) const;


        // Surface redistribution

            //- Calculate face-faces
            static void calcFaceFaces
            (
                const triSurface& s,
                const labelListList& pointFaces,
                labelListList& faceFaces
            );

            //- Finds new bounds based on an independent decomposition.
            List<List<treeBoundBox>> independentlyDistributedBbs
            (
                const triSurface&
            );

            //- Does any part of triangle overlap bb.
            static bool overlaps
            (
                const List<treeBoundBox>& bb,
                const triPointRef& tri
            );

            //- Find points used in subset
            static void subsetMeshMap
            (
                const triSurface& s,
                const boolList& include,
                const label nIncluded,
                labelList& newToOldPoints,
                labelList& oldToNewPoints,
                labelList& newToOldFaces
            );

            //- Construct subsetted surface
            static triSurface subsetMesh
            (
                const triSurface& s,
                const labelList& newToOldPoints,
                const labelList& oldToNewPoints,
                const labelList& newToOldFaces
            );

            //- Subset given marked faces
            static triSurface subsetMesh
            (
                const triSurface& s,
                const boolList& include,
                labelList& newToOldPoints,
                labelList& newToOldFaces
            );

            //- Subset given marked faces
            static triSurface subsetMesh
            (
                const triSurface& s,
                const labelList& newToOldFaces,
                labelList& newToOldPoints
            );

            //- Find triangle otherF in allFaces.
            static label findTriangle
            (
                const List<labelledTri>& allFaces,
                const labelListList& allPointFaces,
                const labelledTri& otherF
            );

            //- Merge triSurface (subTris, subPoints) into allTris, allPoints.
            static void merge
            (
                const scalar mergeDist,
                const List<labelledTri>& subTris,
                const pointField& subPoints,

                List<labelledTri>& allTris,
                pointField& allPoints,

                labelList& faceConstructMap,
                labelList& pointConstructMap
            );

            //- Distribute stored fields
            template<class Type>
            void distributeFields(const mapDistribute& map);


        //- No copy construct
        distributedTriSurfaceMesh(const distributedTriSurfaceMesh&) = delete;

        //- No copy assignment
        void operator=(const distributedTriSurfaceMesh&) = delete;


public:

    //- Runtime type information
    TypeName("distributedTriSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        distributedTriSurfaceMesh
        (
            const IOobject&,
            const triSurface&,
            const dictionary& dict
        );

        //- Construct read. Does findInstance to find io.local().
        distributedTriSurfaceMesh(const IOobject& io);

        //- Construct from dictionary (used by searchableSurface).
        //  Does read. Does findInstance to find io.local().
        distributedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributedTriSurfaceMesh();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Triangle indexing (demand driven)
        const globalIndex& globalTris() const;


        // searchableSurface implementation

            //- Range of global indices that can be returned.
            virtual label globalSize() const
            {
                return globalTris().totalSize();
            }

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            //- Find the nearest locations for the supplied points to a
            //  particular region in the searchable surface.
            virtual void findNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                List<pointIndexHit>& info
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;

            //- Set bounds of surface. Bounds currently set as list of
            //  bounding boxes. Will do redistribution of surface to locally
            //  have all triangles overlapping bounds.
            //  Larger bounds: more triangles (memory), more fully local tests
            //  (quick).
            //  keepNonLocal = true : keep triangles that do not overlap
            //  any processor bounds.
            //  Should really be split into a routine to determine decomposition
            //  and one that does actual distribution but determining
            //  decomposition with duplicate triangle merging requires
            //  same amount as work as actual distribution.
            virtual void distribute
            (
                const List<treeBoundBox>&,
                const bool keepNonLocal,
                autoPtr<mapDistribute>& faceMap,
                autoPtr<mapDistribute>& pointMap
            );


        // Other

            //- WIP. From a set of hits (points and
            //  indices) get the specified field. Misses do not get set.
            virtual void getField(const List<pointIndexHit>&, labelList&) const;

            //- Calculate the triangles that are overlapping bounds.
            static void overlappingSurface
            (
                const triSurface&,
                const List<treeBoundBox>&,
                boolList& includedFace
            );

            //- Subset the part of surface that is overlapping bounds.
            static triSurface overlappingSurface
            (
                const triSurface&,
                const List<treeBoundBox>&,
                labelList& subPointMap,
                labelList& subFaceMap
            );


            //- Obtains global indices from pointIndexHit and swaps them back
            //  to their original processor. Used to calculate local region
            //  and normal.
            virtual autoPtr<mapDistribute> localQueries
            (
                const List<pointIndexHit>&,
                labelList& triangleIndex
            ) const;


            //- Print some stats. Parallel aware version of
            //  triSurface::writeStats.
            void writeStats(Ostream& os) const;


        // regIOobject implementation

            //- Write using stream options
            //  Do not use the triSurfaceMesh::writeObject since it
            //  would filter out empty regions. These need to be preserved
            //  in case we want to make decisions based on the number of
            //  regions.
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool writeOnProc
            ) const;

            //- Is object global
            virtual bool global() const
            {
                return false;
            }

            //- Return complete path + object name if the file exists
            //  either in the case/processor or case otherwise null
            virtual fileName filePath() const
            {
                return searchableSurface::localFilePath(type());
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for distributedTriSurfaceMesh
template<>
struct is_globalIOobject<distributedTriSurfaceMesh> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributedTriSurfaceMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
