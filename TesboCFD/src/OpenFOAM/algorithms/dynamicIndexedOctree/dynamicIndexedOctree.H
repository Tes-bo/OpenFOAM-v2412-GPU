/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::dynamicIndexedOctree

Description
    Non-pointer based hierarchical recursive searching.
    Storage is dynamic, so elements can be deleted.

SourceFiles
    dynamicIndexedOctree.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicIndexedOctree_H
#define Foam_dynamicIndexedOctree_H

#include "indexedOctree.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class dynamicIndexedOctree;
template<class Type>
Ostream& operator<<(Ostream&, const dynamicIndexedOctree<Type>&);


/*---------------------------------------------------------------------------*\
                  Class dynamicIndexedOctreeBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for dynamicIndexedOctree
//  Same type of node bookkeeping as indexedOctree
class dynamicIndexedOctreeBase
:
    public indexedOctreeBase
{
public:

    //- Document that we are using the same types of node
    using node = indexedOctreeBase::node;

    //- Runtime type information
    ClassName("dynamicIndexedOctree");


    // Constructors

        //- Default construct
        dynamicIndexedOctreeBase() = default;
};


/*---------------------------------------------------------------------------*\
                           Class dynamicIndexedOctree Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dynamicIndexedOctree
:
    public dynamicIndexedOctreeBase
{
    // Private Data

        //- Underlying shapes for geometric queries.
        const Type shapes_;

        const treeBoundBox bb_;

        const label maxLevels_;

        //- Current number of levels in the tree.
        label nLevelsMax_;

        const scalar maxLeafRatio_;

        const label minSize_;

        const scalar maxDuplicity_;

        //- List of all nodes
        DynamicList<node> nodes_;

        //- List of all contents (referenced by those nodes that are contents)
        DynamicList<DynamicList<label>> contents_;

        //- Per node per octant whether is fully inside/outside/mixed.
        mutable PackedList<2> nodeTypes_;


    // Private Member Functions

        // Construction

            //- Split list of indices into 8 bins
            //  according to where they are in relation to mid.
            void divide
            (
                const labelUList& indices,
                const treeBoundBox& bb,
                FixedList<DynamicList<label>, 8>& dividedIndices
            ) const;

            //- Subdivide the contents node at position contentI.
            //  Appends to contents.
            node divide
            (
                const treeBoundBox& bb,
                label contentIndex,
                const label parentNodeIndex,
                const label octantToBeDivided
            );

            // Recursively call the divide function
            void recursiveSubDivision
            (
                const treeBoundBox& subBb,
                const label contentI,
                const label parentIndex,
                const label octant,
                label& nLevels
            );

            //- Determine inside/outside per node (mixed if cannot be
            //  determined). Only valid for closed shapes.
            volumeType calcVolumeType(const label nodeI) const;

            //- Search cached volume type.
            volumeType getVolumeType(const label nodeI, const point&) const;


        // Query

            //- Find nearest point to line.
            void findNearest
            (
                const label nodeI,
                const linePointRef& ln,

                treeBoundBox& tightest,
                label& nearestShapeI,
                point& linePoint,
                point& nearestPoint
            ) const;

            //- Return bbox of octant
            treeBoundBox subBbox
            (
                const label parentNodeI,
                const direction octant
            ) const;

            //- Helper: take a point on/close to face of bb and push it
            //  inside or outside of bb.
            static point pushPoint
            (
                const treeBoundBox&,
                const point&,
                const bool pushInside
            );

            //- Helper: take a point on face of bb and push it
            //  inside or outside of bb.
            static point pushPoint
            (
                const treeBoundBox&,
                const direction,
                const point&,
                const bool pushInside
            );

            //- Helper: take point on face(s) of bb and push it away from
            //  edges of face.  If pt is not on a face does nothing.
            static point pushPointIntoFace
            (
                const treeBoundBox& bb,
                const vector& dir,          // end-start
                const point& pt
            );

            //- Walk to parent of node+octant.
            //  Return false if top of tree reached.
            bool walkToParent
            (
                const label nodeI,
                const direction octant,

                label& parentNodeI,
                label& parentOctant
            ) const;

            //- Walk tree to neighbouring node. Return false if edge of tree
            //  hit.
            bool walkToNeighbour
            (
                const point& facePoint,
                const direction faceID,         // direction to walk in
                label& nodeI,
                direction& octant
            ) const;

            //- Debug: return verbose the bounding box faces
            static word faceString(const direction faceID);

            //- Traverse a node. If intersects a triangle return first
            // intersection point.
            // findAny=true : return any intersection
            // findAny=false: return nearest (to start) intersection
            void traverseNode
            (
                const bool findAny,
                const point& treeStart,
                const vector& treeVec,

                const point& start,
                const point& end,
                const label nodeI,
                const direction octantI,

                pointIndexHit& hitInfo,
                direction& faceID
            ) const;

            //- Find any or nearest intersection
            pointIndexHit findLine
            (
                const bool findAny,
                const point& treeStart,
                const point& treeEnd,
                const label startNodeI,
                const direction startOctantI,
                const bool verbose = false
            ) const;

            //- Find any or nearest intersection of line between start and end.
            pointIndexHit findLine
            (
                const bool findAny,
                const point& start,
                const point& end
            ) const;

            //- Find elements intersecting box
            //  Store all results in elements (if non-null), or early exit
            bool findBox
            (
                const label nodeI,
                const treeBoundBox& searchBox,
                labelHashSet* elements
            ) const;

            //- Find elements intersecting sphere.
            //  Store all results in elements (if non-null), or early exit
            bool findSphere
            (
                const label nodeI,
                const point& centre,
                const scalar radiusSqr,
                labelHashSet* elements
            ) const;


            template<class CompareOp>
            static void findNear
            (
                const scalar nearDist,
                const bool okOrder,
                const dynamicIndexedOctree<Type>& tree1,
                const labelBits index1,
                const treeBoundBox& bb1,
                const dynamicIndexedOctree<Type>& tree2,
                const labelBits index2,
                const treeBoundBox& bb2,
                CompareOp& cop
            );


    // Other

        //- Count number of elements on this and sublevels
        label countElements(const labelBits index) const;

        //- Write node treeBoundBoxes in OBJ format
        void writeOBJ
        (
            const label nodeI,
            Ostream& os,
            label& vertIndex,
            const bool leavesOnly,
            const bool writeLinesOnly = false
        ) const;

        //- Dump node+octant to an obj file
        void writeOBJ(const label nodeI, const direction octant) const;


public:

    // Constructors

        //- Construct from shapes
        dynamicIndexedOctree
        (
            const Type& shapes,
            const treeBoundBox& bb,
            const label maxLevels,          // maximum number of levels
            const scalar maxLeafRatio,      // how many elements per leaf
            const scalar maxDuplicity       // in how many leaves is a shape on
                                            // average
        );

        //- Clone
        autoPtr<dynamicIndexedOctree<Type>> clone() const
        {
            return autoPtr<dynamicIndexedOctree<Type>>::New(*this);
        }


    // Member Functions

    // Access

        //- Reference to shape
        const Type& shapes() const noexcept { return shapes_; }

        //- List of all nodes
        const List<node>& nodes() const noexcept { return nodes_; }

        //- List of all contents
        //- (referenced by those nodes that are contents)
        const DynamicList<DynamicList<label>>& contents() const noexcept
        {
            return contents_;
        }

        //- Per node, per octant whether is fully inside/outside/mixed.
        PackedList<2>& nodeTypes() const noexcept
        {
            return nodeTypes_;
        }

        //- Top bounding box
        const treeBoundBox& bb() const
        {
            if (nodes_.empty())
            {
                return treeBoundBox::null();
            }
            return nodes_[0].bb_;
        }


        // Queries

            //- Calculate nearest point on nearest shape.
            //  Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: index in shapes
            //  - point: actual nearest point found
            pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Low level: calculate nearest starting from subnode.
            void findNearest
            (
                const label nodeI,
                const point&,

                scalar& nearestDistSqr,
                label& nearestShapeI,
                point& nearestPoint
            ) const;

            //- Find nearest to line.
            //  Returns
            //  - bool : any point found?
            //  - label: index in shapes
            //  - point: actual nearest point found
            //  sets:
            //  - linePoint : corresponding nearest point on line
            pointIndexHit findNearest
            (
                const linePointRef& ln,
                treeBoundBox& tightest,
                point& linePoint
            ) const;

            //- Find nearest intersection of line between start and end.
            pointIndexHit findLine
            (
                const point& start,
                const point& end
            ) const;

            //- Find any intersection of line between start and end.
            pointIndexHit findLineAny
            (
                const point& start,
                const point& end
            ) const;

            //- True if any shapes overlap the bounding box
            bool overlaps(const treeBoundBox& bb) const;

            //- Find indices of all shapes inside or overlapping
            //- a bounding box (i.e. all shapes not outside box)
            //  \returns the indices (in no particular order)
            labelList findBox
            (
                const treeBoundBox& bb  //!< bound box limits
            ) const;

            //- Find indices of all shapes inside or overlapping
            //- a bounding box (i.e. all shapes not outside box)
            //  \returns the number of elements found
            label findBox
            (
                const treeBoundBox& bb, //!< bound box limits
                labelHashSet& elements  //!< [out] elements found
            ) const;

            //- True if any shapes overlap the bounding sphere
            bool overlaps
            (
                const point& centre,    //!< centre of bound sphere
                const scalar radiusSqr  //!< radius^2 of sphere
            ) const;

            //- Find indices of all shapes inside or overlapping
            //- a bounding sphere (i.e. all shapes not outside a sphere)
            //  \returns the indices (in no particular order)
            labelList findSphere
            (
                const point& centre,    //!< centre of bound sphere
                const scalar radiusSqr  //!< radius^2 of sphere
            ) const;

            //- Find indices of all shapes inside or overlapping
            //- a bounding sphere (i.e. all shapes not outside sphere)
            //  \returns the number of elements found
            label findSphere
            (
                const point& centre,    //!< centre of bound sphere
                const scalar radiusSqr, //!< radius^2 of sphere
                labelHashSet& elements  //!< [out] elements found
            ) const;


            //- Find deepest node (as parent+octant) containing point. Starts
            //  off from starting index in nodes_ (use 0 to start from top)
            //  Use getNode and getOctant to extract info, or call findIndices.
            labelBits findNode(const label nodeI, const point&) const;

            //- Find shape containing point. Only implemented for certain
            //  shapes.
            label findInside(const point&) const;

            //- Find the shape indices that occupy the result of findNode
            const labelList& findIndices(const point&) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            volumeType getVolumeType(const point&) const;

            //- Helper function to return the side. Returns outside if
            //  outsideNormal&vec >= 0, inside otherwise
            static volumeType getSide
            (
                const vector& outsideNormal,
                const vector& vec
            );

            //- Find near pairs and apply CompareOp to them.
            //  tree2 can be *this or different tree.
            template<class CompareOp>
            void findNear
            (
                const scalar nearDist,
                const dynamicIndexedOctree<Type>& tree2,
                CompareOp& cop
            ) const;


        // Edit

            //- Insert a new object into the tree.
            bool insert(label startIndex, label endIndex);

            bool insertIndex
            (
                const label nodIndex,
                const label index,
                label& nLevels
            );

            //- Remove an object from the tree.
            bool remove(const label index);

            label removeIndex(const label nodIndex, const label index);


    // Write

        //- Write all tree boxes as OBJ format
        void writeOBJ(Ostream& os) const;

        //- Print tree. Either print all indices (printContent = true) or
        //  just size of contents nodes.
        void print
        (
            prefixOSstream&,
            const bool printContents,
            const label
        ) const;

        bool write(Ostream& os) const;

        void writeTreeInfo() const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const dynamicIndexedOctree<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicIndexedOctree.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
