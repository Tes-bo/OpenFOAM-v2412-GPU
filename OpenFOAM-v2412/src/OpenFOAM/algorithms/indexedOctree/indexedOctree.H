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
    Foam::indexedOctree

Description
    Non-pointer based hierarchical recursive searching

SourceFiles
    indexedOctree.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_indexedOctree_H
#define Foam_indexedOctree_H

#include "treeBoundBox.H"
#include "pointIndexHit.H"
#include "FixedList.H"
#include "Ostream.H"
#include "HashSet.H"
#include "labelBits.H"
#include "PackedList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;
template<class Type> Ostream& operator<<(Ostream&, const indexedOctree<Type>&);


/*---------------------------------------------------------------------------*\
                      Class indexedOctreeBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for indexedOctree
//  The internal node bookkeeping is encoded follows:
//  - 0: empty
//  - +ve: parent
//  - -ve: leaf
class indexedOctreeBase
{
public:

    // Public Data Types

        //- Tree node. Has up pointer and down pointers.
        class node
        {
        public:

            //- Has exactly 8 sub-nodes (octants)
            static constexpr direction nChildren = 8;

            //- Bounding box of this node
            treeBoundBox bb_;

            //- Parent node (index into flat list addressing for tree)
            label parent_ = -1;

            //- IDs of the 8 nodes on all sides of the mid point
            FixedList<labelBits, 8> subNodes_;


        // Operators

            friend bool operator==(const node& a, const node& b)
            {
                return
                (
                    a.parent_ == b.parent_
                 && a.subNodes_ == b.subNodes_
                 && a.bb_ == b.bb_
                );
            }

            friend bool operator!=(const node& a, const node& b)
            {
                return !(a == b);
            }

            friend Ostream& operator<< (Ostream& os, const node& n)
            {
                return os << n.bb_ << token::SPACE
                    << n.parent_ << token::SPACE << n.subNodes_;
            }

            friend Istream& operator>> (Istream& is, node& n)
            {
                return is >> n.bb_ >> n.parent_ >> n.subNodes_;
            }
        };


    // Static Functions

    // Tests for node handling
    // (0: empty, +ve: parent, -ve: leaf)

        //- An empty node - no content
        static bool isEmpty(labelBits i) noexcept
        {
            return (i.val() == 0);
        }

        //- Node with content (leaf)
        static bool isContent(labelBits i) noexcept
        {
            return (i.val() < 0);
        }

        //- A parent node
        static bool isNode(labelBits i) noexcept
        {
            return (i.val() > 0);
        }


    // Decode/retrieve node addressing

        //- Return real (dereferenced) index for a content node
        static label getContent(labelBits i)
        {
            if (!isContent(i))
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }
            return (-i.val()-1);
        }

        //- Return real (dereferenced) index for a parent node
        static label getNode(const labelBits i)
        {
            if (!isNode(i))
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }
            return (i.val()-1);
        }

        //- Return sub-node direction/octant
        static direction getOctant(labelBits i) noexcept
        {
            return i.bits();
        }


protected:

    // Static Data

        //- Relative perturbation tolerance.
        //  Determines when point is considered to be close to face/edge
        //  of bb of node.
        //  The tolerance is relative to the bounding box of the smallest
        //  node.
        static scalar perturbTol_;


    // Protected Member Functions

    // Encode node addressing
    // (only used when building)

        //- From empty to subNodes_ entry
        static labelBits emptyPlusOctant(direction octant)
        {
            return labelBits(0, octant);
        }

        //- From index into contents_ to subNodes_ entry
        static labelBits contentPlusOctant(label i, direction octant)
        {
            return labelBits(-i-1, octant);
        }

        //- From index into nodes_ to subNodes_ entry
        static labelBits nodePlusOctant(label i, direction octant)
        {
            return labelBits(i+1, octant);
        }


public:

    //- Get the perturbation tolerance
    static scalar& perturbTol() noexcept { return perturbTol_; }

    //- Set the perturbation tolerance, return the old value
    static scalar perturbTol(scalar tol) noexcept
    {
        scalar old(perturbTol_);
        perturbTol_ = tol;
        return old;
    }


    //- Runtime type information
    ClassName("indexedOctree");


    // Constructors

        //- Default construct
        indexedOctreeBase() = default;


    // Output Helpers

        //- Write treeBoundBox in OBJ format
        static void writeOBJ
        (
            Ostream& os,
            const treeBoundBox& bb,
            label& vertIndex,
            const bool writeLinesOnly = false
        );
};


/*---------------------------------------------------------------------------*\
                           Class indexedOctree Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class indexedOctree
:
    public indexedOctreeBase
{
    // Private Data

        //- Underlying shapes for geometric queries.
        const Type shapes_;

        //- List of all nodes
        List<node> nodes_;

        //- List of all contents (referenced by those nodes that are contents)
        List<labelList> contents_;

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
                FixedList<labelList, 8>& dividedIndices
            ) const;

            //- Subdivide the contents node at position contentI.
            //  Appends to contents.
            node divide
            (
                const treeBoundBox& bb,
                DynamicList<labelList>& contents,
                label contentIndex
            ) const;

            //- Split any contents node with more than minSize elements.
            void splitNodes
            (
                const label minSize,
                DynamicList<node>& nodes,
                DynamicList<labelList>& contents
            ) const;

            //- Reorder contents to be in same order as nodes.
            //  Returns number of nodes on the compactLevel.
            static label compactContents
            (
                DynamicList<node>& nodes,
                DynamicList<labelList>& contents,
                const label compactLevel,
                const label nodeI,
                const label level,
                List<labelList>& compactedContents,
                label& compactI
            );

            //- Determine inside/outside per node (mixed if cannot be
            //  determined). Only valid for closed shapes.
            volumeType calcVolumeType(const label nodeI) const;

            //- Search cached volume type.
            volumeType getVolumeType(const label nodeI, const point&) const;


        // Query

            //- Find nearest point to line.
            template<class FindNearestOp>
            void findNearest
            (
                const label nodeI,
                const linePointRef& ln,

                treeBoundBox& tightest,
                label& nearestShapeI,
                point& linePoint,
                point& nearestPoint,

                const FindNearestOp& fnOp
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
            //  edges of face.
            //  Guarantees that if pt is on a face it gets perturbed
            //  so it is away from the face edges.
            //  If pt is not on a face does nothing.
            static point pushPointIntoFace
            (
                const treeBoundBox& bb,
                const vector& dir,          // end-start
                const point& pt
            );

            //- Walk to parent of node+octant.
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
            template<class FindIntersectOp>
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
                direction& faceID,

                const FindIntersectOp& fiOp
            ) const;

            //- Find any or nearest intersection
            template<class FindIntersectOp>
            pointIndexHit findLine
            (
                const bool findAny,
                const point& treeStart,
                const point& treeEnd,
                const label startNodeI,
                const direction startOctantI,
                const FindIntersectOp& fiOp,
                const bool verbose = false
            ) const;

            //- Find any or nearest intersection of line between start and end.
            template<class FindIntersectOp>
            pointIndexHit findLine
            (
                const bool findAny,
                const point& start,
                const point& end,
                const FindIntersectOp& fiOp
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
                const indexedOctree<Type>& tree1,
                const labelBits index1,
                const treeBoundBox& bb1,
                const indexedOctree<Type>& tree2,
                const labelBits index2,
                const treeBoundBox& bb2,
                CompareOp& cop
            );


    // Other

        //- Count number of elements on this and sublevels
        label countElements(const labelBits index) const;

        //- Number of leafs below given node
        label countLeafs(const label nodeI) const;

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

        //- Construct null
        indexedOctree(const Type& shapes);

        //- Construct from components
        indexedOctree
        (
            const Type& shapes,
            const List<node>& nodes,
            const List<labelList>& contents
        );

        //- Construct from shapes
        indexedOctree
        (
            const Type& shapes,
            const treeBoundBox& bb,
            const label maxLevels,          // maximum number of levels
            const scalar maxLeafRatio,      // how many elements per leaf
            const scalar maxDuplicity       // in how many leaves is a shape on
                                            // average
        );

        //- Construct from Istream
        indexedOctree(const Type& shapes, Istream& is);

        //- Clone
        autoPtr<indexedOctree<Type>> clone() const
        {
            return autoPtr<indexedOctree<Type>>::New(*this);
        }


    // Member Functions

    // Access

        //- Reference to shape
        const Type& shapes() const noexcept { return shapes_; }

        //- List of all nodes
        const List<node>& nodes() const noexcept { return nodes_; }

        //- List of all contents
        //- (referenced by those nodes that are contents)
        const List<labelList>& contents() const noexcept
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

        //- Return the number of leaf nodes
        label nLeafs() const;


        // Queries

            pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Calculate nearest point on nearest shape.
            //  Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: index in shapes
            //  - point: actual nearest point found
            template<class FindNearestOp>
            pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr,

                const FindNearestOp& fnOp
            ) const;

            //- Low level: calculate nearest starting from subnode.
            template<class FindNearestOp>
            void findNearest
            (
                const label nodeI,
                const point&,

                scalar& nearestDistSqr,
                label& nearestShapeI,
                point& nearestPoint,

                const FindNearestOp& fnOp
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

            template<class FindNearestOp>
            pointIndexHit findNearest
            (
                const linePointRef& ln,
                treeBoundBox& tightest,
                point& linePoint,

                const FindNearestOp& fnOp
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

            //- Find nearest intersection of line between start and end.
            template<class FindIntersectOp>
            pointIndexHit findLine
            (
                const point& start,
                const point& end,
                const FindIntersectOp& fiOp
            ) const;

            //- Find any intersection of line between start and end.
            template<class FindIntersectOp>
            pointIndexHit findLineAny
            (
                const point& start,
                const point& end,
                const FindIntersectOp& fiOp
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
                const indexedOctree<Type>& tree2,
                CompareOp& cop
            ) const;


    // Write

        //- Write (non-empty) tree boxes in OBJ format
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


    // IOstream Operators

        friend Ostream& operator<< <Type>(Ostream&, const indexedOctree<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "indexedOctree.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
