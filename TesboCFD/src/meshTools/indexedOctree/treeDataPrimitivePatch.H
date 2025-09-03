/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::treeDataPrimitivePatch

Description
    Encapsulation of data needed to search on PrimitivePatches

SourceFiles
    treeDataPrimitivePatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDataPrimitivePatch_H
#define Foam_treeDataPrimitivePatch_H

#include "treeBoundBoxList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;


/*---------------------------------------------------------------------------*\
                 Class treeDataPrimitivePatchName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(treeDataPrimitivePatch);


/*---------------------------------------------------------------------------*\
                   Class treeDataPrimitivePatch Declaration
\*---------------------------------------------------------------------------*/

template<class PatchType>
class treeDataPrimitivePatch
:
    public treeDataPrimitivePatchName
{
    // Private Data

        //- Underlying geometry
        const PatchType& patch_;

        //- Whether to precalculate and store face bounding box
        const bool cacheBb_;

        //- Tolerance to use for intersection tests
        const scalar planarTol_;

        //- Face bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Get face bounding box
        inline treeBoundBox getBounds(const label patchFacei) const;

        //- Initialise all member data
        void update();


public:


    class findNearestOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataPrimitivePatch>& tree);

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

        //- Calculates nearest (to line) point in shape.
        //  Returns point and distance (squared)
        void operator()
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;
    };


    class findIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

    public:

        findIntersectOp(const indexedOctree<treeDataPrimitivePatch>& tree);

        //- Calculate intersection of any face with ray.
        //  Sets result accordingly. Used to find first intersection.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    class findAllIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

        DynamicList<label>& shapeMask_;

    public:

        findAllIntersectOp
        (
            const indexedOctree<treeDataPrimitivePatch>& tree,
            DynamicList<label>& shapeMask
        );

        //- Calculate intersection of unique face with ray.
        //  Sets result accordingly. Used to find all faces.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    class findSelfIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

        const label edgeID_;

    public:

        findSelfIntersectOp
        (
            const indexedOctree<treeDataPrimitivePatch>& tree,
            const label edgeID
        );

        //- Calculate intersection of face with edge of patch.
        //  Excludes faces that use edgeID. Used to find self intersection.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };



    // Constructors (cachable versions)

        //- Construct from patch.
        treeDataPrimitivePatch
        (
            const bool cacheBb,
            const PatchType& patch,
            const scalar planarTol
        );


    // Constructors (non-caching)

        //- Construct from patch.
        treeDataPrimitivePatch(const PatchType& patch, const scalar planarTol)
        :
            treeDataPrimitivePatch(false, patch, planarTol)
        {}


    // Static Functions

        //- Calculate and return bounding boxes for each patch face
        static treeBoundBoxList boxes(const PatchType& patch);


    // Member Functions

        //- Object dimension == 2 (face element)
        int nDim() const noexcept { return 2; }

        //- Return bounding box for the specified face indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- The underlying patch
        const PatchType& patch() const noexcept { return patch_; }

        //TDB //- The subset of face ids to use
        ///const labelList& faceLabels() const noexcept { labelList::null(); }

        //- Use a subset of the patch
        bool useSubset() const noexcept { return false; }

        //- Is the patch empty (no faces)?
        bool empty() const { return !patch_.size(); }

        //- The patch size
        label size() const { return patch_.size(); }

        //- Map from shape index to original (non-subset) face label
        label objectIndex(const label index) const noexcept { return index; }

        //TBD //- Face at specified shape index
        ///const typename PatchType::face_type&
        ///operator[](label index) const { return patch_[index]; }

        //- Representative point (face centre) at shape index
        const point& centre(const label index) const
        {
            return patch_.faceCentres()[index];
        }

        //- Representative point cloud for contained shapes.
        //- One point per shape, corresponding to the face centres.
        tmp<pointField> centres() const
        {
            return patch_.faceCentres();
        }


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const indexedOctree<treeDataPrimitivePatch<PatchType>>&,
            const point&
        ) const;

        //- Does shape at index overlap searchBox
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Does shape at index overlap sphere
        bool overlaps
        (
            const label index,
            const point& centre,
            const scalar radiusSqr
        ) const;

        //- Helper: find intersection of line with shapes
        static bool findIntersection
        (
            const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        );

        //- Calculates nearest (to sample) point in shape.
        //  Returns actual point and distance (squared)
        void findNearest
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& nearestIndex,
            point& nearestPoint
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "treeDataPrimitivePatch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
