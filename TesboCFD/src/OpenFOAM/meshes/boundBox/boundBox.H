/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::boundBox

Description
    A bounding box defined in terms of min/max extrema points.

Note
    When a bounding box is created without any points, it creates an inverted
    bounding box. Points can be added later and the bounding box will grow to
    include them.

SeeAlso
    Foam::treeBoundBox

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundBox_H
#define Foam_boundBox_H

#include "pointField.H"
#include "faceList.H"
#include "Pair.H"
#include "triangleFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBox;
class plane;
class Random;
template<class T> class MinMax;

Istream& operator>>(Istream& is, boundBox& bb);
Ostream& operator<<(Ostream& os, const boundBox& bb);


/*---------------------------------------------------------------------------*\
                          Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private Data

        //- Minimum and maximum points describing the bounding box
        point min_, max_;


protected:

    // Protected Member Functions

        //- Test for overlap of box and box (inclusive check)
        inline static bool box_box_overlaps
        (
            const point& minA,  // boxA (min)
            const point& maxA,  // boxA (max)
            const point& minB,  // boxB (min)
            const point& maxB   // boxB (max)
        );

        //- Test for overlap of box and boundingSphere (centre + sqr(radius))
        //  Note: ordering of corners is irrelevant
        inline static bool box_sphere_overlaps
        (
            const point& corner0,   // box corner
            const point& corner1,   // box corner
            const point& centre,    // sphere centre
            const scalar radiusSqr  // sqr(radius)
        );


public:

    // Static Data Members

        //- Bits used for (x/y/z) directional encoding.
        enum directionBit : direction
        {
            XDIR = 1,  //!< 1: x-direction. Same as (1 << vector::X)
            YDIR = 2,  //!< 2: y-direction. Same as (1 << vector::Y)
            ZDIR = 4   //!< 4: z-direction. Same as (1 << vector::Z)
        };

        //- A large boundBox: min/max == -/+ ROOTVGREAT
        static const boundBox greatBox;

        //- A large inverted boundBox: min/max == +/- ROOTVGREAT
        static const boundBox invertedBox;

        //- The unit normal per face
        static const FixedList<vector, 6> faceNormals;


    // Static Methods

        //- The null boundBox is the same as an inverted box
        static const boundBox& null() noexcept
        {
            return invertedBox;
        }

        //- The boundBox faces as a hexCell, using hexCorner points.
        //- Same as hexCell::modelFaces()
        static const Foam::faceList& hexFaces();

        //- Number of points for boundBox and HEX
        static constexpr label nPoints() noexcept { return 8; }

        //- Number of edges for boundBox and HEX
        static constexpr label nEdges() noexcept { return 12; }

        //- Number of faces for boundBox and HEX
        static constexpr label nFaces() noexcept { return 6; }


    // Standard (Generated) Methods

        //- Default construct: an inverted bounding box
        inline boundBox();

        //- Copy construct
        boundBox(const boundBox&) = default;

        //- Copy assignment
        boundBox& operator=(const boundBox&) = default;


    // Constructors

        //- Copy construct with specified global reduction
        boundBox(const boundBox& bb, bool doReduce);

        //- Construct a bounding box containing a single initial point
        inline explicit boundBox(const point& p);

        //- Construct a 0/1 unit bounding box
        inline explicit boundBox(const Foam::zero_one);

        //- Construct from bound box min/max points
        inline boundBox(const point& min, const point& max);

        ///TBD: Construct from bound box min/max points
        ///inline explicit boundBox(const MinMax<point>& bb);

        //- Construct from bound box min/max points
        inline explicit boundBox(const Pair<point>& bb);

        //- Construct as the bounding box of the given points
        //  Does parallel communication (doReduce = true)
        explicit boundBox(const UList<point>& points, bool doReduce = true);

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (doReduce = true)
        explicit boundBox(const tmp<pointField>& tpoints, bool doReduce = true);

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from cell/face etc.
        //  Does parallel communication (doReduce = true)
        boundBox
        (
            const UList<point>& points,
            const labelUList& indices,
            bool doReduce = true
        );

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from edge/triFace etc.
        //  Does parallel communication (doReduce = true)
        template<unsigned N>
        boundBox
        (
            const UList<point>& points,
            const FixedList<label, N>& indices,
            bool doReduce = true
        );

        //- Construct from Istream
        inline explicit boundBox(Istream& is);


    // Member Functions

        // Access

            //- Bounding box is inverted, contains no points.
            inline bool empty() const;

            //- Bounding box is non-inverted.
            inline bool good() const;

            //- Bounding box is non-inverted - same as good().
            bool valid() const { return good(); }

            //- Minimum describing the bounding box
            inline const point& min() const noexcept;

            //- Maximum describing the bounding box
            inline const point& max() const noexcept;

            //- Minimum describing the bounding box, non-const access
            inline point& min() noexcept;

            //- Maximum describing the bounding box, non-const access
            inline point& max() noexcept;

            //- The centre (midpoint) of the bounding box
            inline point centre() const;

            //- The bounding box span (from minimum to maximum)
            inline vector span() const;

            //- The magnitude/length of the bounding box diagonal
            inline scalar mag() const;

            //- The magnitude/length squared of bounding box diagonal
            inline scalar magSqr() const;

            //- The volume of the bound box
            inline scalar volume() const;

            //- Smallest length/height/width dimension
            inline scalar minDim() const;

            //- Largest length/height/width dimension
            inline scalar maxDim() const;

            //- Average length/height/width dimension
            inline scalar avgDim() const;

            //- Direction (X/Y/Z) of the smallest span (for empty box: 0)
            inline direction minDir() const;

            //- Direction (X/Y/Z) of the largest span (for empty box: 0)
            inline direction maxDir() const;

            //- Count the number of positive, non-zero dimensions.
            //  \return
            //   - -1 : if any dimensions are negative
            //   -  0 : 0D (point)
            //   -  1 : 1D (line aligned with an axis)
            //   -  2 : 2D (plane aligned with an axis)
            //   -  3 : 3D (box)
            inline int nDim() const;

            //- Return corner point [0..7] corresponding to a 'hex' cell
            template<direction CornerNumber>
            inline point hexCorner() const;

            //- Corner points in an order corresponding to a 'hex' cell
            tmp<pointField> hexCorners() const;

            //- Corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const { return hexCorners(); }

            //- Face midpoints
            tmp<pointField> faceCentres() const;

            //- Face centre of given face index
            point faceCentre(const direction facei) const;


        // Manipulate

            //- Reset to an inverted box
            inline void reset();

            //- Reset to a 0/1 unit bounding box
            inline void reset(const Foam::zero_one);

            //- Reset min/max to be identical to the specified point
            inline void reset(const point& pt);

            //- Reset min/max to specified values
            inline void reset(const point& min, const point& max);

            //- Same as reset() - reset to an inverted box
            void clear() { reset(); }

            //- Extend to include the second box.
            inline void add(const boundBox& bb);

            //- Extend to include the point.
            inline void add(const point& pt);

            //- Extend to include two additional points.
            inline void add(const point& pt0, const point& pt1);

            //- Extend to include two additional points.
            inline void add(const Pair<point>& points);

            //- Extend to include the points.
            inline void add(const UList<point>& points);

            //- Extend to include the points from the temporary point field.
            inline void add(const tmp<pointField>& tpoints);

            //- Extend to include the points.
            template<unsigned N>
            void add(const FixedList<point, N>& points);

            //- Extend to include a (subsetted) point field.
            //  The indices could be from edge/triFace etc.
            template<unsigned N>
            void add
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            );

            //- Extend to include a (subsetted) point field.
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            void add
            (
                const UList<point>& points,
                const IntContainer& indices
            );

            //- Expand box by adjusting min/max by specified amount
            //- in each dimension
            inline void grow(const scalar delta);

            //- Expand box by adjusting min/max by specified amounts
            inline void grow(const vector& delta);

            //- Expand box by factor*mag(span) in all dimensions
            inline void inflate(const scalar factor);

            //- Expand box slightly by expanding all dimensions with
            //- factor*span*(random 0-1) and guarantees
            //- factor*mag(span) minimum width in any direction
            void inflate(Random& rndGen, const scalar factor);

            //- As per two parameter version but with additional
            //- grow() by given amount in each dimension
            void inflate(Random& r, const scalar factor, const scalar delta);


        // Parallel

            //- Inplace parallel reduction of min/max values
            void reduce();

            //- Perform a reduction on a copy and return the result
            static boundBox returnReduce(const boundBox& bb);


        // Query

            //- Does plane intersect this bounding box.
            //  There is an intersection if the plane segments the corner points
            //  \note the results are unreliable when plane coincides almost
            //  exactly with a box face
            bool intersects(const plane& pln) const;

            //- Does triangle intersect this bounding box
            //- or is contained within this bounding box.
            //  \note results may be unreliable when it coincides almost
            //  exactly with a box face
            bool intersects(const triPointRef& tri) const;

            //- Overlaps/touches boundingBox?
            inline bool overlaps(const boundBox& bb) const;

            //- Overlaps boundingSphere (centre + sqr(radius))?
            inline bool overlaps
            (
                const point& centre,
                const scalar radiusSqr
            ) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point& pt) const;

            //- Fully contains other boundingBox?
            inline bool contains(const boundBox& bb) const;

            //- Contains point? (inside only)
            inline bool containsInside(const point& pt) const;

            //- Contains all points? (inside or on edge)
            bool contains(const UList<point>& points) const;

            //- Contains all of the (subsetted) points? (inside or on edge)
            template<unsigned N>
            bool contains
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            ) const;


            //- Contains all of the (subsetted) points? (inside or on edge)
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            bool contains
            (
                const UList<point>& points,
                const IntContainer& indices
            ) const;


            //- Contains any of the points? (inside or on edge)
            bool containsAny(const UList<point>& points) const;

            //- Contains any of the (subsetted) points? (inside or on edge)
            template<unsigned N>
            bool containsAny
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            ) const;

            //- Contains any of the (subsetted) points? (inside or on edge)
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            bool containsAny
            (
                const UList<point>& points,
                const IntContainer& indices
            ) const;


            //- Return the nearest point on the boundBox to the supplied point.
            //  If point is inside the boundBox then the point is returned
            //  unchanged.
            point nearest(const point& p) const;


    // Member Operators

        //- Restrict min/max to union with other box.
        void operator&=(const boundBox& bb);

        //- Extend box to include the second box, as per the add() method.
        //  Can be used in a reduction operation.
        void operator+=(const boundBox& bb) { add(bb); }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, boundBox& bb);
        friend Ostream& operator<<(Ostream& os, const boundBox& bb);


    // Housekeeping

        //- Deprecated(2022-10) - use 'operator&=' to avoid confusion with
        //- other intersects() methods.
        //  \deprecated(2022-10) - use 'operator&='
        FOAM_DEPRECATED_FOR(2022-10, "boundBox::operator&=(const boundBox&)")
        bool intersect(const boundBox& bb) { *this &= bb; return good(); }

        //- Identical to centre()
        point midpoint() const { return centre(); }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for boundBox
template<> struct is_contiguous<boundBox> : is_contiguous<point> {};

//- Contiguous scalar data for boundBox
template<> struct is_contiguous_scalar<boundBox>
:
    is_contiguous_scalar<point>
{};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const boundBox& a, const boundBox& b)
{
    return (a.min() == b.min() && a.max() == b.max());
}


inline bool operator!=(const boundBox& a, const boundBox& b)
{
    return !(a == b);
}


inline bool operator<(const boundBox& a, const boundBox& b)
{
    return
    (
        a.min() < b.min()
     || (!(b.min() < a.min()) && a.max() < b.max())
    );
}


inline bool operator<=(const boundBox& a, const boundBox& b)
{
    return !(b < a);
}


inline bool operator>(const boundBox& a, const boundBox& b)
{
    return (b < a);
}


inline bool operator>=(const boundBox& a, const boundBox& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "boundBoxI.H"

#ifdef NoRepository
    #include "boundBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
