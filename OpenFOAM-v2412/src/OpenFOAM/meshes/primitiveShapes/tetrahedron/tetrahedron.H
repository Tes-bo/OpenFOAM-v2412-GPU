/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::tetrahedron

Description
    A tetrahedron primitive.

    Ordering of edges needs to be the same for a tetrahedron
    class, a tetrahedron cell shape model and a tetCell.

SourceFiles
    tetrahedronI.H
    tetrahedron.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tetrahedron_H
#define Foam_tetrahedron_H

#include "point.H"
#include "pointHit.H"
#include "primitiveFieldsFwd.H"
#include "Random.H"
#include "FixedList.H"
#include "UList.H"
#include "triangle.H"
#include "treeBoundBox.H"
#include "barycentric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class plane;

template<class Point, class PointRef> class tetrahedron;

template<class Point, class PointRef>
inline Istream& operator>>(Istream&, tetrahedron<Point, PointRef>&);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream&, const tetrahedron<Point, PointRef>&);


// Common Typedefs

//- A tetrahedron using referred points
typedef tetrahedron<point, const point&> tetPointRef;


/*---------------------------------------------------------------------------*\
                         Class tetPoints Declaration
\*---------------------------------------------------------------------------*/

//- Tet point storage. Default constructable (tetrahedron is not)
class tetPoints
:
    public FixedList<point, 4>
{
public:

    // Generated Methods

        //- Default construct
        tetPoints() = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Construct from four points
        inline tetPoints
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );

        //- Construct from point references
        inline explicit tetPoints(const tetPointRef& pts);

        //- Construct from four points
        inline tetPoints(const FixedList<point, 4>& pts);

        //- Copy construct from subset of points
        inline tetPoints
        (
            const UList<point>& points,
            const FixedList<label, 4>& indices
        );


    // Member Functions

        //- The first vertex
        const point& a() const noexcept { return get<0>(); }

        //- The second vertex
        const point& b() const noexcept { return get<1>(); }

        //- The third vertex
        const point& c() const noexcept { return get<2>(); }

        //- The fourth vertex
        const point& d() const noexcept { return get<3>(); }

        //- The first vertex
        point& a() noexcept { return get<0>(); }

        //- The second vertex
        point& b() noexcept { return get<1>(); }

        //- The third vertex
        point& c() noexcept { return get<2>(); }

        //- The fourth vertex
        point& d() noexcept { return get<3>(); }

        //- Return as tetrahedron reference
        inline tetPointRef tet() const;

        //- Invert tetrahedron by swapping third and fourth vertices
        inline void flip();

        //- The enclosing (bounding) box for the tetrahedron
        inline Pair<point> box() const;

        //- The bounding box for the tetrahedron
        inline treeBoundBox bounds() const;
};


/*---------------------------------------------------------------------------*\
                         Class tetrahedron Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class tetrahedron
{
public:

    // Public Typedefs

        //- The point type
        typedef Point point_type;

        //- Storage type for tets originating from intersecting tets.
        //  (can possibly be smaller than 200)
        typedef FixedList<tetPoints, 200> tetIntersectionList;


        // Classes for use in sliceWithPlane. What to do with decomposition
        // of tet.

            //- Dummy
            class dummyOp
            {
            public:
                inline void operator()(const tetPoints&);
            };

            //- Sum resulting volumes
            class sumVolOp
            {
            public:
                scalar vol_;

                inline sumVolOp();

                inline void operator()(const tetPoints&);
            };

            //- Store resulting tets
            class storeOp
            {
                tetIntersectionList& tets_;
                label& nTets_;

            public:
                inline storeOp(tetIntersectionList&, label&);

                inline void operator()(const tetPoints&);
            };

private:

    // Private Data

        PointRef a_, b_, c_, d_;


    // Private Member Functions

        inline static point planeIntersection
        (
            const FixedList<scalar, 4>&,
            const tetPoints&,
            const label,
            const label
        );

        template<class TetOp>
        inline static void decomposePrism
        (
            const FixedList<point, 6>& points,
            TetOp& op
        );

        template<class AboveTetOp, class BelowTetOp>
        inline static void tetSliceWithPlane
        (
            const plane& pl,
            const tetPoints& tet,
            AboveTetOp& aboveOp,
            BelowTetOp& belowOp
        );


public:

    // Member Constants

        enum
        {
            nVertices = 4,  // Number of vertices in tetrahedron
            nEdges = 6      // Number of edges in tetrahedron
        };


    // Constructors

        //- Construct from four points
        inline tetrahedron
        (
            const Point& p0,
            const Point& p1,
            const Point& p2,
            const Point& p3
        );

        //- Construct from four points
        inline tetrahedron(const FixedList<Point, 4>& pts);

        //- Construct from four points in the list of points
        inline tetrahedron
        (
            const UList<Point>& points,
            const FixedList<label, 4>& indices
        );

        //- Construct from Istream
        inline explicit tetrahedron(Istream&);


    // Member Functions

    // Access

        //- Return vertex a
        const Point& a() const noexcept { return a_; }

        //- Return vertex b
        const Point& b() const noexcept { return b_; }

        //- Return vertex c
        const Point& c() const noexcept { return c_; }

        //- Return vertex d
        const Point& d() const noexcept { return d_; }

        //- Return i-th face
        inline triPointRef tri(const label facei) const;


    // Tet properties (static calculations)

        //- The enclosing (bounding) box for four points
        inline static Pair<Point> box
        (
            const Point& p0,
            const Point& p1,
            const Point& p2,
            const Point& p3
        );


        // Properties

            //- Face area normal for side a
            inline vector Sa() const;

            //- Face area normal for side b
            inline vector Sb() const;

            //- Face area normal for side c
            inline vector Sc() const;

            //- Face area normal for side d
            inline vector Sd() const;

            //- Return centre (centroid)
            inline Point centre() const;

            //- Return volume
            inline scalar mag() const;

            //- The enclosing (bounding) box for the tetrahedron
            inline Pair<Point> box() const;

            //- The bounding box for the tetrahedron
            inline treeBoundBox bounds() const;

            //- Return circum-centre
            inline Point circumCentre() const;

            //- Return circum-radius
            inline scalar circumRadius() const;

            //- Return quality: Ratio of tetrahedron and circum-sphere
            //- volume, scaled so that a regular tetrahedron has a
            //- quality of 1
            inline scalar quality() const;

            //- Return a random point in the tetrahedron from a
            //- uniform distribution
            inline Point randomPoint(Random& rndGen) const;

            //- Calculate the point from the given barycentric coordinates.
            inline Point barycentricToPoint(const barycentric& bary) const;

            //- Calculate the barycentric coordinates from the given point
            inline barycentric pointToBarycentric(const point& pt) const;

            //- Calculate the barycentric coordinates from the given point.
            //  Returns the determinant.
            inline scalar pointToBarycentric
            (
                const point& pt,
                barycentric& bary
            ) const;

            //- Return nearest point to p on tetrahedron. Is p itself
            //  if inside.
            inline pointHit nearestPoint(const point& p) const;

            //- Return true if point is inside tetrahedron
            inline bool inside(const point& pt) const;

            //- Decompose tet into tets above and below plane
            template<class AboveTetOp, class BelowTetOp>
            inline void sliceWithPlane
            (
                const plane& pl,
                AboveTetOp& aboveOp,
                BelowTetOp& belowOp
            ) const;

            //- Decompose tet into tets inside and outside other tet
            inline void tetOverlap
            (
                const tetrahedron<Point, PointRef>& tetB,
                tetIntersectionList& insideTets,
                label& nInside,
                tetIntersectionList& outsideTets,
                label& nOutside
            ) const;


            //- Return (min)containment sphere, i.e. the smallest sphere with
            //  all points inside. Returns pointHit with:
            //  - hit         : if sphere is equal to circumsphere
            //                  (biggest sphere)
            //  - point       : centre of sphere
            //  - distance    : radius of sphere
            //  - eligiblemiss: false
            // Tol (small compared to 1, e.g. 1e-9) is used to determine
            // whether point is inside: mag(pt - ctr) < (1+tol)*radius.
            pointHit containmentSphere(const scalar tol) const;

            //- Fill buffer with shape function products
            void gradNiSquared(scalarField& buffer) const;

            void gradNiDotGradNj(scalarField& buffer) const;

            void gradNiGradNi(tensorField& buffer) const;

            void gradNiGradNj(tensorField& buffer) const;


    // IOstream Operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            tetrahedron&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const tetrahedron&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetrahedronI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tetrahedron.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
