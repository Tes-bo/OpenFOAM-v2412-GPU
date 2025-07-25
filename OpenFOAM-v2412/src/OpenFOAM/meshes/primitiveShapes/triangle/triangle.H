/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::triangle

Description
    A triangle primitive used to calculate face normals and swept volumes.
    Uses referred points.

SourceFiles
    triangleI.H
    triangle.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triangle_H
#define Foam_triangle_H

#include "triangleFwd.H"
#include "intersection.H"
#include "vector.H"
#include "tensor.H"
#include "pointHit.H"
#include "Random.H"
#include "FixedList.H"
#include "UList.H"
#include "line.H"
#include "Pair.H"
#include "Tuple2.H"
#include "barycentric2D.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class plane;

template<class Point, class PointRef>
inline Istream& operator>>(Istream&, triangle<Point, PointRef>&);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream&, const triangle<Point, PointRef>&);


/*---------------------------------------------------------------------------*\
                          Class triPoints Declaration
\*---------------------------------------------------------------------------*/

//- Triangle point storage. Default constructable (triangle is not)
class triPoints
:
    public FixedList<point, 3>
{
public:

    // Generated Methods

        //- Default construct
        triPoints() = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Construct from three points
        inline triPoints(const point& p0, const point& p1, const point& p2);

        //- Construct from point references
        inline explicit triPoints(const triPointRef& pts);

        //- Construct from three points
        inline triPoints(const FixedList<point, 3>& pts);

        //- Copy construct from subset of points
        inline triPoints
        (
            const UList<point>& points,
            const FixedList<label, 3>& indices
        );

        //- Copy construct from subset of points
        inline triPoints
        (
            const UList<point>& points,
            const label p0,
            const label p1,
            const label p2
        );


    // Member Functions

        //- The first vertex
        const point& a() const noexcept { return get<0>(); }

        //- The second vertex
        const point& b() const noexcept { return get<1>(); }

        //- The third vertex
        const point& c() const noexcept { return get<2>(); }

        //- The first vertex
        point& a() noexcept { return get<0>(); }

        //- The second vertex
        point& b() noexcept { return get<1>(); }

        //- The third vertex
        point& c() noexcept { return get<2>(); }

        //- Flip triangle orientation by swapping second and third vertices
        inline void flip();

        //- Return as triangle reference
        inline triPointRef tri() const;


    // Properties

        //- Return centre (centroid)
        inline point centre() const;

        //- The area normal - with magnitude equal to area of triangle
        inline vector areaNormal() const;

        //- Return unit normal
        inline vector unitNormal() const;

        //- The magnitude of the triangle area
        inline scalar mag() const;

        //- The enclosing (bounding) box for the triangle
        inline Pair<point> box() const;

        //- Edge vector opposite point a(): from b() to c()
        inline vector vecA() const;

        //- Edge vector opposite point b(): from c() to a()
        inline vector vecB() const;

        //- Edge vector opposite point c(): from a() to b()
        inline vector vecC() const;
};


/*---------------------------------------------------------------------------*\
                          Class triangle Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class triangle
{
public:

    // Public Typedefs

        //- The point type
        typedef Point point_type;

        //- Storage type for triangles originating from intersecting triangle
        //- with another triangle
        typedef FixedList<triPoints, 27> triIntersectionList;

        //- Proximity classifications
        enum proxType
        {
            NONE = 0,   //!< Unknown proximity
            POINT,      //!< Close to point
            EDGE        //!< Close to edge
        };


    // Public Classes

        // Classes for use in sliceWithPlane.
        // What to do with decomposition of triangle.

            //- Dummy
            class dummyOp
            {
            public:
                inline void operator()(const triPoints&);
            };

            //- Sum resulting areas
            class sumAreaOp
            {
            public:
                scalar area_;

                inline sumAreaOp();

                inline void operator()(const triPoints&);
            };

            //- Store resulting tris
            class storeOp
            {
            public:
                triIntersectionList& tris_;
                label& nTris_;

                inline storeOp(triIntersectionList&, label&);

                inline void operator()(const triPoints&);
            };


private:

    // Private Data

        //- Reference to the first triangle point
        PointRef a_;

        //- Reference to the second triangle point
        PointRef b_;

        //- Reference to the third triangle point
        PointRef c_;


   // Private Member Functions

        //- Helper: calculate intersection point
        inline static point planeIntersection
        (
            const FixedList<scalar, 3>& d,
            const triPoints& t,
            const label negI,
            const label posI
        );

        //- Helper: slice triangle with plane
        template<class AboveOp, class BelowOp>
        inline static void triSliceWithPlane
        (
            const plane& pln,
            const triPoints& tri,
            AboveOp& aboveOp,
            BelowOp& belowOp
        );


public:

    // Constructors

        //- Construct from three points
        inline triangle(const Point& p0, const Point& p1, const Point& p2);

        //- Construct from three points
        inline triangle(const FixedList<Point, 3>& pts);

        //- Construct from three points out of the list of points
        //  The indices could be from triFace etc.
        inline triangle
        (
            const UList<Point>& points,
            const FixedList<label, 3>& indices
        );

        //- Construct from three points out of the list of points
        inline triangle
        (
            const UList<Point>& points,
            const label p0,
            const label p1,
            const label p2
        );

        //- Construct from Istream
        inline explicit triangle(Istream& is);


    // Member Functions

    // Access

        //- The first vertex
        const Point& a() const noexcept { return a_; }

        //- The second vertex
        const Point& b() const noexcept { return b_; }

        //- The third vertex
        const Point& c() const noexcept { return c_; }


    // Triangle properties (static calculations)

        //- The centre (centroid) of three points
        inline static Point centre
        (
            const Point& p0,
            const Point& p1,
            const Point& p2
        );

        //- The area normal for a triangle defined by three points
        //- (right-hand rule). Magnitude equal to area of the triangle
        inline static vector areaNormal
        (
            const Point& p0,
            const Point& p1,
            const Point& p2
        );

        //- The unit normal for a triangle defined by three points
        //- (right-hand rule).
        inline static vector unitNormal
        (
            const Point& p0,
            const Point& p1,
            const Point& p2
        );

        //- The enclosing (bounding) box for three points
        inline static Pair<Point> box
        (
            const Point& p0,
            const Point& p1,
            const Point& p2
        );


    // Properties

        //- Return centre (centroid)
        inline Point centre() const;

        //- The area normal - with magnitude equal to area of triangle
        inline vector areaNormal() const;

        //- Return unit normal
        inline vector unitNormal() const;

        //- Legacy name for areaNormal().
        //  \deprecated(2018-06) Deprecated for new use
        FOAM_DEPRECATED_FOR(2018-12, "areaNormal() or unitNormal()")
        vector normal() const
        {
            return areaNormal();
        }

        //- The magnitude of the triangle area
        inline scalar mag() const;

        //- The enclosing (bounding) box for the triangle
        inline Pair<Point> box() const;

        //- Edge vector opposite point a(): from b() to c()
        inline Point vecA() const;

        //- Edge vector opposite point b(): from c() to a()
        inline Point vecB() const;

        //- Edge vector opposite point c(): from a() to b()
        inline Point vecC() const;


    // Properties

            //- Return circum-centre
            inline Point circumCentre() const;

            //- Return circum-radius
            inline scalar circumRadius() const;

            //- Return quality: Ratio of triangle and circum-circle
            //  area, scaled so that an equilateral triangle has a
            //  quality of 1
            inline scalar quality() const;

            //- Return swept-volume
            inline scalar sweptVol(const triangle& t) const;

            //- Return the inertia tensor, with optional reference
            //  point and density specification
            inline tensor inertia
            (
                PointRef refPt = Zero,
                scalar density = 1.0
            ) const;

            //- Return a random point on the triangle from a uniform
            //- distribution
            inline Point randomPoint(Random& rndGen) const;

            //- Calculate the point from the given barycentric coordinates.
            inline Point barycentricToPoint(const barycentric2D& bary) const;

            //- Calculate the barycentric coordinates from the given point
            inline barycentric2D pointToBarycentric(const point& pt) const;

            //- Calculate the barycentric coordinates from the given point.
            //  Returns the determinant.
            inline scalar pointToBarycentric
            (
                const point& pt,
                barycentric2D& bary
            ) const;

            //- Fast intersection detection with a plane.
            inline bool intersects
            (
                const point& origin,
                const vector& normal
            ) const;

            //- Fast intersection detection with an \b axis plane.
            inline bool intersects
            (
                //! Origin of the plane
                const point& origin,
                //! Normal of the plane (vector::X, vector::Y, vector::Z)
                const vector::components axis
            ) const;

            //- Return point intersection with a ray.
            //  For a hit, the distance is signed. Positive number
            //  represents the point in front of triangle.
            //  In case of miss pointHit is set to nearest point
            //  on triangle and its distance to the distance between
            //  the original point and the plane intersection point
            inline pointHit ray
            (
                const point& p,
                const vector& q,
                const intersection::algorithm = intersection::FULL_RAY,
                const intersection::direction dir = intersection::VECTOR
            ) const;

            //- Fast intersection with a ray.
            //  For a hit, the pointHit.distance() is the line parameter t :
            //  intersection=p+t*q. Only defined for VISIBLE, FULL_RAY or
            //  HALF_RAY. tol increases the virtual size of the triangle
            // by a relative factor.
            inline pointHit intersection
            (
                const point& p,
                const vector& q,
                const intersection::algorithm alg,
                const scalar tol = 0.0
            ) const;

            //- Find the nearest point to p on the triangle and classify it:
            //  + near point (nearType=POINT, nearLabel=0, 1, 2)
            //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
            //    Note: edges are counted from starting
            //    vertex so e.g. edge 2 is from f[2] to f[0]
            pointHit nearestPointClassify
            (
                const point& p,
                label& nearType,
                label& nearLabel
            ) const;

            //- Return nearest point to p on triangle
            inline pointHit nearestPoint(const point& p) const;

            //- Classify nearest point to p in triangle plane
            //  w.r.t. triangle edges and points.  Returns inside
            //  (true)/outside (false).
            bool classify
            (
                const point& p,
                label& nearType,
                label& nearLabel
            ) const;

            //- Return nearest point to line on triangle. Returns hit if
            //  point is inside triangle. Sets edgePoint to point on edge
            //  (hit if nearest is inside line)
            inline pointHit nearestPoint
            (
                const linePointRef& edge,
                pointHit& edgePoint
            ) const;

            //- The sign for which side of the face plane the point is on.
            //  Uses the supplied tolerance for rounding around zero.
            //  \return
            //  -  0: on plane
            //  - +1: above plane
            //  - -1: below plane
            inline int sign(const point& p, const scalar tol = SMALL) const;

            //- Decompose triangle into triangles above and below plane
            template<class AboveOp, class BelowOp>
            inline void sliceWithPlane
            (
                const plane& pln,
                AboveOp& aboveOp,
                BelowOp& belowOp
            ) const;

            //- Decompose triangle into triangles inside and outside
            //  (with respect to user provided normal) other
            //  triangle.
            template<class InsideOp, class OutsideOp>
            inline void triangleOverlap
            (
                const vector& n,
                const triangle<Point, PointRef>& tri,
                InsideOp& insideOp,
                OutsideOp& outsideOp
            ) const;


    // IOstream Operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            triangle&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const triangle&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triangleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "triangle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
