/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::line

Description
    A line primitive.

SourceFiles
    lineI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_line_H
#define Foam_line_H

#include "point.H"
#include "point2D.H"
#include "vector.H"
#include "pointHit.H"
#include "FixedList.H"
#include "UList.H"
#include "Pair.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Point, class PointRef> class line;

template<class Point, class PointRef>
inline Istream& operator>>(Istream& is, line<Point, PointRef>& l);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream& os, const line<Point, PointRef>& l);


// Common Typedefs

//- A line using referred points
typedef line<point, const point&> linePointRef;


/*---------------------------------------------------------------------------*\
                         Class linePoints Declaration
\*---------------------------------------------------------------------------*/

//- Line point storage. Default constructable (line is not)
class linePoints
:
    public Pair<point>
{
public:

    // Generated Methods

        //- Default construct
        linePoints() = default;

        //- Inherit constructors
        using Pair<point>::Pair;


    // Constructors

        //- Construct from point references
        inline explicit linePoints(const linePointRef& pts);

        //- Copy construct from subset of points
        inline linePoints
        (
            const UList<point>& points,
            const FixedList<label, 2>& indices
        );


    // Member Functions

        //- The first vertex
        const point& a() const noexcept { return Pair<point>::first(); }

        //- The second vertex
        const point& b() const noexcept { return Pair<point>::second(); }

        //- The first vertex
        point& a() noexcept { return Pair<point>::first(); }

        //- The second vertex
        point& b() noexcept { return Pair<point>::second(); }

        //- Return as line reference
        inline linePointRef ln() const;


    // Properties

        //- Return centre (centroid)
        inline point centre() const;

        //- The magnitude (length) of the line
        inline scalar mag() const;

        //- Return start-to-end vector
        inline vector vec() const;

        //- Return the unit vector (start-to-end)
        inline vector unitVec() const;

        //- The enclosing (bounding) box for the line
        inline Pair<point> box() const;
};


/*---------------------------------------------------------------------------*\
                           Class line Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class line
{
    // Private Data

        //- Reference to the first line point
        PointRef a_;

        //- Reference to the second line point
        PointRef b_;


public:

    // Constructors

        //- Construct from two points
        inline line(const Point& from, const Point& to);

        //- Construct from two points in the list of points
        //  The indices could be from edge etc.
        inline line
        (
            const UList<Point>& points,
            const FixedList<label, 2>& indices
        );

        //- Construct from Istream
        inline explicit line(Istream& is);


    // Member Functions

    // Access

        //- The first point
        PointRef a() const noexcept { return a_; }

        //- The second point
        PointRef b() const noexcept { return b_; }

        //- The first point
        PointRef first() const noexcept { return a_; }

        //- The second (last) point
        PointRef second() const noexcept { return b_; }

        //- The start (first) point
        PointRef start() const noexcept { return a_; }

        //- The end (second) point
        PointRef end() const noexcept { return b_; }

        //- The last (second) point
        PointRef last() const noexcept { return b_; }


    // Line properties (static calculations)

        //- The enclosing (bounding) box for two points
        inline static Pair<Point> box(const Point& p0, const Point& p1);


    // Properties

        //- Return centre (centroid)
        inline Point centre() const;

        //- The magnitude (length) of the line
        inline scalar mag() const;

        //- Return start-to-end vector
        inline Point vec() const;

        //- Return the unit vector (start-to-end)
        inline Point unitVec() const;

        //- The enclosing (bounding) box for the line
        inline Pair<Point> box() const;

        //- Return nearest distance to line from a given point
        //  If the nearest point is on the line, return a hit
        PointHit<Point> nearestDist(const Point& p) const;

        //- Return nearest distance from line to line.
        //- Returns distance and sets both points
        //- (one on *this, one on the provided linePointRef.
        scalar nearestDist
        (
            const line<Point, const Point&>& edge,
            Point& thisPoint,
            Point& edgePoint
        ) const;


    // IOstream Operators

        friend Istream& operator>> <Point, PointRef>(Istream&, line&);
        friend Ostream& operator<< <Point, PointRef>(Ostream&, const line&);
};


//- 2D specialisation
template<>
scalar line<point2D, const point2D&>::nearestDist
(
    const line<point2D, const point2D&>& edge,
    point2D& thisPoint,
    point2D& edgePoint
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
