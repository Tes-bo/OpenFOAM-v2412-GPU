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

\*---------------------------------------------------------------------------*/

#include "zero.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::linePoints::linePoints(const linePointRef& pts)
:
    Pair<point>(pts.a(), pts.b())
{}


inline Foam::linePoints::linePoints
(
    const UList<point>& points,
    const FixedList<label, 2>& indices
)
:
    Pair<point>(points[indices.get<0>()], points[indices.get<1>()])
{}


template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line
(
    const Point& from,
    const Point& to
)
:
    a_(from),
    b_(to)
{}


template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line
(
    const UList<Point>& points,
    const FixedList<label, 2>& indices
)
:
    a_(points[indices.template get<0>()]),
    b_(points[indices.template get<1>()])
{}


template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::linePointRef Foam::linePoints::ln() const
{
    return linePointRef(a(), b());
}


template<class Point, class PointRef>
inline Point Foam::line<Point, PointRef>::centre() const
{
    return 0.5*(a_ + b_);
}


inline Foam::point Foam::linePoints::centre() const
{
    return 0.5*(a() + b());
}


template<class Point, class PointRef>
inline Foam::scalar Foam::line<Point, PointRef>::mag() const
{
    return ::Foam::mag(b() - a());
}


inline Foam::scalar Foam::linePoints::mag() const
{
    return ::Foam::mag(b() - a());
}


template<class Point, class PointRef>
inline Point Foam::line<Point, PointRef>::vec() const
{
    return (b() - a());
}


inline Foam::vector Foam::linePoints::vec() const
{
    return (b() - a());
}


template<class Point, class PointRef>
inline Point Foam::line<Point, PointRef>::unitVec() const
{
    const Point v = (b_ - a_);

    #ifdef __clang__
    volatile  // Use volatile to avoid aggressive branch optimization
    #endif
    const scalar s(::Foam::mag(v));

    return s < ROOTVSMALL ? Zero : v/s;
}


inline Foam::vector Foam::linePoints::unitVec() const
{
    return normalised(b() - a());
}


template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::line<Point, PointRef>::box
(
    const Point& p0,
    const Point& p1
)
{
    return Pair<Point>(min(p0, p1), max(p0, p1));
}


template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::line<Point, PointRef>::box() const
{
    return line<Point, PointRef>::box(a_, b_);
}


inline Foam::Pair<Foam::point> Foam::linePoints::box() const
{
    return linePointRef::box(a(), b());
}


template<class Point, class PointRef>
Foam::PointHit<Point> Foam::line<Point, PointRef>::nearestDist
(
    const Point& p
) const
{
    Point v = vec();

    Point w(p - a_);

    const scalar c1 = v & w;

    if (c1 <= 0)
    {
        return PointHit<Point>(false, a_, Foam::mag(p - a_), true);
    }

    const scalar c2 = v & v;

    if (c2 <= c1)
    {
        return PointHit<Point>(false, b_, Foam::mag(p - b_), true);
    }

    const scalar b = c1/c2;

    Point pb(a_ + b*v);

    return PointHit<Point>(true, pb, Foam::mag(p - pb), false);
}


template<class Point, class PointRef>
Foam::scalar Foam::line<Point, PointRef>::nearestDist
(
    const line<Point, const Point&>& edge,
    Point& thisPt,
    Point& edgePt
) const
{
    // From Mathworld Line-Line distance/(Gellert et al. 1989, p. 538).
    Point a(end() - start());
    Point b(edge.end() - edge.start());
    Point c(edge.start() - start());

    Point crossab = a ^ b;
    const scalar magCrossSqr = Foam::magSqr(crossab);

    if (magCrossSqr > VSMALL)
    {
        scalar s = ((c ^ b) & crossab)/magCrossSqr;
        scalar t = ((c ^ a) & crossab)/magCrossSqr;

        // Check for end points outside of range 0..1
        if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
        {
            // Both inside range 0..1
            thisPt = start() + a*s;
            edgePt = edge.start() + b*t;
        }
        else
        {
            // Do brute force. Distance of everything to everything.
            // Can quite possibly be improved!

            // From edge endpoints to *this
            PointHit<Point> this0(nearestDist(edge.start()));
            PointHit<Point> this1(nearestDist(edge.end()));
            scalar thisDist = min(this0.distance(), this1.distance());

            // From *this to edge
            PointHit<Point> edge0(edge.nearestDist(start()));
            PointHit<Point> edge1(edge.nearestDist(end()));
            scalar edgeDist = min(edge0.distance(), edge1.distance());

            if (thisDist < edgeDist)
            {
                if (this0.distance() < this1.distance())
                {
                    thisPt = this0.point();
                    edgePt = edge.start();
                }
                else
                {
                    thisPt = this1.point();
                    edgePt = edge.end();
                }
            }
            else
            {
                if (edge0.distance() < edge1.distance())
                {
                    thisPt = start();
                    edgePt = edge0.point();
                }
                else
                {
                    thisPt = end();
                    edgePt = edge1.point();
                }
            }
        }
    }
    else
    {
        // Parallel lines. Find overlap of both lines by projecting onto
        // direction vector (now equal for both lines).

        scalar edge0 = edge.start() & a;
        scalar edge1 = edge.end() & a;
        bool edgeOrder = edge0 < edge1;

        scalar minEdge = (edgeOrder ? edge0 : edge1);
        scalar maxEdge = (edgeOrder ? edge1 : edge0);
        const Point& minEdgePt = (edgeOrder ? edge.start() : edge.end());
        const Point& maxEdgePt = (edgeOrder ? edge.end() : edge.start());

        scalar this0 = start() & a;
        scalar this1 = end() & a;
        bool thisOrder = this0 < this1;

        scalar minThis = min(this0, this1);
        scalar maxThis = max(this1, this0);
        const Point& minThisPt = (thisOrder ? start() : end());
        const Point& maxThisPt = (thisOrder ? end() : start());

        if (maxEdge < minThis)
        {
            // edge completely below *this
            edgePt = maxEdgePt;
            thisPt = minThisPt;
        }
        else if (maxEdge < maxThis)
        {
            // maxEdge inside interval of *this
            edgePt = maxEdgePt;
            thisPt = nearestDist(edgePt).point();
        }
        else
        {
            // maxEdge outside. Check if minEdge inside.
            if (minEdge < minThis)
            {
                // Edge completely envelops this. Take any this point and
                // determine nearest on edge.
                thisPt = minThisPt;
                edgePt = edge.nearestDist(thisPt).point();
            }
            else if (minEdge < maxThis)
            {
                // minEdge inside this interval.
                edgePt = minEdgePt;
                thisPt = nearestDist(edgePt).point();
            }
            else
            {
                // minEdge outside this interval
                edgePt = minEdgePt;
                thisPt = maxThisPt;
            }
        }
    }

    return Foam::mag(thisPt - edgePt);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    line<Point, PointRef>& l
)
{
    is.readBegin("line");
    is  >> l.a_ >> l.b_;
    is.readEnd("line");

    is.check(FUNCTION_NAME);
    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const line<Point, PointRef>& l
)
{
    os  << token::BEGIN_LIST
        << l.a_ << token::SPACE
        << l.b_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
