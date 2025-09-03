/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef, class PolygonRef>
inline Foam::pyramid<Point, PointRef, PolygonRef>::pyramid
(
    PolygonRef base,
    const Point& apex
)
:
    base_(base),
    apex_(apex)
{}


template<class Point, class PointRef, class PolygonRef>
inline Foam::pyramid<Point, PointRef, PolygonRef>::pyramid(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef, class PolygonRef>
inline Point Foam::pyramid<Point, PointRef, PolygonRef>::centre
(
    const UList<point>& points
) const
{
    return (3.0/4.0)*base_.centre(points) + (1.0/4.0)*apex_;
}


template<class Point, class PointRef, class PolygonRef>
inline Foam::vector Foam::pyramid<Point, PointRef, PolygonRef>::height
(
    const UList<point>& points
) const
{
    // Height = apex - baseCentroid
    return (apex_ - base_.centre(points));
}


template<class Point, class PointRef, class PolygonRef>
inline Foam::scalar Foam::pyramid<Point, PointRef, PolygonRef>::mag
(
    const UList<point>& points
) const
{
    return (1.0/3.0)*(base_.areaNormal(points) & (height(points)));
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Point, class PointRef, class PolygonRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    pyramid<Point, PointRef, PolygonRef>& p
)
{
    is  >> p.base_ >> p.apex_;
    is.check(FUNCTION_NAME);
    return is;
}


template<class Point, class PointRef, class PolygonRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pyramid<Point, PointRef, PolygonRef>& p
)
{
    os  << p.base_ << tab << p.apex_ << nl;
    return os;
}


// ************************************************************************* //
