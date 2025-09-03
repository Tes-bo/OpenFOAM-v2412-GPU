/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::plane::plane()
:
    normal_(Zero),
    origin_(Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::plane::normal() const noexcept
{
    return normal_;
}


inline const Foam::point& Foam::plane::origin() const noexcept
{
    return origin_;
}


inline Foam::point& Foam::plane::origin() noexcept
{
    return origin_;
}


inline void Foam::plane::flip()
{
    normal_ = -normal_;
}


inline Foam::point Foam::plane::nearestPoint(const point& p) const
{
    return p - normal_*((p - origin_) & normal_);
}


inline Foam::scalar Foam::plane::distance(const point& p) const
{
    return mag(signedDistance(p));
}


inline Foam::scalar Foam::plane::signedDistance(const point& p) const
{
    return ((p - origin_) & normal_);
}


inline Foam::plane::side Foam::plane::whichSide(const point& p) const
{
    const scalar dist = signedDistance(p);

    return (dist < 0 ? side::BACK : side::FRONT);
}


inline int Foam::plane::sign(const point& p, const scalar tol) const
{
    const scalar dist = signedDistance(p);

    return ((dist < -tol) ? -1 : (dist > tol) ? +1 : 0);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const plane& a, const plane& b)
{
    return (a.origin() == b.origin() && a.normal() == b.normal());
}


inline bool Foam::operator!=(const plane& a, const plane& b)
{
    return !(a == b);
}


inline bool Foam::operator<(const plane& a, const plane& b)
{
    return (a.origin() < b.origin());
}


// ************************************************************************* //
