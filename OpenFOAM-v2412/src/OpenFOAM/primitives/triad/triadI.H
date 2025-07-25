/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

inline Foam::triad::triad()
:
    Vector<vector>(triad::unset)
{}


inline Foam::triad::triad(const Vector<vector>& vv)
:
    Vector<vector>(vv)
{}


inline Foam::triad::triad(const vector& x, const vector& y, const vector& z)
:
    Vector<vector>(x, y, z)
{}


inline Foam::triad::triad(const tensor& t)
:
    Vector<vector>(t.x(), t.y(), t.z())
{}


inline Foam::triad::triad(const vector& pa)
{
    operator=(triad::unset);
    operator[](primaryDirection(pa)) = pa;
}


inline Foam::triad::triad(Istream& is)
:
    Vector<vector>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::triad::set(const direction d) const
{
    return operator[](d)[0] < GREAT;
}


inline bool Foam::triad::set() const
{
    return set(0) && set(1) && set(2);
}


inline Foam::direction Foam::triad::primaryDirection(const vector& v)
{
    if
    (
        Foam::mag(v.x()) > Foam::mag(v.y())
     && Foam::mag(v.x()) > Foam::mag(v.z())
    )
    {
        return triad::X;
    }
    else if (Foam::mag(v.y()) > Foam::mag(v.z()))
    {
        return triad::Y;
    }
    else
    {
        return triad::Z;
    }
}


inline Foam::vector Foam::triad::orthogonal
(
    const vector& v1,
    const vector& v2
)
{
    vector v3 = v1 ^ v2;

    scalar magV3 = Foam::mag(v3);

    if (magV3 > 0.5)
    {
        return v3/magV3;
    }
    else
    {
        return triad::unset[0];
    }
}


inline void Foam::triad::normalise()
{
    if (set(0)) (*this)[0].normalise();
    if (set(1)) (*this)[1].normalise();
    if (set(2)) (*this)[2].normalise();
}


inline Foam::vector Foam::triad::cx() const
{
    return vector(x().x(), y().x(), z().x());
}


inline Foam::vector Foam::triad::cy() const
{
    return vector(x().y(), y().y(), z().y());
}


inline Foam::vector Foam::triad::cz() const
{
    return vector(x().z(), y().z(), z().z());
}


inline Foam::triad Foam::triad::T() const
{
    return triad(cx(), cy(), cz());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::triad::operator=(const Vector<vector>& vv)
{
    Vector<vector>::operator=(vv);
}


inline void Foam::triad::operator=(const tensor& t)
{
    x() = t.x();
    y() = t.y();
    z() = t.z();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, triad& t)
{
    is >> static_cast<Vector<vector>&>(t);
    return is;
}


inline Foam::Ostream& Foam::operator<<(Ostream& os, const triad& t)
{
    os << static_cast<const Vector<vector>&>(t);
    return os;
}


// ************************************************************************* //
