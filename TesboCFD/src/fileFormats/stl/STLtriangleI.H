/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "triangle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::STLtriangle::STLtriangle
(
    const STLpoint& normal,
    const STLpoint& p0,
    const STLpoint& p1,
    const STLpoint& p2,
    uint16_t attrib
)
:
    normal_(normal),
    a_(p0),
    b_(p1),
    c_(p2),
    attrib_(attrib)
{}


inline Foam::STLtriangle::STLtriangle(std::istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::STLtriangle::read(std::istream& is)
{
    is.read(reinterpret_cast<char*>(&normal_), 4*sizeof(STLpoint));
    is.read(reinterpret_cast<char*>(&attrib_), sizeof(STLattrib));
}


inline void Foam::STLtriangle::write(std::ostream& os) const
{
    os.write(reinterpret_cast<const char*>(&normal_), 4*sizeof(STLpoint));
    os.write(reinterpret_cast<const char*>(&attrib_), sizeof(STLattrib));
}


inline Foam::Ostream& Foam::STLtriangle::print(Ostream& os) const
{
    os  << " facet normal "
        << normal_.x() << ' ' << normal_.y() << ' ' << normal_.z() << nl
        << "  outer loop" << nl
        << "   vertex " << a_.x() << ' ' << a_.y() << ' ' << a_.z() << nl
        << "   vertex " << b_.x() << ' ' << b_.y() << ' ' << b_.z() << nl
        << "   vertex " << c_.x() << ' ' << c_.y() << ' ' << c_.z() << nl
        << "  endloop" << nl
        << " endfacet" << nl;

    return os;
}


inline void Foam::STLtriangle::write
(
    Ostream& os,
    const vector& norm,
    const point& p0,
    const point& p1,
    const point& p2
)
{
    os  << " facet normal "
        << norm.x() << ' ' << norm.y() << ' ' << norm.z() << nl
        << "  outer loop" << nl
        << "   vertex " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl
        << "   vertex " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl
        << "   vertex " << p2.x() << ' ' << p2.y() << ' ' << p2.z() << nl
        << "  endloop" << nl
        << " endfacet" << nl;
}


inline void Foam::STLtriangle::write
(
    Ostream& os,
    const point& p0,
    const point& p1,
    const point& p2
)
{
    // Calculate the normal ourselves
    const vector norm = triPointRef::unitNormal(p0, p1, p2);

    write(os, norm, p0, p1, p2);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const STLtriangle& tri)
{
    os  << tri.normal() << token::SPACE
        << tri.a() << token::SPACE
        << tri.b() << token::SPACE
        << tri.c() << token::SPACE
        << tri.attrib();

    return os;
}


// ************************************************************************* //
