/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "surfZone.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZone::surfZone()
:
    surfZoneIdentifier(),
    size_(0),
    start_(0)
{}


Foam::surfZone::surfZone(const word& name, const label size)
:
    surfZoneIdentifier(name, 0, word::null),
    size_(size),
    start_(0)
{}


Foam::surfZone::surfZone
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    surfZoneIdentifier(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfZone::surfZone
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfZoneIdentifier(name, dict, index),
    size_(dict.get<label>("nFaces")),
    start_(dict.get<label>("startFace"))
{}


Foam::surfZone::surfZone(const surfZone& zone, const label index)
:
    surfZone(zone)
{
    surfZoneIdentifier::index() = index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfZone::write(Ostream& os) const
{
    os.beginBlock(name());

    surfZoneIdentifier::write(os);
    os.writeEntry("nFaces", size());
    os.writeEntry("startFace", start());

    os.endBlock();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const surfZone& a, const surfZone& b)
{
    return
    (
        a.size()  == b.size()
     && a.start() == b.start()
     && a.geometricType() == b.geometricType()
    );
}


bool Foam::operator!=(const surfZone& a, const surfZone& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfZone& obj)
{
    const word name(is);
    const dictionary dict(is);

    // Could also leave index untouched?
    obj = surfZone(name, dict, 0);

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZone& obj)
{
    obj.write(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
