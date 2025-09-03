/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "surfacePatch.H"
#include "surfZone.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfacePatch::surfacePatch()
:
    surfacePatch(-1)
{}


Foam::surfacePatch::surfacePatch(const label index)
:
    geometricSurfacePatch(word::null, index, word::null),
    size_(0),
    start_(0)
{}


Foam::surfacePatch::surfacePatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    geometricSurfacePatch(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfacePatch::surfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricSurfacePatch(name, dict, index),
    size_(dict.get<label>("nFaces")),
    start_(dict.get<label>("startFace"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfacePatch::write(Ostream& os) const
{
    os.beginBlock(name());

    geometricSurfacePatch::write(os);

    os.writeEntry("nFaces", size());
    os.writeEntry("startFace", start());

    os.endBlock();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::surfacePatch::operator Foam::surfZone() const
{
    return surfZone
    (
        this->name(),
        this->size(),
        this->start(),
        this->index(),
        this->geometricType()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==
(
    const surfacePatch& a,
    const surfacePatch& b
)
{
    return
    (
        (a.size() == b.size())
     && (a.start() == b.start())
     && (a.geometricType() == b.geometricType())
    );
}


bool Foam::operator!=
(
    const surfacePatch& a,
    const surfacePatch& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfacePatch& obj)
{
    os  << static_cast<const geometricSurfacePatch&>(obj) << token::SPACE
        << obj.size() << token::SPACE
        << obj.start();

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
