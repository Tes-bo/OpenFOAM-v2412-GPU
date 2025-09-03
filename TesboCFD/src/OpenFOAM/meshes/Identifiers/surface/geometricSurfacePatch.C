/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "geometricSurfacePatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricSurfacePatch::geometricSurfacePatch()
:
    geometricSurfacePatch(0)
{}


Foam::geometricSurfacePatch::geometricSurfacePatch(const label index)
:
    name_("patch"),
    index_(index),
    geometricType_()
{}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const label index
)
:
    name_(name),
    index_(index),
    geometricType_()
{}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const label index,
    const word& geometricType
)
:
    name_(name),
    index_(index),
    geometricType_(geometricType)
{}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricSurfacePatch(name, index)
{
    dict.readIfPresent("geometricType", geometricType_);
}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const surfZoneIdentifier& ident
)
:
    geometricSurfacePatch(ident.name(), ident.index(), ident.geometricType())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geometricSurfacePatch::write(Ostream& os) const
{
    if (!geometricType_.empty())
    {
        os.writeEntry("geometricType", geometricType_);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==
(
    const geometricSurfacePatch& a,
    const geometricSurfacePatch& b
)
{
    return
    (
        (a.name() == b.name())
     && (a.geometricType() == b.geometricType())
    );
}


bool Foam::operator!=
(
    const geometricSurfacePatch& a,
    const geometricSurfacePatch& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, geometricSurfacePatch& obj)
{
    // Also read "" for empty words
    obj.name() = word::validate(token(is).stringToken());
    obj.geometricType() = word::validate(token(is).stringToken());

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const geometricSurfacePatch& obj)
{
    // Force unconditional line-breaks on list output.
    // We otherwise risk extremely unreadable entries
    os << nl;

    // Empty words are double-quoted so they are treated as 'string'
    // and not simply lost

    os.writeQuoted(obj.name(), obj.name().empty()) << token::SPACE;
    os.writeQuoted(obj.geometricType(), obj.geometricType().empty());

    // Equivalent for OpenFOAM-1912 and earlier:
    // if (obj.geometricType().empty())
    // {
    //     os << emptyType;
    // }
    // else
    // {
    //     os << obj.geometricType();
    // }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
