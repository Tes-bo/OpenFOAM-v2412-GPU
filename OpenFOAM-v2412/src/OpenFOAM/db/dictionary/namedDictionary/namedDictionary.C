/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "namedDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::namedDictionary::namedDictionary(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::namedDictionary::clear()
{
    first().clear();
    second().clear();
}


bool Foam::namedDictionary::empty() const noexcept
{
    return (first().empty() && second().empty());
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, namedDictionary& obj)
{
    obj.clear();

    // Three possible inputs:
    // - key
    // - key { ... }
    // - { ... }

    // Minor consistency with primitiveEntry, also accept the following:
    // - key ;

    token tok(is);
    is.putBack(tok);

    if (!tok.isPunctuation(token::BEGIN_BLOCK))
    {
        is >> obj.keyword();
        is >> tok;

        // Discards possible trailing ';'
        if (!tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    if (tok.isPunctuation(token::BEGIN_BLOCK))
    {
        obj.dict().read(is);
        // Provide a name as context for error messages etc.
        obj.dict().name() = is.relativeName();
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const namedDictionary& obj)
{
    // Three possible outputs:
    // - key
    // - key { ... }
    // - { ... }
    // No distinction between a missing and an empty dictionary

    if (obj.keyword().empty() || !obj.dict().empty())
    {
        // Never allow empty output.
        // Otherwise cannot re-read for streaming
        obj.dict().writeEntry(obj.keyword(), os);
    }
    else
    {
        os << obj.keyword();
    }

    return os;
}


// ************************************************************************* //
