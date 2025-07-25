/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "keyType.H"
#include "regExp.H"
#include "token.H"
#include "IOstreams.H"
#include <utility>  // For std::swap

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::keyType Foam::keyType::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::keyType::keyType(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::keyType::swap(keyType& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-swap is a no-op
    }

    word::swap(static_cast<word&>(rhs));
    std::swap(type_, rhs.type_);
}


bool Foam::keyType::match(const std::string& text, bool literal) const
{
    if (!literal && isPattern())
    {
        return regExp(*this).match(text);  // Match as regex
    }

    return !compare(text);  // Compare as literal
}


bool Foam::keyType::assign(const token& tok)
{
    if (tok.isWord())
    {
        // Assign from word - literal
        assign(tok.wordToken());
        uncompile();
        return true;
    }
    else if (tok.isQuotedString())
    {
        // Assign from quoted string - regular expression
        assign(tok.stringToken());
        compile();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, keyType& val)
{
    token tok(is);

    if (val.assign(tok))
    {
        if (val.empty())
        {
            // Empty strings are an error
            FatalIOErrorInFunction(is)
                << "Zero-length regex"
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
    }
    else
    {
        FatalIOErrorInFunction(is);
        if (tok.good())
        {
            FatalIOError
                << "Wrong token type - expected word or string, found "
                << tok.info();
        }
        else
        {
            FatalIOError
                << "Bad token - could not get keyType";
        }
        FatalIOError << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const keyType& val)
{
    os.writeQuoted(val, val.isPattern());
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
