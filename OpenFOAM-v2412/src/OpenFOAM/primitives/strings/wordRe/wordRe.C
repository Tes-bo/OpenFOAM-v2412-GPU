/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "wordRe.H"
#include "keyType.H"
#include "token.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::wordRe Foam::wordRe::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wordRe::wordRe(const keyType& str)
:
    word(str, false),  // No stripping
    regexPtr_(nullptr)
{
    if (str.isPattern())
    {
        compile();
    }
}


Foam::wordRe::wordRe(Istream& is)
:
    word(),
    regexPtr_(nullptr)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wordRe::assign(const token& tok)
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
        // Assign from quoted string - auto-detect regex
        assign(tok.stringToken());
        compile(wordRe::DETECT);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::wordRe::operator=(const keyType& str)
{
    assign(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        uncompile();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, wordRe& val)
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
                << "Bad token - could not get wordRe";
        }
        FatalIOError << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const wordRe& val)
{
    os.writeQuoted(val, val.isPattern());
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
