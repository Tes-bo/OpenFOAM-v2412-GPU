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

#include "word.H"
#include "token.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::word::word(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, word& val)
{
    token tok(is);

    if (tok.isWord())
    {
        val = tok.wordToken();
    }
    else if (tok.isQuotedString() || tok.isVariable())
    {
        // Try a bit harder, convert some string types to word
        // - accept "quoted" or $tag, but not verbatim/expression

        const auto& str = tok.stringToken();
        val = word::validate(str);

        // Flag empty strings and bad chars as an error
        if (val.empty() || val.length() != str.length())
        {
            FatalIOErrorInFunction(is)
                << "Empty word or non-word characters "
                << tok.info() << exit(FatalIOError);
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
                << "Wrong token type - expected word, found "
                << tok.info();
        }
        else
        {
            FatalIOError
                << "Bad token - could not get word";
        }
        FatalIOError << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const word& val)
{
    os.write(val);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
