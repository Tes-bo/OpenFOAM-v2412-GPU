/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Istream.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The current get position (std::istream only)
inline std::streampos stream_tellg(Foam::Istream* isptr)
{
    auto* sptr = dynamic_cast<Foam::ISstream*>(isptr);
    return sptr ? sptr->stdStream().tellg() : std::streampos(0);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::token& Foam::Istream::peekBack() const noexcept
{
    return (putBackAvail_ ? putBackToken_ : token::undefinedToken);
}

// Return the putback token if available or fetch a new token
// from the stream.
//
// Foam::token& Foam::Istream::peekToken()
// {
//     if (!putBackAvail_)
//     {
//         putBackToken_.reset();
//         token tok;
//         this->read(tok);
//         putBackToken_ = std::move(tok);
//     }
//
//     return putBackToken_;
// }


void Foam::Istream::putBackClear()
{
    putBackAvail_ = false;
    putBackToken_.reset();
}


void Foam::Istream::putBack(const token& tok)
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back onto bad stream"
            << exit(FatalIOError);
    }
    else if (putBackAvail_)
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back another token"
            << exit(FatalIOError);
    }
    else
    {
        putBackAvail_ = true;
        putBackToken_ = tok;
    }
}


void Foam::Istream::putBack(token&& tok)
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back onto bad stream"
            << exit(FatalIOError);
    }
    else if (putBackAvail_)
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back another token"
            << exit(FatalIOError);
    }
    else
    {
        putBackAvail_ = true;
        putBackToken_ = std::move(tok);
    }
}


bool Foam::Istream::getBack(token& tok)
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to get back from bad stream"
            << exit(FatalIOError);
    }
    else if (putBackAvail_)
    {
        putBackAvail_ = false;
        tok = std::move(putBackToken_);
        return true;
    }

    return false;
}


bool Foam::Istream::readBegin(const char* funcName)
{
    const token delimiter(*this);

    if (delimiter != token::BEGIN_LIST)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::BEGIN_LIST
            << "' while reading " << funcName
            << ", found " << delimiter.info() << nl
            << exit(FatalIOError);
    }

    return true;
}


bool Foam::Istream::readEnd(const char* funcName)
{
    const token delimiter(*this);

    if (delimiter != token::END_LIST)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::END_LIST
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << " at stream position " << stream_tellg(this) << nl
            << exit(FatalIOError);
    }

    return true;
}


char Foam::Istream::readBeginList(const char* funcName)
{
    const token delimiter(*this);

    if (delimiter != token::BEGIN_LIST && delimiter != token::BEGIN_BLOCK)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::BEGIN_LIST
            << "' or a '" << token::BEGIN_BLOCK
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


char Foam::Istream::readEndList(const char* funcName)
{
    const token delimiter(*this);

    if (delimiter != token::END_LIST && delimiter != token::END_BLOCK)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::END_LIST
            << "' or a '" << token::END_BLOCK
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << " at stream position " << stream_tellg(this) << nl
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


Foam::Istream& Foam::Istream::operator()() const
{
    if (!good())
    {
        check(FUNCTION_NAME);
        FatalIOError.exit();
    }

    return const_cast<Istream&>(*this);
}


// ************************************************************************* //
