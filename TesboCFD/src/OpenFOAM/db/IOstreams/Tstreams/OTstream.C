/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "error.H"
#include "OTstream.H"
#include <cctype>

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::OTstream::write(const token& tok)
{
    if (tok.good())
    {
        tokens().push_back(tok);
        return true;
    }

    return false;
}


Foam::Ostream& Foam::OTstream::write(const char c)
{
    if (!std::isspace(c) && std::isprint(c))
    {
        // Should generally work, but need to verify corner cases
        tokens().push_back(token(token::punctuationToken(c)));
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::writeQuoted
(
    const char* str,
    std::streamsize len,
    const bool quoted
)
{
    if (quoted)
    {
        // tokenType::STRING
        tokens().emplace_back() = string(str, len);
    }
    else if (len > 0)
    {
        // Create from std::string with specified type never strips
        tokens().emplace_back
        (
            token::tokenType::WORD,  // or perhaps tokenType::CHAR_DATA ?
            std::string(str, len)
        );
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* str)
{
    word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        // Like punctuation
        write(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        // As a word
        tokens().emplace_back() = std::move(nonWhiteChars);  // Move assign
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const word& str)
{
    // tokenType::WORD
    tokens().emplace_back() = str;  // Copy assign

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const std::string& str)
{
    // tokenType::STRING
    tokens().emplace_back() = Foam::string(str);  // Move assign

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int32_t val)
{
    tokens().push_back(token(label(val))); // tokenType::LABEL

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int64_t val)
{
    tokens().push_back(token(label(val))); // tokenType::LABEL

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const float val)
{
    tokens().push_back(token(val)); // tokenType::FLOAT

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const double val)
{
    tokens().push_back(token(val)); // tokenType::DOUBLE

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* data, std::streamsize count)
{
    // if (format() != IOstreamOption::BINARY)
    // {
    //     FatalErrorInFunction
    //         << "stream format not binary"
    //         << Foam::abort(FatalError);
    // }

    NotImplemented;
    return *this;
}


Foam::Ostream& Foam::OTstream::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for IOstreamOption::BINARY since this is either done in the
    // beginRawWrite() method, or the caller knows what they are doing.

    NotImplemented;
    return *this;
}


bool Foam::OTstream::beginRawWrite(std::streamsize count)
{
    // if (format() != IOstreamOption::BINARY)
    // {
    //     FatalErrorInFunction
    //         << "stream format not binary"
    //         << Foam::abort(FatalError);
    // }

    NotImplemented;
    return true;
}


void Foam::OTstream::print(Ostream& os) const
{
    os  << "OTstream : " << name().c_str() << ", " << size() << " tokens, ";
    IOstream::print(os);
}


// ************************************************************************* //
