/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "SHA1Digest.H"
#include "IOstreams.H"
#include <cstring>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::SHA1Digest Foam::SHA1Digest::null;

static const char hexChars[] = "0123456789abcdef";

// The char '0' == 0
static constexpr int offsetZero = int('0');

// The char 'A' (or 'a') == 10
static constexpr int offsetAlpha = int('A') - 10;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read hexadecimal value, ignoring leading or intermediate '_'
static unsigned char readHexDigit(Istream& is)
{
    // Silently ignore leading or intermediate '_'
    char c = 0;
    do
    {
        is.read(c);
    }
    while (c == '_');

    if (isdigit(c))
    {
        return int(c) - offsetZero;
    }
    else if (!isxdigit(c))
    {
        FatalIOErrorInFunction(is)
            << "Illegal hex digit: '" << c << "'"
            << exit(FatalIOError);
    }

    return toupper(c) - offsetAlpha;
}

} // End namespace Foam


namespace
{

// Copy assign digest from content
bool assign
(
    std::array<unsigned char, 20>& digest,
    const unsigned char* content,
    std::size_t len
)
{
    if (!content || !len)
    {
        return false;
    }

    if (len == digest.size())
    {
        // ie, std::copy
        for (auto& val : digest)
        {
            val = *content;
            ++content;
        }

        return true;
    }

    // Skip possible '_' prefix
    if (*content == '_')
    {
        ++content;
        --len;
    }

    // Incorrect length - can never assign
    if (len != 2*digest.size())
    {
        return false;
    }

    for (auto& val : digest)
    {
        const unsigned char upp = *content++;
        const unsigned char low = *content++;

        val = (upp << 4) + low;
    }

    return true;
}


// Byte-wise compare digest contents
bool isEqual
(
    const std::array<unsigned char, 20>& digest,
    const char* hexdigits,
    std::size_t len
)
{
    // Skip possible '_' prefix
    if (*hexdigits == '_')
    {
        ++hexdigits;
        --len;
    }

    // Incorrect length - can never match
    if (len != 2*digest.size())
    {
        return false;
    }

    for (const auto& byteVal : digest)
    {
        const char upp = hexChars[((byteVal >> 4) & 0xF)];
        const char low = hexChars[(byteVal & 0xF)];

        if (upp != *hexdigits++) return false;
        if (low != *hexdigits++) return false;
    }

    return true;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SHA1Digest::SHA1Digest()
{
    clear();
}


Foam::SHA1Digest::SHA1Digest(const char* content, std::size_t len)
{
    clear();
    assign(dig_, reinterpret_cast<const unsigned char*>(content), len);
}


Foam::SHA1Digest::SHA1Digest(const unsigned char* content, std::size_t len)
{
    clear();
    assign(dig_, content, len);
}


Foam::SHA1Digest::SHA1Digest(Istream& is)
{
    clear();
    read(is);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SHA1Digest::clear()
{
    dig_.fill(0);  // Same as memset(dig_.data(), 0, dig_.size());
}


bool Foam::SHA1Digest::empty() const
{
    for (const auto& byteVal : dig_)
    {
        if (byteVal)
        {
            return false;
        }
    }

    return true;
}


Foam::Istream& Foam::SHA1Digest::read(Istream& is)
{
    for (auto& byteVal : dig_)
    {
        const unsigned char upp = readHexDigit(is);
        const unsigned char low = readHexDigit(is);

        byteVal = (upp << 4) + low;
    }

    is.check(FUNCTION_NAME);
    return is;
}


std::string Foam::SHA1Digest::str(const bool prefixed) const
{
    std::string buf;
    std::size_t nChar = 0;

    if (prefixed)
    {
        buf.resize(1 + 2*dig_.size());
        buf[nChar++] = '_';
    }
    else
    {
        buf.resize(2*dig_.size());
    }

    for (const auto& byteVal : dig_)
    {
        buf[nChar++] = hexChars[((byteVal >> 4) & 0xF)];  // Upper nibble
        buf[nChar++] = hexChars[(byteVal & 0xF)];         // Lower nibble
    }

    return buf;
}


Foam::Ostream& Foam::SHA1Digest::write(Ostream& os, const bool prefixed) const
{
    if (prefixed)
    {
        os.write('_');
    }

    for (const auto& byteVal : dig_)
    {
        os.write(hexChars[((byteVal >> 4) & 0xF)]);  // Upper nibble
        os.write(hexChars[(byteVal & 0xF)]);         // Lower nibble
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::SHA1Digest::operator==(const SHA1Digest& rhs) const
{
    return (dig_ == rhs.dig_);
}


bool Foam::SHA1Digest::operator==(const std::string& hexdigits) const
{
    // Interpret empty string as '0000..'
    size_t len = hexdigits.length();

    return len ? isEqual(dig_, hexdigits.data(), len) : empty();
}


bool Foam::SHA1Digest::operator==(const char* hexdigits) const
{
    // Interpret nullptr or empty string as '0000..'
    size_t len = (hexdigits ? strlen(hexdigits) : 0);

    return len ? isEqual(dig_, hexdigits, len) : empty();
}


bool Foam::SHA1Digest::operator!=(const SHA1Digest& rhs) const
{
    return !this->operator==(rhs);
}


bool Foam::SHA1Digest::operator!=(const std::string& hexdigits) const
{
    return !this->operator==(hexdigits);
}


bool Foam::SHA1Digest::operator!=(const char* hexdigits) const
{
    return !this->operator==(hexdigits);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, SHA1Digest& dig)
{
    return dig.read(is);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1Digest& dig)
{
    // Write with prefixed = false
    return dig.write(os, false);
}


// ************************************************************************* //
