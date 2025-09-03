/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
#include "debug.H"
#include <cctype>
#include <cstdint>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::word::typeName = "word";

int Foam::word::debug(Foam::debug::debugSwitch(word::typeName, 0));

const Foam::word Foam::word::null;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::word::validate(const std::string& s, const bool prefix)
{
    word out;
    out.resize(s.size() + (prefix ? 1 : 0));

    std::string::size_type len = 0;

    // As per validate, but optionally detect if the first character
    // is a digit, which we'd like to avoid having since this will
    // cause parse issues when read back later.
    for (auto iter = s.cbegin(); iter != s.cend(); ++iter)
    {
        const char c = *iter;

        if (word::valid(c))
        {
            if (!len && prefix && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                out[len++] = '_';
            }

            out[len++] = c;
        }
    }

    out.erase(len);

    return out;
}


Foam::word Foam::word::validate
(
    const char* first,
    const char* last,
    const bool prefix
)
{
    std::string::size_type len = (last - first) + (prefix ? 1 : 0);

    word out;
    out.resize(len);

    for (len=0; first != last; ++first)
    {
        const char c = *first;

        if (word::valid(c))
        {
            if (!len && prefix && isdigit(c))
            {
                // First valid character was a digit - prefix with '_'
                out[len++] = '_';
            }

            out[len++] = c;
        }
    }

    out.erase(len);

    return out;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::word Foam::operator&(const word& a, const word& b)
{
    if (a.size())
    {
        if (b.size())
        {
            // Two non-empty words: can concatenate and perform camel case
            word camelCase(a + b);
            camelCase[a.size()] = char(toupper(b[0]));

            return camelCase;
        }

        return a;
    }

    // Or, if the first string is empty (or both are empty)

    return b;
}


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::word Foam::name(const void* ptr)
{
    std::ostringstream buf;
    buf.setf(std::ios_base::hex, std::ios_base::basefield);

    buf << "0x";  // Same as setf(std::ios::showbase)
    buf << uintptr_t(ptr);

    return word(buf.str(), false);  // Needs no stripping
}


// ************************************************************************* //
