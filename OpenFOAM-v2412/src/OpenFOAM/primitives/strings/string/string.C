/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "string.H"
#include "stringOps.H"
#include "word.H"
#include "wordRe.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

const char* const Foam::string::typeName = "string";

int Foam::string::debug(Foam::debug::debugSwitch(Foam::string::typeName, 0));

const Foam::string Foam::string::null;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::string::ext() const
{
    const auto i = find_ext();

    if (i == npos)
    {
        return word();
    }

    return substr(i+1);
}


bool Foam::string::ext(const word& ending)
{
    if (ending.empty() || empty() || back() == '/')
    {
        return false;
    }
    else if (ending[0] == '.')
    {
        if (ending.size() == 1)
        {
            return false;
        }
    }
    else
    {
        append(1u, '.');
    }
    append(ending);

    return true;
}


bool Foam::string::has_ext(const wordRe& ending) const
{
    if (ending.empty() || ending.isLiteral())
    {
        return has_ext(static_cast<const std::string&>(ending));
    }

    const auto i = find_ext();
    if (i == npos)
    {
        return false;
    }

    // Regex match - compare *after* the dot
    return ending.match(substr(i+1));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::string::size_type Foam::string::count(const char c) const
{
    return stringOps::count(*this, c);
}


Foam::string& Foam::string::replace
(
    const std::string& s1,
    const std::string& s2,
    size_type pos
)
{
    if ((pos = find(s1, pos)) != npos)
    {
        std::string::replace(pos, s1.size(), s2);
    }

    return *this;
}


Foam::string& Foam::string::replaceAll
(
    const std::string& s1,
    const std::string& s2,
    size_type pos
)
{
    const auto n1 = s1.length();
    const auto n2 = s2.length();

    if (n1)
    {
        while ((pos = find(s1, pos)) != npos)
        {
            std::string::replace(pos, n1, s2);
            pos += n2;
        }
    }

    return *this;
}


Foam::string& Foam::string::replaceAny
(
    const std::string& s1,
    const char c2,
    size_type pos
)
{
    if (s1.length())
    {
        while ((pos = find_first_of(s1, pos)) != npos)
        {
            if (c2)
            {
                operator[](pos) = c2;
                ++pos;
            }
            else
            {
                erase(pos, 1);
            }
        }
    }

    return *this;
}


Foam::string& Foam::string::expand(const bool allowEmpty)
{
    stringOps::inplaceExpand(*this, allowEmpty);
    return *this;
}


bool Foam::string::removeRepeated(const char character)
{
    bool changed = false;

    if (character && find(character) != npos)
    {
        string::size_type nChar = 0;
        iterator outIter = begin();

        char prev = 0;

        for (auto iter = cbegin(); iter != cend(); ++iter)
        {
            const char c = *iter;

            if (prev == c && c == character)
            {
                changed = true;
            }
            else
            {
                *outIter = prev = c;
                ++outIter;
                ++nChar;
            }
        }

        erase(nChar);
    }

    return changed;
}


bool Foam::string::removeStart(const std::string& text)
{
    const auto txtLen = text.length();
    const auto strLen = length();

    if (txtLen && strLen >= txtLen && !compare(0, txtLen, text))
    {
        erase(0, txtLen);
        return true;
    }

    return false;
}


bool Foam::string::removeEnd(const std::string& text)
{
    const auto txtLen = text.length();
    const auto strLen = length();

    if (txtLen && strLen >= txtLen && !compare(strLen - txtLen, npos, text))
    {
        erase(strLen - txtLen);
        return true;
    }

    return false;
}


bool Foam::string::removeStart(const char c)
{
    if (length() > 1 && front() == c)
    {
        erase(0, 1);
        return true;
    }

    return false;
}


bool Foam::string::removeEnd(const char c)
{
    const auto n = length();
    if (n > 1 && back() == c)
    {
        erase(n-1);
        return true;
    }

    return false;
}


// ************************************************************************* //
