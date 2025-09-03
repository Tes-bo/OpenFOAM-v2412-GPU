/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "regExpPosix.H"
#include "SubStrings.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::regExpPosix::grammar(0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Matched entire length
static inline bool fullMatch(const regmatch_t& m, const regoff_t len)
{
    return (m.rm_so == 0 && m.rm_eo == len);
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::regExpPosix::set_pattern
(
    const char* pattern,
    size_t len,
    bool ignoreCase
)
{
    clear();  // Also sets ctrl_ = 0

    const char* pat = pattern;
    bool doNegate = false;

    // Handle known embedded prefixes
    if (len > 2 && pat[0] == '(' && pat[1] == '?')
    {
        pat += 2;
        len -= 2;

        for (bool done = false; !done && len; ++pat, --len)
        {
            switch (*pat)
            {
                case '!':
                {
                    // Negated (inverted) match
                    doNegate = true;
                    break;
                }
                case 'i':
                {
                    // Ignore-case
                    ignoreCase = true;
                    break;
                }
                case ')':
                {
                    // End of prefix parsing
                    done = true;
                    break;
                }
            }
        }
    }

    // Avoid zero-length patterns
    if (len)
    {
        int flags = REG_EXTENDED;
        if (ignoreCase)
        {
            flags |= REG_ICASE;
        }

        {
            preg_ = new regex_t;
            int err = regcomp(preg_, pat, flags);

            if (err == 0)
            {
                ctrl_ = (doNegate ? ctrlType::NEGATED : ctrlType::NORMAL);
                return true;
            }
            else
            {
                char errbuf[200];
                regerror(err, preg_, errbuf, sizeof(errbuf));

                FatalErrorInFunction
                    << "Failed to compile regular expression '"
                    << pattern << "'\n" << errbuf
                    << exit(FatalError);
            }
        }
    }

    return false;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::regExpPosix::clear()
{
    ctrl_ = 0;

    if (preg_)
    {
        regfree(preg_);
        delete preg_;
        preg_ = nullptr;
        return true;
    }

    return false;
}


std::string::size_type Foam::regExpPosix::find(const std::string& text) const
{
    // Find with negated is probably not very reliable...
    if (!preg_ || !ctrl_)
    {
        // Undefined: never matches
        return std::string::npos;
    }
    else if (text.empty())
    {
        if (ctrl_ == ctrlType::NEGATED)
        {
            return 0;  // No match - pretend it starts at position 0
        }
        else
        {
            return std::string::npos;
        }
    }
    else
    {
        const size_t nmatch = 1;
        regmatch_t pmatch[1];

        const bool ok = (regexec(preg_, text.c_str(), nmatch, pmatch, 0) == 0);

        if (ctrl_ == ctrlType::NEGATED)
        {
            if (!ok)
            {
                return 0;  // No match - claim that is starts at position 0
            }
        }
        else if (ok)
        {
            return pmatch[0].rm_so;
        }
    }

    return std::string::npos;
}


bool Foam::regExpPosix::match(const std::string& text) const
{
    bool ok = false;

    if (!preg_ || !ctrl_)
    {
        // Undefined: never matches
        return false;
    }

    const auto len = text.length();

    if (len)
    {
        const size_t nmatch = 1;
        regmatch_t pmatch[1];

        // Verify that the entire string was matched
        // - [0] is the entire match result
        ok =
        (
            regexec(preg_, text.c_str(), nmatch, pmatch, 0) == 0
         && fullMatch(pmatch[0], len)
        );
    }

    return (ctrl_ == ctrlType::NEGATED ? !ok : ok);
}


bool Foam::regExpPosix::match
(
    const std::string& text,
    SubStrings<std::string>& matches
) const
{
    matches.clear();

    // Probably does not make sense for negated pattern...
    if (negated())
    {
        return match(text);
    }

    const auto len = text.size();
    if (preg_ && len)
    {
        const size_t nmatch = ngroups() + 1;
        regmatch_t pmatch[nmatch];

        // Verify that the entire string was matched
        // - [0] is the entire match result
        // - [1..] are the match groups (1..)
        if
        (
            regexec(preg_, text.c_str(), nmatch, pmatch, 0) != 0
         || !fullMatch(pmatch[0], len)
        )
        {
            return false;
        }

        matches.reserve(nmatch);

        for (size_t matchi = 0; matchi < nmatch; ++matchi)
        {
            const auto& mat = pmatch[matchi];

            if (mat.rm_so != -1 && mat.rm_eo != -1)
            {
                matches.append
                (
                    text.cbegin() + mat.rm_so,
                    text.cbegin() + mat.rm_eo
                );
            }
            else
            {
                // This may be misleading...
                matches.append(text.cbegin(), text.cbegin());
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
