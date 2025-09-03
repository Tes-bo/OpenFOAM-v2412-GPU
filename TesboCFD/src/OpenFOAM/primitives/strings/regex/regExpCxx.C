/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "regExpCxx.H"
#include "debug.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::regExpCxx::grammar(Foam::debug::optimisationSwitch("regExpCxx", 0));


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static std::string error_string(const std::regex_error& err)
{
    switch (err.code())
    {
        case std::regex_constants::error_collate :
            return "invalid collating element name";
            break;

        case std::regex_constants::error_ctype :
            return "invalid character class name";
            break;

        case std::regex_constants::error_escape :
            return "invalid escaped character or a trailing escape";
            break;

        case std::regex_constants::error_backref :
            return "invalid back reference";
            break;

        case std::regex_constants::error_brack :
            return "mismatched [ and ]";
            break;

        case std::regex_constants::error_paren :
            return "mismatched ( and )";
            break;

        case std::regex_constants::error_brace :
            return "mismatched { and }";
            break;

        case std::regex_constants::error_badbrace :
            return "invalid range in a {..}";
            break;

        case std::regex_constants::error_range :
            return "invalid [..] character range";
            break;

        case std::regex_constants::error_space :
            return "memory error";
            break;

        case std::regex_constants::error_badrepeat :
            return "bad '*?+{' repeat";
            break;

        case std::regex_constants::error_complexity :
            return "expression too complex";
            break;

        case std::regex_constants::error_stack :
            return "memory stack error";
            break;

        default:
            break;

    }

    return "";
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::regExpCxx::set_pattern
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
        std::regex::flag_type flags = syntax();
        if (ignoreCase)
        {
            flags |= std::regex::icase;
        }

        try
        {
            re_.assign(pat, len, flags);
            ctrl_ = (doNegate ? ctrlType::NEGATED : ctrlType::NORMAL);
            return true;
        }
        catch (const std::regex_error& err)
        {
            FatalErrorInFunction
                << "Failed to compile regular expression '"
                << pattern << "'" << nl
                << err.what() << ": " << error_string(err).c_str() << nl
                << exit(FatalError);
        }
    }

    return false;
}


// ************************************************************************* //
