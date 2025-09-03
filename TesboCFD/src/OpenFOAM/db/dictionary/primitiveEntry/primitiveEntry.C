/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "primitiveEntry.H"
#include "dictionary.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//  Find the type/position of the ":-" or ":+" alternative values
//  Returns 0, '-', '+' corresponding to not-found or ':-' or ':+'
static inline char findParameterAlternative
(
    const std::string& s,
    std::string::size_type& pos,
    std::string::size_type endPos = std::string::npos
)
{
    while (pos != std::string::npos)
    {
        pos = s.find(':', pos);
        if (pos != std::string::npos)
        {
            if (pos < endPos)
            {
                // in-range: check for '+' or '-' following the ':'
                const char altType = s[pos+1];
                if (altType == '+' || altType == '-')
                {
                    return altType;
                }

                ++pos;    // unknown/unsupported - continue at next position
            }
            else
            {
                // out-of-range: abort
                pos = std::string::npos;
            }
        }
    }

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveEntry::expandVariable
(
    const string& varName,
    const dictionary& dict
)
{
    char altType = 0; // Type ('-' or '+') for ":-" or ":+" alternatives
    word expanded;
    string altValue;

    // Any ${{ expr }} entries have been trapped and processed elsewhere

    if (varName[0] == token::BEGIN_BLOCK && varName.size() > 1)
    {
        // Replace content between {} with string expansion and
        // handle ${parameter:-word} or ${parameter:+word}

        // Copy into a word without stripping
        expanded.assign(varName, 1, varName.size()-2);

        // Substitute dictionary and environment variables.
        // - Allow environment.
        // - No empty substitutions.
        // - No sub-dictionary lookups

        stringOps::inplaceExpand(expanded, dict, true, false, false);

        // Position of ":-" or ":+" alternative values
        std::string::size_type altPos = 0;

        // Check for parameter:-word or parameter:+word
        altType = findParameterAlternative(expanded, altPos);

        if (altType)
        {
            altValue = expanded.substr(altPos + 2);
            expanded.erase(altPos);
        }

        // Catch really bad expansions and let them die soon after.
        // Eg, ${:-other} should not be allowed.
        if (expanded.empty())
        {
            altType = 0;
            altValue.clear();
        }

        // Fallthrough for further processing
    }


    // Lookup variable name in the given dictionary WITHOUT pattern matching.
    // Having a pattern match means that in this example:
    // {
    //     internalField XXX;
    //     boundaryField { ".*" {YYY;} movingWall {value $internalField;}
    // }
    // The $internalField would be matched by the ".*" !!!

    // Recursive, non-patterns

    const word& lookupName = (expanded.empty() ? varName : expanded);

    const entry* eptr =
        dict.findScoped(lookupName, keyType::LITERAL_RECURSIVE);

    if (!eptr)
    {
        // Not found - revert to environment variable
        // and parse into a series of tokens.

        // We wish to fail if the environment variable returns
        // an empty string and there is no alternative given.
        //
        // Always allow empty strings as alternative parameters,
        // since the user provided them for a reason.

        string str(Foam::getEnv(lookupName));

        if (str.empty() ? (altType == '-') : (altType == '+'))
        {
            // Not found or empty:  use ":-" alternative value
            // Found and not empty: use ":+" alternative value
            str = std::move(altValue);
        }
        else if (str.empty())
        {
            FatalIOErrorInFunction(dict)
                << "Illegal dictionary entry or environment variable name "
                << lookupName << nl
                << "Known dictionary entries: " << dict.toc() << nl
                << exit(FatalIOError);

            return false;
        }

        // Parse string into a series of tokens

        tokenList toks(ITstream::parse(str));   // ASCII

        ITstream::add_tokens(std::move(toks));  // Add at tokenIndex
    }
    else if (eptr->isDict())
    {
        // Found dictionary entry

        tokenList toks(eptr->dict().tokens());

        if (toks.empty() ? (altType == '-') : (altType == '+'))
        {
            // Not found or empty:  use ":-" alternative value
            // Found and not empty: use ":+" alternative value

            toks = ITstream::parse(altValue);   // ASCII
        }

        ITstream::add_tokens(std::move(toks));  // Add at tokenIndex
    }
    else
    {
        // Found primitive entry - copy tokens

        if (eptr->stream().empty() ? (altType == '-') : (altType == '+'))
        {
            // Not found or empty:  use ":-" alternative value
            // Found and not empty: use ":+" alternative value

            tokenList toks(ITstream::parse(altValue));  // ASCII

            ITstream::add_tokens(std::move(toks));  // Add at tokenIndex
        }
        else
        {
            ITstream::add_tokens(eptr->stream());   // Add at tokenIndex
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveEntry::primitiveEntry(const keyType& key)
:
    entry(key),
    ITstream(Foam::zero{}, key)
{}


Foam::primitiveEntry::primitiveEntry(const keyType& key, const token& tok)
:
    entry(key),
    ITstream(key, tokenList(Foam::one{}, tok))
{}


Foam::primitiveEntry::primitiveEntry(const keyType& key, token&& tok)
:
    entry(key),
    ITstream(key, tokenList(Foam::one{}, std::move(tok)))
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const UList<token>& tokens
)
:
    entry(key),
    ITstream(key, tokens)
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    List<token>&& tokens
)
:
    entry(key),
    ITstream(key, std::move(tokens))
{}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const ITstream& is
)
:
    entry(key),
    ITstream(is)
{
    ITstream::name() = fileName::concat(ITstream::name(), key, '/');
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ITstream* Foam::primitiveEntry::streamPtr() const noexcept
{
    ITstream* ptr = const_cast<primitiveEntry*>(this);
    ptr->seek(0);
    return ptr;
}


Foam::ITstream& Foam::primitiveEntry::stream() const
{
    ITstream& is = const_cast<primitiveEntry&>(*this);
    is.seek(0);
    return is;
}


const Foam::dictionary& Foam::primitiveEntry::dict() const
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return dictionary::null;
}


Foam::dictionary& Foam::primitiveEntry::dict()
{
    FatalErrorInFunction
        << "Attempt to return primitive entry " << info()
        << " as a sub-dictionary"
        << abort(FatalError);

    return const_cast<dictionary&>(dictionary::null);
}


// ************************************************************************* //
