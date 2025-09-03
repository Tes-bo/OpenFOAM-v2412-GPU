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

#include "entry.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "functionEntry.H"
#include "includeEntry.H"
#include "stringOps.H"
#include "dictionaryListEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::entry::getKeyword(keyType& keyword, token& keyToken, Istream& is)
{
    // Read the next valid token discarding spurious ';'s
    do
    {
        if
        (
            is.read(keyToken).bad()
         || is.eof()
         || !keyToken.good()
        )
        {
            return false;
        }
    }
    while (keyToken == token::END_STATEMENT);

    // If the token is a valid keyword set 'keyword' return true...
    if (keyToken.isWord())
    {
        keyword = keyToken.wordToken();
        return true;
    }

    if (keyToken.isString())
    {
        // Enable wildcards
        keyword = keyToken.stringToken();
        return true;
    }

    return false;
}


bool Foam::entry::getKeyword(keyType& keyword, Istream& is)
{
    token keyToken;
    const bool valid = getKeyword(keyword, keyToken, is);

    if (valid)
    {
        return true;
    }

    // Mark as invalid, but allow for some more checking
    if (keyToken == token::END_BLOCK || is.eof())
    {
        return false;
    }

    // Otherwise the token is invalid
    std::cerr
        << "--> FOAM Warning :" << nl
        << "    From function " << FUNCTION_NAME << nl
        << "    in file " << __FILE__ << " at line " << __LINE__ << nl
        << "    Reading " << is.relativeName() << nl
        << "    found " << keyToken << nl
        << "    expected either " << token::END_BLOCK << " or EOF"
        << std::endl;

    return false;
}


bool Foam::entry::New
(
    dictionary& parentDict,
    Istream& is,
    const entry::inputMode inpMode,
    const int endChar
)
{
    // The inputMode for dealing with duplicate entries
    const entry::inputMode mode =
    (
        inpMode == inputMode::GLOBAL
      ? globalInputMode
      : inpMode
    );

    // If somehow the global itself is 'global' - this is a severe logic error.
    if (mode == inputMode::GLOBAL)
    {
        FatalIOErrorInFunction(is)
            << "Cannot use 'GLOBAL' as an inputMode"
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);

    keyType keyword;
    token keyToken;

    // Get the next keyword and if a valid keyword return true
    const bool valid = getKeyword(keyword, keyToken, is);

    // Can accept a list of entries too
    if (keyToken.isLabel() || keyToken.isPunctuation(token::BEGIN_LIST))
    {
        is.putBack(keyToken);
        return parentDict.add
        (
            new dictionaryListEntry(parentDict, is),
            false
        );
    }

    if (!valid)
    {
        // Error processing for invalid or unexpected input

        // Do some more checking
        if (keyToken == token::END_BLOCK)
        {
            if (token::END_BLOCK != endChar)
            {
                FatalIOErrorInFunction(is)
                    << "Unexpected '}' while reading dictionary entry"
                    << exit(FatalIOError);
            }
            return false;
        }
        if (is.eof())
        {
            if (endChar)
            {
                FatalIOErrorInFunction(is)
                    << "Unexpected EOF while reading dictionary entry"
                    << exit(FatalIOError);
            }
            return false;
        }


        if (endChar)
        {
            FatalIOErrorInFunction(is)
                << "Found " << keyToken
                << " but expected " << char(endChar)
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Found " << keyToken
                << " but expected EOF, or perhaps a '}' char"
                << exit(FatalIOError);
        }

        return false;
    }


    if (keyword[0] == token::HASH)
    {
        // Function entry - #function

        if (disableFunctionEntries)
        {
            return parentDict.add
            (
                new functionEntry
                (
                    keyword,
                    parentDict,
                    is
                ),
                false
            );
        }

        const word functionName(keyword.substr(1), false);
        return functionEntry::execute(functionName, parentDict, is);
    }


    if (!disableFunctionEntries && keyword[0] == token::DOLLAR)
    {
        // Substitution entry - $variable

        token nextToken(is);
        is.putBack(nextToken);

        if (keyword.size() > 2 && keyword[1] == token::BEGIN_BLOCK)
        {
            // Recursive substitution mode.
            // Content between {} is replaced with expansion.
            // Then let standard variable expansion deal with rest.
            string expanded = keyword.substr(2, keyword.size()-3);

            // Substitute dictionary and environment variables.
            // Do not allow empty substitutions.
            stringOps::inplaceExpand(expanded, parentDict, true, false);

            // Restore the '$' prefix.
            // Use replace since operator= is private
            keyword.std::string::replace(1, keyword.size()-1, expanded);
        }

        if (nextToken == token::BEGIN_BLOCK)
        {
            const word varName = keyword.substr(1);

            // Lookup the variable name in the given dictionary
            const auto finder =
                parentDict.csearchScoped(varName, keyType::REGEX_RECURSIVE);

            if (finder.good())
            {
                // Read as primitiveEntry
                const keyType newKeyword(finder.stream());

                return parentDict.add
                (
                    new dictionaryEntry(newKeyword, parentDict, is),
                    false
                );
            }

            FatalIOErrorInFunction(is)
                << "Attempt to use undefined variable " << varName
                << " as keyword"
                << exit(FatalIOError);
            return false;
        }
        else
        {
            // Deal with duplicate entries (at least partially)
            const bool mergeEntry =
            (
                mode == inputMode::MERGE
             || mode == inputMode::OVERWRITE
            );

            parentDict.substituteScopedKeyword(keyword, mergeEntry);
        }

        return true;
    }


    // Normal or scoped entry
    {
        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::END_LIST)
        {
            FatalIOErrorInFunction(is)
                << "Unexpected token encountered for "
                << keyword << " - " << nextToken.info()
                << exit(FatalIOError);
            return false;
        }

        const bool scoped =
        (
            !disableFunctionEntries
         && keyword.contains('/')
        );

        // See (using exact match) if entry already present
        auto finder =
        (
            scoped
          ? parentDict.searchScoped(keyword, keyType::LITERAL)
          : parentDict.search(keyword, keyType::LITERAL)
        );

        // How to manage duplicate entries
        bool mergeEntry = false;

        if (finder.good())
        {
            // Use keyword from the found entry (ie, eliminate scoping chars)
            const keyType key = finder.ref().keyword();

            if (mode == inputMode::PROTECT || keyword == "FoamFile")
            {
                // Read and discard if existing element should be protected,
                // or would potentially alter the "FoamFile" header.

                // Disable function/variable expansion to avoid side-effects
                const int oldFlag = entry::disableFunctionEntries;
                entry::disableFunctionEntries = 1;

                if (nextToken == token::BEGIN_BLOCK)
                {
                    dictionaryEntry dummy("dummy", finder.context(), is);
                }
                else
                {
                    primitiveEntry  dummy("dummy", finder.context(), is);
                }

                entry::disableFunctionEntries = oldFlag;
                return true;
            }

            if (mode == inputMode::ERROR)
            {
                FatalIOErrorInFunction(is)
                    << "duplicate entry: " << key
                    << exit(FatalIOError);

                return false;
            }

            if (mode == inputMode::MERGE)
            {
                mergeEntry = true;
            }
            else if (mode == inputMode::OVERWRITE)
            {
                // Clear existing dictionary so merge acts like overwrite
                if (finder.isDict())
                {
                    finder.dict().clear();
                }
                mergeEntry = true;
            }

            // Merge/overwrite data entry

            if (nextToken == token::BEGIN_BLOCK)
            {
                return finder.context().add
                (
                    new dictionaryEntry(key, finder.context(), is),
                    mergeEntry
                );
            }
            else
            {
                return finder.context().add
                (
                    new primitiveEntry(key, finder.context(), is),
                    mergeEntry
                );
            }
        }
        else if (scoped)
        {
            // A slash-scoped entry - did not previously exist

            string fullPath(keyword);
            fileName::clean(fullPath);

            // Get or create the dictionary-path.
            // fileName::path == dictionary-path
            dictionary* subDictPtr =
                parentDict.makeScopedDict(fileName::path(fullPath));

            if (subDictPtr)
            {
                // fileName::name == keyword-name
                string keyName = fileName::name(fullPath);
                keyType key;

                // Patterns allowed for the final element.
                // - use if key name begins with a (single|double) quote

                if (keyName.find_first_of("\"'") == 0)
                {
                    // Begins with a quote - treat as pattern
                    key = keyType
                    (
                        string::validate<keyType>(keyName),
                        keyType::REGEX
                    );
                }
                else
                {
                    // Treat as a word
                    key = word::validate(keyName, false);
                }

                if (nextToken == token::BEGIN_BLOCK)
                {
                    return subDictPtr->add
                    (
                        new dictionaryEntry(key, *subDictPtr, is),
                        mergeEntry
                    );
                }
                else
                {
                    return subDictPtr->add
                    (
                        new primitiveEntry(key, *subDictPtr, is),
                        mergeEntry
                    );
                }
            }

            // Some error finding/creating intermediate dictionaries
            return false;
        }
        else
        {
            // A non-scoped entry - did not previously exist

            if (nextToken == token::BEGIN_BLOCK)
            {
                return parentDict.add
                (
                    new dictionaryEntry(keyword, parentDict, is),
                    mergeEntry
                );
            }
            else
            {
                return parentDict.add
                (
                    new primitiveEntry(keyword, parentDict, is),
                    mergeEntry
                );
            }
        }
    }
}


Foam::autoPtr<Foam::entry> Foam::entry::New(Istream& is)
{
    is.fatalCheck(FUNCTION_NAME);

    autoPtr<entry> ptr;

    // Get the next keyword and if invalid return false
    keyType keyword;
    if (getKeyword(keyword, is))
    {
        // Keyword starts entry ...
        token nextToken(is);
        is.putBack(nextToken);

        if (nextToken == token::BEGIN_BLOCK)
        {
            // A sub-dictionary
            ptr.reset(new dictionaryEntry(keyword, dictionary::null, is));
        }
        else
        {
            ptr.reset(new primitiveEntry(keyword, is));
        }
    }

    return ptr;
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const entry& e)
{
    e.write(os);
    return os;
}


// ************************************************************************* //
