/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "wordDirective.H"
#include "dictionary.H"
#include "stringOps.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        wordDirective,
        execute,
        dictionaryIstream,
        word
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        wordDirective,
        execute,
        primitiveEntryIstream,
        word
    );

} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::token Foam::functionEntries::wordDirective::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    token tok(is);

    // The string to evaluate
    string str;

    if (tok.isStringType())  // Also accepts a single bare word
    {
        // - #word expr
        // - #word "expr"
        // - #word #{ expr #}
        str = tok.stringToken();
    }
    else if (tok.isPunctuation(token::BEGIN_BLOCK))
    {
        // - #word { expr }
        // strip comments
        if (!continueReadUntilRightBrace(is, str, true))
        {
            reportReadWarning
            (
                is,
                "Premature end while reading #word - missing '}'?"
            );
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Invalid input for #word."
               " Expecting a string or block to expand, but found" << nl
            << tok.info() << endl
            << exit(FatalIOError);
    }

    stringOps::inplaceExpand(str, parentDict);

    word result(word::validate(str));  // Includes trimming etc.

    if (!result.empty())
    {
        tok = std::move(result);
        return tok;
    }

    // Expanded to nothing - treat as a no-op
    return token::undefinedToken;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::wordDirective::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    token tok(evaluate(parentDict, is));

    if (tok.good())
    {
        // Can add evaluated value directly into primitiveEntry
        entry.add_tokens(std::move(tok));   // Add at tokenIndex
    }

    return true;
}


bool Foam::functionEntries::wordDirective::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    token tok(evaluate(parentDict, is));

    if (tok.good())
    {
        // Using putBack to insert evaluated value into stream
        is.putBack(tok);
    }

    return true;
}


// ************************************************************************* //
