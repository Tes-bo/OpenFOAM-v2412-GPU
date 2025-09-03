/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "messageDirective.H"
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
        messageDirective,
        execute,
        dictionaryIstream,
        message
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        messageDirective,
        execute,
        primitiveEntryIstream,
        message
    );

} // End namespace functionEntry
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::functionEntries::messageDirective::evaluate
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
        // - #message expr
        // - #message "expr"
        // - #message #{ expr #}
        str = tok.stringToken();
    }
    else if (tok.isPunctuation(token::BEGIN_BLOCK))
    {
        // - #message { expr }
        // strip comments
        if (!continueReadUntilRightBrace(is, str, true))
        {
            reportReadWarning
            (
                is,
                "Premature end while reading #message - missing '}'?"
            );
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Invalid input for #message."
               " Expecting a string or block to expand, but found" << nl
            << tok.info() << endl
            << exit(FatalIOError);
    }

    stringOps::inplaceExpand(str, parentDict);
    stringOps::inplaceTrim(str);

    if (!str.empty() && error::master())
    {
        // Use stderr directly, in case message should be part of startup
        std::cerr
            << str << " (file: \""
            << parentDict.relativeName() << "\" line: "
            << tok.lineNumber() << ")\n" << std::flush;
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::messageDirective::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    evaluate(parentDict, is);

    return true;
}


bool Foam::functionEntries::messageDirective::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    evaluate(parentDict, is);

    return true;
}

// ************************************************************************* //
