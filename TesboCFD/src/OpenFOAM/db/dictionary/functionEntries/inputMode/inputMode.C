/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "inputMode.H"
#include "dictionary.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputMode,
        execute,
        dictionaryIstream,
        inputMode
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeDefault,
        execute,
        dictionaryIstream,
        default
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeMerge,
        execute,
        dictionaryIstream,
        merge
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeOverwrite,
        execute,
        dictionaryIstream,
        overwrite
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeWarn,
        execute,
        dictionaryIstream,
        warn
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeError,
        execute,
        dictionaryIstream,
        error
    );
} // End namespace functionEntries
} // End namespace Foam


const Foam::Enum
<
    Foam::entry::inputMode
>
Foam::functionEntries::inputMode::selectableNames
({
    { entry::inputMode::MERGE,  "merge" },
    { entry::inputMode::OVERWRITE, "overwrite" },
    { entry::inputMode::PROTECT, "protect" },
    { entry::inputMode::WARN, "warn" },
    { entry::inputMode::ERROR, "error" },
    // Aliases
    { entry::inputMode::MERGE, "default" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::functionEntries::inputMode::execute
(
    dictionary& unused,
    Istream& is
)
{
    const word modeName(is);

    // Like Enum::getOrDefault() with failsafe behaviour
    if (selectableNames.found(modeName))
    {
        entry::globalInputMode = selectableNames.get(modeName);
    }
    else
    {
        WarningInFunction
            << "Unsupported inputMode '" << modeName
            << "' ... defaulting to 'merge'"
            << endl;

        entry::resetInputMode();
    }

    return true;
}


bool Foam::functionEntries::inputModeDefault::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return entry::New(parentDict, is, entry::inputMode::PROTECT);
}


bool Foam::functionEntries::inputModeMerge::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return entry::New(parentDict, is, entry::inputMode::MERGE);
}


bool Foam::functionEntries::inputModeOverwrite::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return entry::New(parentDict, is, entry::inputMode::OVERWRITE);
}


bool Foam::functionEntries::inputModeWarn::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return entry::New(parentDict, is, entry::inputMode::WARN);
}


bool Foam::functionEntries::inputModeError::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    return entry::New(parentDict, is, entry::inputMode::ERROR);
}


// ************************************************************************* //
