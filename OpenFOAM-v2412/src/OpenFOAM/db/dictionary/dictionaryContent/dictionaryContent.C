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

#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class UnaryPredicate>
static dictionary copyFilteredDict
(
    const dictionary& input,
    const UnaryPredicate& pred
)
{
    dictionary dict;
    dict.name() = input.name();  // rename

    for (const entry& e : input)
    {
        const keyType& key = e.keyword();

        bool accept = false;

        if (key.isPattern())
        {
            // Non-trivial to filter a regex itself - so just accept it
            // - could also have a "pruneRegex" flag (for example)
            accept = true;
        }
        else
        {
            accept = pred(key);
        }

        if (accept)
        {
            dict.add(e, false);  // No merge - entries are already unique
        }
    }

    return dict;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::dictionary
Foam::dictionaryContent::copyDict
(
    const dictionary& input,
    const wordList& allow,
    const wordList& deny
)
{
    if (allow.empty())
    {
        if (deny.empty())
        {
            // Pass-through
            return dictionary(input);
        }

        // Deny only
        return copyFilteredDict
        (
            input,
            [&](const std::string& key) { return !deny.found(key); }
        );
    }

    // Allow is non-empty

    // No general case with deny as well
    return copyFilteredDict
    (
        input,
        [&](const std::string& key) { return allow.found(key); }
    );
}


Foam::dictionary
Foam::dictionaryContent::copyDict
(
    const dictionary& input,
    const wordRes& allow,
    const wordRes& deny
)
{
    if (allow.empty())
    {
        if (deny.empty())
        {
            // Pass-through
            return dictionary(input);
        }

        // Deny only
        return copyFilteredDict
        (
            input,
            [&](const std::string& key) { return !deny.match(key); }
        );
    }

    // Allow is non-empty

    if (deny.empty())
    {
        return copyFilteredDict
        (
            input,
            [&](const std::string& key) { return allow.match(key); }
        );
    }

    // General case - have both deny and allow
    return copyFilteredDict
    (
        input,
        [&](const std::string& key)
        {
            const auto result = allow.matched(key);
            return
            (
                result == wordRe::LITERAL
              ? true
              : (result == wordRe::REGEX && !deny.match(key))
            );
        }
    );
}


// ************************************************************************* //
