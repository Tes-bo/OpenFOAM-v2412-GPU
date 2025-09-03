/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "fileFormats.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Extract and merge 'default' + formatName from list of dictionaries
//
// \returns dictionary of merged options
static dictionary combineFormatOptions
(
    const word& formatName,
    std::initializer_list<const dictionary*> dicts
)
{
    dictionary options;

    // Default specification. Merge from all levels
    // - literal search only
    for (const dictionary* dict : dicts)
    {
        if
        (
            dict
         && (dict = dict->findDict("default", keyType::LITERAL)) != nullptr
        )
        {
            options.merge(*dict);
        }
    }

    // Format specification. Merge from all levels
    // - allow REGEX search
    if (!formatName.empty())
    {
        for (const dictionary* dict : dicts)
        {
            if
            (
                dict
             && (dict = dict->findDict(formatName)) != nullptr
            )
            {
                options.merge(*dict);
            }
        }
    }

    return options;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::fileFormats::getFormatOptions
(
    const dictionary& dict,
    const word& formatName,
    const word& entryName
)
{
    return combineFormatOptions
    (
        formatName,
        {
            dict.findDict(entryName, keyType::LITERAL)
        }
    );
}


Foam::dictionary Foam::fileFormats::getFormatOptions
(
    const dictionary& dict,
    const dictionary& altDict,
    const word& formatName,
    const word& entryName
)
{
    return combineFormatOptions
    (
        formatName,
        {
            dict.findDict(entryName, keyType::LITERAL),
            altDict.findDict(entryName, keyType::LITERAL)
        }
    );
}


// ************************************************************************* //
