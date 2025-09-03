/*---------------------------------------------------------------------------*\
  =========                |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "userValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace substitutionModels
{
    defineTypeNameAndDebug(userValue, 0);
    addToRunTimeSelectionTable
    (
        substitutionModel,
        userValue,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::substitutionModels::userValue::userValue
(
    const dictionary& dict,
    const Time& time
)
:
    substitutionModel(dict, time),
    entries_()
{
    // Populate entries
    const dictionary& entriesDict = dict.subDict("entries");
    for (const auto& e : entriesDict)
    {
        entries_.insert(cleanKey(e.keyword()), string(e.stream()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::substitutionModels::userValue::valid
(
    const word& keyName
) const
{
    return entries_.found(keyName);
}


bool Foam::substitutionModels::userValue::apply
(
    const word& key,
    string& buffer
) const
{
    if (!valid(key)) return false;

    buffer.replaceAll(keyify(key), entries_[key]);

    return true;
}


Foam::wordList Foam::substitutionModels::userValue::keys() const
{
    return entries_.sortedToc();
}


// ************************************************************************* //