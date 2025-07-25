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

#include "functionObjectValue.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace substitutionModels
{
    defineTypeNameAndDebug(functionObjectValue, 0);
    addToRunTimeSelectionTable
    (
        substitutionModel,
        functionObjectValue,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::substitutionModels::functionObjectValue::getValue
(
    OStringStream& oss,
    const word& lookup
) const
{
    const auto& foProps = time_.functionObjects().propsDict();

    Type result;
    if (foProps.getObjectResult(functionObject_, lookup, result))
    {
        oss << result;
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::substitutionModels::functionObjectValue::functionObjectValue
(
    const dictionary& dict,
    const Time& time
)
:
    substitutionModel(dict, time),
    functionObject_(dict.get<word>("functionObject")),
    entries_(),
    debugValues_(dict.getOrDefault<bool>("debugValues", false))
{
    // Populate entries
    const dictionary& entriesDict = dict.subDict("entries");
    for (const auto& e : entriesDict)
    {
        entries_.insert(cleanKey(e.keyword()), word(e.stream()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::substitutionModels::functionObjectValue::update()
{
    if (debugValues_)
    {
        Info<< nl << "Function object results:" << nl;
        time_.functionObjects().propsDict().writeAllResultEntries(Info);
    }

    return true;
}


bool Foam::substitutionModels::functionObjectValue::valid
(
    const word& keyName
) const
{
    return entries_.found(keyName);
}


bool Foam::substitutionModels::functionObjectValue::apply
(
    const word& key,
    string& buffer
) const
{
    if (!valid(key)) return false;

    OStringStream oss;

    const word& lookup = entries_[key];

    bool ok =
        getValue<label>(oss, lookup)
     || getValue<scalar>(oss, lookup)
     || getValue<vector>(oss, lookup)
     || getValue<sphericalTensor>(oss, lookup)
     || getValue<symmTensor>(oss, lookup)
     || getValue<tensor>(oss, lookup);

    if (!ok) return false;

    buffer.replaceAll(keyify(key), oss.str());

    return true;
}


Foam::wordList Foam::substitutionModels::functionObjectValue::keys() const
{
    return entries_.sortedToc();
}


// ************************************************************************* //