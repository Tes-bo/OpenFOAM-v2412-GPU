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

#include "dictionaryValue.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace substitutionModels
{
    defineTypeNameAndDebug(dictionaryValue, 0);
    addToRunTimeSelectionTable(substitutionModel, dictionaryValue, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::substitutionModels::dictionaryValue::processDict
(
    const dictionary& dict,
    const word& key,
    string& buffer
) const
{
    const string& lookup = entries_[key];

    OStringStream oss;
    if (lookup.empty())
    {
        // Add complete dictionary
        oss << dict;
    }
    else
    {
        const entry* ePtr = dict.findScoped(lookup);

        if (!ePtr)
        {
            WarningInFunction
                << "Unable to find entry " << lookup
                << endl;
            return false;
        }

        if (ePtr->isDict())
        {
            const dictionary& de = ePtr->dict();

            // Write dictionary contents
            oss << de.dictName() << de;
        }
        else
        {
            for (const auto& t : ePtr->stream())
            {
                if (oss.count()) oss << separator_;
                oss << t;
            }
        }
    }

    buffer.replaceAll(keyify(key), oss.str());

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::substitutionModels::dictionaryValue::dictionaryValue
(
    const dictionary& dict,
    const Time& time
)
:
    substitutionModel(dict, time),
    object_(),
    region_(polyMesh::defaultRegion),
    path_(),
    separator_(dict.getOrDefault<word>("separator", " ")),
    entries_()
{
    const auto* oPtr = dict.findEntry("object");
    const auto* pPtr = dict.findEntry("path");

    if (oPtr && pPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Specify either 'object' or 'path' but not both"
            << exit(FatalIOError);
    }

    if (oPtr)
    {
        // Optionally read the region
        dict.readIfPresent<word>("region", region_);

        // Must read the object name to look up
        object_ = dict.get<word>("object");
    }

    if (pPtr)
    {
        path_ = dict.get<fileName>("path").expand();
    }

    // Populate entries
    const dictionary& entriesDict = dict.subDict("entries");
    for (const auto& e : entriesDict)
    {
        entries_.insert(cleanKey(e.keyword()), string(e.stream()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::substitutionModels::dictionaryValue::valid(const word& keyName) const
{
    return entries_.found(keyName);
}


bool Foam::substitutionModels::dictionaryValue::apply
(
    const word& key,
    string& buffer
) const
{
    if (!valid(key)) return false;

    if (path_.size())
    {
        fileName path(path_);
        if (replaceBuiltin(path))
        {
            path.clean();
        }

        IFstream is(path);

        if (!is.good())
        {
            WarningInFunction
                << "Unable to find dictionary at " << path
                << ". Deactivating." << endl;

            return false;
        }

        return processDict(dictionary(is), key, buffer);
    }
    else
    {
        const auto* obrPtr = time_.cfindObject<objectRegistry>(region_);

        if (!obrPtr)
        {
            WarningInFunction
                << "Unable to find region " << region_
                << ". Deactivating." << endl;

            return false;
        }

        // Find object; recursive lookup into parent
        const auto* dictPtr = obrPtr->cfindObject<IOdictionary>(object_, true);

        if (!dictPtr)
        {
            WarningInFunction
                << "Unable find dictionary " << object_
                << " on region " << region_
                << ". Deactivating." << endl;

            return false;
        }

        return processDict(*dictPtr, key, buffer);
    }
}


Foam::wordList Foam::substitutionModels::dictionaryValue::keys() const
{
    return entries_.sortedToc();
}


// ************************************************************************* //