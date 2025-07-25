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

#include "fileRegEx.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace substitutionModels
{
    defineTypeNameAndDebug(fileRegEx, 0);
    addToRunTimeSelectionTable(substitutionModel, fileRegEx, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::substitutionModels::fileRegEx::fileRegEx
(
    const dictionary& dict,
    const Time& time
)
:
    substitutionModel(dict, time),
    path_(dict.get<fileName>("path")),
    entries_(),
    sectionSeparator_
    (
        dict.getOrDefault<string>
        (
            "sectionSeparator",
            "Time ="
        )
    ),
    matchSeparator_(dict.getOrDefault<string>("matchSeparator", " ")),
    lastMatch_(dict.getOrDefault<bool>("lastMatch", true))
{
    // Populate entries
    const dictionary& entriesDict = dict.subDict("entries");
    for (const auto& e : entriesDict)
    {
        entries_.insert(cleanKey(e.keyword()), string(e.stream()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::substitutionModels::fileRegEx::valid(const word& keyName) const
{
    return entries_.found(keyName);
}


bool Foam::substitutionModels::fileRegEx::apply
(
    const word& key,
    string& buffer
) const
{
    if (!valid(key)) return false;

    fileName path(path_);
    replaceBuiltin(path);
    IFstream is(path);

    if (!is.good())
    {
        WarningInFunction
            << "Unable to find file at " << path_
            << ". Deactivating." << endl;

        return false;
    }

    Info<< "Scanning for sections beginning with "
        << sectionSeparator_ << endl;

    // For log files containing multiple time steps
    // - put strings for last time step into a string list
    DynamicList<string> lines(96);
    string line;
    bool started = sectionSeparator_.empty() ? true : false;
    while (is.good())
    {
        is.getLine(line);
        if (line.starts_with(sectionSeparator_))
        {
            started = true;
            lines.clear();
        }
        if (started)
        {
            lines.append(line);
        }
    }

    Info<< "Cached " << lines.size() << " lines" << endl;

    OStringStream oss;
    regExp re(entries_[key].c_str());

    for (const string& data : lines)
    {
        regExp::results_type match;
        if (re.match(data, match))
        {
            oss.reset();

            for (size_t i = 1; i < match.size(); ++i)
            {
                if (i > 1) oss << matchSeparator_;
                oss << match[i].str().c_str();
            }

            if (!lastMatch_) break;
        }
    }

    if (oss.count())
    {
        buffer.replaceAll(keyify(key), oss.str());
        return true;
    }

    return false;
}


Foam::wordList Foam::substitutionModels::fileRegEx::keys() const
{
    return entries_.sortedToc();
}


// ************************************************************************* //