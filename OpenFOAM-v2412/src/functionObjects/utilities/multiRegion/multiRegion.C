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

#include "multiRegion.H"
#include "fvMesh.H"
#include "timeControlFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(multiRegion, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        multiRegion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::multiRegion::multiRegion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::multiRegion::read(const dictionary& dict)
{
    if (functionObject::read(dict))
    {
        Info<< type() << ' ' << name() << ':' << nl;

        wordRes regionMatcher;

        const bool matchAny =
            !dict.readIfPresent<wordRes>("regions", regionMatcher);

        const dictionary& functionDict = dict.subDict("function");

        const wordList allRegions(time_.sortedNames<fvMesh>());
        functions_.resize(allRegions.size());

        const bool needsTimeControl = timeControl::entriesPresent(functionDict);

        label functioni = 0;
        for (const word& regionName : allRegions)
        {
            if (matchAny || regionMatcher.match(regionName))
            {
                const word localName(IOobject::scopedName(name(), regionName));

                dictionary regionDict(functionDict);
                regionDict.add("region", regionName);

                if (needsTimeControl)
                {
                    functions_.set
                    (
                        functioni,
                        new timeControl(localName, time_, regionDict)
                    );
                }
                else
                {
                    functions_.set
                    (
                        functioni,
                        functionObject::New(localName, time_, regionDict)
                    );
                }

                ++functioni;
            }
        }

        functions_.resize(functioni);

        if (functions_.empty())
        {
            WarningInFunction
                << "No regions applied"
                << endl;

            return false;
        }

        Info<< "    Spawned additional object(s):" << nl;
        for (const auto& f : functions_)
        {
            Info<< "        " << f.name() << nl;
        }

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::multiRegion::execute()
{
    bool result = true;

    for (auto& f : functions_)
    {
        result = f.execute() && result;
    }

    return result;
}


bool Foam::functionObjects::multiRegion::write()
{
    bool result = true;

    for (auto& f : functions_)
    {
        result = f.write() && result;
    }

    return result;
}


// ************************************************************************* //
