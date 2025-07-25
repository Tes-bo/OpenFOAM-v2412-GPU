/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "setToPointZone.H"
#include "polyMesh.H"
#include "pointZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setToPointZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setToPointZone, word);
    addToRunTimeSelectionTable(topoSetSource, setToPointZone, istream);

    addToRunTimeSelectionTable(topoSetPointZoneSource, setToPointZone, word);
    addToRunTimeSelectionTable(topoSetPointZoneSource, setToPointZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::setToPointZone::usage_
(
    setToPointZone::typeName,
    "\n    Usage: setToPointZone <pointSet>\n\n"
    "    Select all points in the pointSet.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetPointZoneSource(mesh),
    setName_(setName)
{}


Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointZoneSource(mesh, dict),
    setName_(dict.get<word>("set"))
{}


Foam::setToPointZone::setToPointZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointZoneSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToPointZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<pointZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a pointZoneSet." << endl;
        return;
    }
    else
    {
        pointZoneSet& zoneSet = refCast<pointZoneSet>(set);

        if (action == topoSetSource::ADD || action == topoSetSource::NEW)
        {
            if (verbose_)
            {
                Info<< "    Adding all points from point set: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            pointSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);
            const labelHashSet& pointLabels = loadedSet;

            // Start off from copy
            DynamicList<label> newAddressing(zoneSet.addressing());

            for (const label pointi : pointLabels)
            {
                if (!zoneSet.found(pointi))
                {
                    newAddressing.push_back(pointi);
                }
            }

            zoneSet.addressing().transfer(newAddressing);
            zoneSet.updateSet();
        }
        else if (action == topoSetSource::SUBTRACT)
        {
            if (verbose_)
            {
                Info<< "    Removing all points from point set: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            pointSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);

            // Start off empty
            DynamicList<label> newAddressing(zoneSet.addressing().size());

            for (const label pointi : zoneSet.addressing())
            {
                if (!loadedSet.found(pointi))
                {
                    newAddressing.push_back(pointi);
                }
            }
            zoneSet.addressing().transfer(newAddressing);
            zoneSet.updateSet();
        }
    }
}


// ************************************************************************* //
