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

#include "setToCellZone.H"
#include "polyMesh.H"
#include "cellZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setToCellZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setToCellZone, word);
    addToRunTimeSelectionTable(topoSetSource, setToCellZone, istream);

    addToRunTimeSelectionTable(topoSetCellZoneSource, setToCellZone, word);
    addToRunTimeSelectionTable(topoSetCellZoneSource, setToCellZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::setToCellZone::usage_
(
    setToCellZone::typeName,
    "\n    Usage: setToCellZone <cellSet>\n\n"
    "    Select all cells in the cellSet.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setToCellZone::setToCellZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetCellZoneSource(mesh),
    setName_(setName)
{}


Foam::setToCellZone::setToCellZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellZoneSource(mesh, dict),
    setName_(dict.get<word>("set"))
{}


Foam::setToCellZone::setToCellZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellZoneSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToCellZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<cellZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a cellZoneSet." << endl;
        return;
    }
    else
    {
        cellZoneSet& zoneSet = refCast<cellZoneSet>(set);

        if (action == topoSetSource::ADD || action == topoSetSource::NEW)
        {
            if (verbose_)
            {
                Info<< "    Adding all cells from cell set: "
                    << setName_ << " ..." << endl;
            }

            // Load the sets
            cellSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);

            // Start off from copy
            DynamicList<label> newAddressing(zoneSet.addressing());

            for (const label celli : loadedSet)
            {
                if (!zoneSet.found(celli))
                {
                    newAddressing.push_back(celli);
                }
            }

            zoneSet.addressing().transfer(newAddressing);
            zoneSet.updateSet();
        }
        else if (action == topoSetSource::SUBTRACT)
        {
            if (verbose_)
            {
                Info<< "    Removing all cells from cell set: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            cellSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);

            // Start off empty
            DynamicList<label> newAddressing(zoneSet.addressing().size());

            for (const label celli : zoneSet.addressing())
            {
                if (!loadedSet.found(celli))
                {
                    newAddressing.push_back(celli);
                }
            }
            zoneSet.addressing().transfer(newAddressing);
            zoneSet.updateSet();
        }
    }
}


// ************************************************************************* //
