/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "setToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, setToFaceZone, istream);

    addToRunTimeSelectionTable(topoSetFaceZoneSource, setToFaceZone, word);
    addToRunTimeSelectionTable(topoSetFaceZoneSource, setToFaceZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::setToFaceZone::usage_
(
    setToFaceZone::typeName,
    "\n    Usage: setToFaceZone <faceSet>\n\n"
    "    Select all faces in the faceSet."
    " Sets flipMap.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setToFaceZone::setToFaceZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetFaceZoneSource(mesh),
    setName_(setName)
{}


Foam::setToFaceZone::setToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    setName_(dict.get<word>("faceSet"))
{
    if (dict.found("cellSet"))
    {
        IOWarningInFunction(dict)
            << "Ignoring entry 'cellSet' - maybe use setsToFaceZone instead ?"
            << endl;
    }
}


Foam::setToFaceZone::setToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setToFaceZone::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!isA<faceZoneSet>(set))
    {
        WarningInFunction
            << "Operation only allowed on a faceZoneSet." << endl;
        return;
    }
    else
    {
        faceZoneSet& zoneSet = refCast<faceZoneSet>(set);

        if (action == topoSetSource::ADD || action == topoSetSource::NEW)
        {
            if (verbose_)
            {
                Info<< "    Adding all faces from face set: "
                    << setName_ << " ..." << endl;
            }

            // Load the sets
            faceSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);
            const labelHashSet& faceLabels = loadedSet;

            // Start off from copy
            DynamicList<label> newAddressing(zoneSet.addressing());
            DynamicList<bool> newFlipMap(zoneSet.flipMap());

            for (const label facei : faceLabels)
            {
                if (!zoneSet.found(facei))
                {
                    newAddressing.append(facei);
                    newFlipMap.append(false);
                }
            }

            zoneSet.addressing().transfer(newAddressing);
            zoneSet.flipMap().transfer(newFlipMap);
            zoneSet.updateSet();
        }
        else if (action == topoSetSource::SUBTRACT)
        {
            if (verbose_)
            {
                Info<< "    Removing all faces from face set: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            faceSet loadedSet(mesh_, setName_, IOobject::NO_REGISTER);

            // Start off empty
            DynamicList<label> newAddressing(zoneSet.addressing().size());
            DynamicList<bool> newFlipMap(zoneSet.flipMap().size());

            forAll(zoneSet.addressing(), i)
            {
                if (!loadedSet.found(zoneSet.addressing()[i]))
                {
                    newAddressing.append(zoneSet.addressing()[i]);
                    newFlipMap.append(zoneSet.flipMap()[i]);
                }
            }
            zoneSet.addressing().transfer(newAddressing);
            zoneSet.flipMap().transfer(newFlipMap);
            zoneSet.updateSet();
        }
    }
}


// ************************************************************************* //
