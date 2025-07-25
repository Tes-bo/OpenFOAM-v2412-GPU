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

#include "faceZoneToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToFaceZone, istream);

    addToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        faceZoneToFaceZone,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        faceZoneToFaceZone,
        istream
    );
}


Foam::topoSetSource::addToUsageTable Foam::faceZoneToFaceZone::usage_
(
    faceZoneToFaceZone::typeName,
    "\n    Usage: faceZoneToFaceZone <faceZone>\n\n"
    "    Select all faces in the faceZone\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetFaceZoneSource(mesh),
    setName_(setName)
{}


Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    setName_(dict.get<word>("zone"))
{}


Foam::faceZoneToFaceZone::faceZoneToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    setName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneToFaceZone::applyToSet
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
                Info<< "    Adding all faces from face zone: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            faceZoneSet loadedSet(mesh_, setName_);

            DynamicList<label> newAddressing(zoneSet.addressing());
            DynamicList<bool> newFlipMap(zoneSet.flipMap());

            forAll(loadedSet.addressing(), i)
            {
                if (!zoneSet.found(loadedSet.addressing()[i]))
                {
                    newAddressing.append(loadedSet.addressing()[i]);
                    newFlipMap.append(loadedSet.flipMap()[i]);
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
                Info<< "    Removing all faces from face zone: "
                    << setName_ << " ..." << endl;
            }

            // Load the set
            faceZoneSet loadedSet(mesh_, setName_);

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
