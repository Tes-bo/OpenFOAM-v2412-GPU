/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "setAndNormalToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setAndNormalToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setAndNormalToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, setAndNormalToFaceZone, istream);

    addToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        setAndNormalToFaceZone,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetFaceZoneSource,
        setAndNormalToFaceZone,
        istream
    );
}


Foam::topoSetSource::addToUsageTable Foam::setAndNormalToFaceZone::usage_
(
    setAndNormalToFaceZone::typeName,
    "\n    Usage: setAndNormalToFaceZone <faceSet> <normal>\n\n"
    "    Select all faces in the faceSet and orient using normal.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setAndNormalToFaceZone::setAndNormalToFaceZone
(
    const polyMesh& mesh,
    const word& setName,
    const vector& normal
)
:
    topoSetFaceZoneSource(mesh),
    setName_(setName),
    normal_(normal)
{}


Foam::setAndNormalToFaceZone::setAndNormalToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    setName_(dict.get<word>("faceSet")),
    normal_(dict.get<vector>("normal"))
{}


Foam::setAndNormalToFaceZone::setAndNormalToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    setName_(checkIs(is)),
    normal_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setAndNormalToFaceZone::applyToSet
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

            const faceList& faces = mesh_.faces();
            const pointField& points = mesh_.points();

            for (const label facei : faceLabels)
            {
                if (!zoneSet.found(facei))
                {
                    newAddressing.append(facei);

                    const vector n = faces[facei].areaNormal(points);
                    if ((n & normal_) > 0)
                    {
                        newFlipMap.append(false);
                    }
                    else
                    {
                        newFlipMap.append(true);
                    }
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
