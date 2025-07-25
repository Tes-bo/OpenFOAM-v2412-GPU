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

#include "setsToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setsToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, setsToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, setsToFaceZone, istream);

    addToRunTimeSelectionTable(topoSetFaceZoneSource, setsToFaceZone, word);
    addToRunTimeSelectionTable(topoSetFaceZoneSource, setsToFaceZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::setsToFaceZone::usage_
(
    setsToFaceZone::typeName,
    "\n    Usage: setsToFaceZone <faceSet> <slaveCellSet>\n\n"
    "    Select all faces in the faceSet."
    " Orientated so slave side is in cellSet.\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setsToFaceZone::setsToFaceZone
(
    const polyMesh& mesh,
    const word& faceSetName,
    const word& cellSetName,
    const bool flip
)
:
    topoSetFaceZoneSource(mesh),
    faceSetName_(faceSetName),
    cellSetName_(cellSetName),
    flip_(flip)
{}


Foam::setsToFaceZone::setsToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    faceSetName_(dict.get<word>("faceSet")),
    cellSetName_(dict.get<word>("cellSet")),
    flip_(dict.getOrDefault("flip", false))
{}


Foam::setsToFaceZone::setsToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    faceSetName_(checkIs(is)),
    cellSetName_(checkIs(is)),
    flip_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setsToFaceZone::applyToSet
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
                    << faceSetName_ << "; orientation pointing ";

                if (flip_)
                {
                    Info<< "into";
                }
                else
                {
                    Info<< "away from";
                }

                Info<< " cell set: " << cellSetName_ << " ..." << endl;
            }

            // Load the sets
            faceSet fSet(mesh_, faceSetName_, IOobject::NO_REGISTER);
            cellSet cSet(mesh_, cellSetName_, IOobject::NO_REGISTER);

            // Start off from copy
            DynamicList<label> newAddressing(zoneSet.addressing());
            DynamicList<bool> newFlipMap(zoneSet.flipMap());

            for (const label facei : fSet)
            {
                if (!zoneSet.found(facei))
                {
                    bool flipFace = false;

                    const label own = mesh_.faceOwner()[facei];
                    const bool ownFound = cSet.found(own);

                    if (mesh_.isInternalFace(facei))
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        bool neiFound = cSet.found(nei);

                        if (ownFound && !neiFound)
                        {
                            flipFace = false;
                        }
                        else if (!ownFound && neiFound)
                        {
                            flipFace = true;
                        }
                        else
                        {
                            WarningInFunction
                                << "One of owner or neighbour of internal face "
                                << facei << " should be in cellSet "
                                << cSet.name()
                                << " to be able to determine orientation."
                                << endl
                                << "Face:" << facei << " own:" << own
                                << " OwnInCellSet:" << ownFound
                                << " nei:" << nei
                                << " NeiInCellSet:" << neiFound
                                << endl;
                        }
                    }
                    else
                    {
                        flipFace = !ownFound;
                    }


                    if (flip_)
                    {
                        flipFace = !flipFace;
                    }

                    newAddressing.append(facei);
                    newFlipMap.append(flipFace);
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
                    << faceSetName_ << " ..." << endl;
            }

            // Load the set
            faceZoneSet loadedSet(mesh_, faceSetName_);

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
