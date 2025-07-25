/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "cellToFaceZone.H"
#include "polyMesh.H"
#include "faceZoneSet.H"
#include "cellSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToFaceZone, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToFaceZone, word);
    addToRunTimeSelectionTable(topoSetSource, cellToFaceZone, istream);

    addToRunTimeSelectionTable(topoSetFaceZoneSource, cellToFaceZone, word);
    addToRunTimeSelectionTable(topoSetFaceZoneSource, cellToFaceZone, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cellToFaceZone::usage_
(
    cellToFaceZone::typeName,
    "\n    Usage: cellToFaceZone <slaveCellSet>\n\n"
    "    Select all outside faces in the cellSet."
    " Orientated so slave side is in cellSet.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFaceZone::selectFaces
(
    const bitSet& whichCells,
    bitSet& selectedFace,
    bitSet& doFlip
) const
{
    selectedFace.resize_nocopy(mesh_.nFaces());
    selectedFace = false;

    doFlip.resize_nocopy(mesh_.nFaces());
    doFlip = false;


    // Add all faces whose both neighbours are in set.

    const label nInt = mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Check all internal faces
    for (label facei = 0; facei < nInt; ++facei)
    {
        const bool ownFound = whichCells.test(own[facei]);
        const bool neiFound = whichCells.test(nei[facei]);

        if (ownFound && !neiFound)
        {
            selectedFace.set(facei);
            doFlip.set(facei, flip_);
        }
        else if (!ownFound && neiFound)
        {
            selectedFace.set(facei);
            doFlip.set(facei, !flip_);
        }
    }

    // Get coupled cell status
    boolList neiInSet(mesh_.nBoundaryFaces(), false);

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                neiInSet[facei-nInt] = whichCells.test(own[facei]);
                ++facei;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiInSet);


    // Check all boundary faces
    for (const polyPatch& pp : patches)
    {
        label facei = pp.start();
        forAll(pp, i)
        {
            const bool ownFound = whichCells.test(own[facei]);
            const bool neiFound = neiInSet[facei-nInt];

            if (ownFound && !neiFound)
            {
                selectedFace.set(facei);
                doFlip.set(facei, flip_);
            }
            else if (!ownFound && neiFound)
            {
                selectedFace.set(facei);
                doFlip.set(facei, !flip_);
            }
            ++facei;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFaceZone::cellToFaceZone
(
    const polyMesh& mesh,
    const word& setName,
    const bool flip
)
:
    topoSetFaceZoneSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false),
    flip_(flip)
{}


Foam::cellToFaceZone::cellToFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceZoneSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_)),
    flip_(dict.getOrDefault("flip", false))
{}


Foam::cellToFaceZone::cellToFaceZone
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceZoneSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false),
    flip_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFaceZone::applyToSet
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
                Info<< "    Adding all faces on outside of cell "
                    << (isZone_ ? "zones:" : "sets: ")
                    << flatOutput(names_) << "; orientation pointing ";

                if (flip_)
                {
                    Info<< "into cell sets" << endl;
                }
                else
                {
                    Info<< "away from cell sets" << endl;
                }
            }

            bitSet selectedFace(mesh_.nFaces());
            bitSet doFlip(mesh_.nFaces());

            for (const word& setName : names_)
            {
                bitSet whichCells(mesh_.nCells());
                if (isZone_)
                {
                    whichCells.set(mesh_.cellZones()[setName]);
                }
                else
                {
                    // Load the sets
                    cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
                    whichCells.setMany(loadedSet.begin(), loadedSet.end());
                }

                // Select outside faces
                selectFaces(whichCells, selectedFace, doFlip);
            }

            // Start off from copy
            DynamicList<label> newAddressing(zoneSet.addressing());
            DynamicList<bool> newFlipMap(zoneSet.flipMap());

            for (const label facei : selectedFace)
            {
                if (!zoneSet.found(facei))
                {
                    newAddressing.append(facei);
                    newFlipMap.append(doFlip[facei]);
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
                Info<< "    Removing all faces on outside of cell "
                    << (isZone_ ? "zones:" : "sets: ")
                    << flatOutput(names_) << " ..." << endl;
            }

            bitSet selectedFace(mesh_.nFaces());
            bitSet doFlip(mesh_.nFaces());
            for (const word& setName : names_)
            {
                bitSet whichCells(mesh_.nCells());
                if (isZone_)
                {
                    whichCells.set(mesh_.cellZones()[setName]);
                }
                else
                {
                    // Load the sets
                    cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
                    whichCells.setMany(loadedSet.begin(), loadedSet.end());
                }
                // Select outside faces
                selectFaces(whichCells, selectedFace, doFlip);
            }

            // Start off empty
            DynamicList<label> newAddressing(zoneSet.addressing().size());
            DynamicList<bool> newFlipMap(zoneSet.flipMap().size());

            for (const label facei : selectedFace)
            {
                newAddressing.append(facei);
                newFlipMap.append(doFlip[facei]);
            }
            zoneSet.addressing().transfer(newAddressing);
            zoneSet.flipMap().transfer(newFlipMap);
            zoneSet.updateSet();
        }
    }
}


// ************************************************************************* //
