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

#include "cellToCell.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToCell, word);
    addToRunTimeSelectionTable(topoSetSource, cellToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, cellToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, cellToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cellToCell::usage_
(
    cellToCell::typeName,
    "\n    Usage: cellToCell <cellSet>\n\n"
    "    Select all cells in the cellSet\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToCell::cellToCell
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false)
{}


Foam::cellToCell::cellToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_))
{}


Foam::cellToCell::cellToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all elements of cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.addSet(mesh_.cellZones()[setName]);
            }
            else
            {
                cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
                set.addSet(loadedSet);
            }
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all elements of cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.subtractSet(mesh_.cellZones()[setName]);
            }
            else
            {
                cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
                set.subtractSet(loadedSet);
            }
        }
    }
}


// ************************************************************************* //
