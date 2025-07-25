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

#include "cellToPoint.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, cellToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, cellToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, cellToPoint, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cellToPoint::usage_
(
    cellToPoint::typeName,
    "\n    Usage: cellToPoint <cellSet> all\n\n"
    "    Select all points of cells in the cellSet\n\n"
);

const Foam::Enum
<
    Foam::cellToPoint::cellAction
>
Foam::cellToPoint::cellActionNames_
({
    { cellAction::ALL, "all" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Selector>
void Foam::cellToPoint::combineImpl
(
    topoSet& set,
    const bool add,
    const Selector& cellLabels
) const
{
    // Add all points from given cells
    for (const label celli : cellLabels)
    {
        const labelList& cFaces = mesh_.cells()[celli];

        for (const label facei : cFaces)
        {
            const face& f = mesh_.faces()[facei];

            addOrDelete(set, f, add);
        }
    }
}


void Foam::cellToPoint::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    if (isZone_)
    {
        const labelList& cellLabels = mesh_.cellZones()[setName];

        combineImpl(set, add, cellLabels);
    }
    else
    {
        // Load the set
        cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
        const labelHashSet& cellLabels = loadedSet;

        combineImpl(set, add, cellLabels);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    const word& setName,
    const cellAction option
)
:
    topoSetPointSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false),
    option_(option)
{}


Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_)),
    option_(cellActionNames_.get("option", dict))
{}


Foam::cellToPoint::cellToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false),
    option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding points in cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, true, setName);
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing points in cell "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, false, setName);
        }
    }
}


// ************************************************************************* //
