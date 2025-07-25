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

#include "faceToCell.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, faceToCell, word);
    addToRunTimeSelectionTable(topoSetSource, faceToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, faceToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, faceToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::faceToCell::usage_
(
    faceToCell::typeName,
    "\n    Usage: faceToCell <faceSet> neighbour|owner|any|all\n\n"
    "    Select cells that are the owner|neighbour|any"
    " of the faces in the faceSet or where all faces are in the faceSet\n\n"
);

const Foam::Enum
<
    Foam::faceToCell::faceAction
>
Foam::faceToCell::faceActionNames_
({
    { faceAction::ANY, "any" },
    { faceAction::ALL, "all" },
    { faceAction::OWNER, "owner" },
    { faceAction::NEIGHBOUR, "neighbour" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Selector>
void Foam::faceToCell::combineImpl
(
    topoSet& set,
    const bool add,
    const Selector& faceLabels
) const
{
    // Handle owner/neighbour/any selection
    for (const label facei : faceLabels)
    {
        if ((option_ == OWNER) || (option_ == ANY))
        {
            const label celli = mesh_.faceOwner()[facei];

            addOrDelete(set, celli, add);
        }

        if (mesh_.isInternalFace(facei))
        {
            if ((option_ == NEIGHBOUR) || (option_ == ANY))
            {
                const label celli = mesh_.faceNeighbour()[facei];

                addOrDelete(set, celli, add);
            }
        }
    }

    // Handle all selection.
    if (option_ == ALL)
    {
        // Count number of selected faces per cell.

        Map<label> facesPerCell(faceLabels.size());

        for (const label facei : faceLabels)
        {
            // Count faces on owner
            ++(facesPerCell(mesh_.faceOwner()[facei], 0));

            if (mesh_.isInternalFace(facei))
            {
                // Count faces on neighbour
                ++(facesPerCell(mesh_.faceNeighbour()[facei], 0));
            }
        }

        // Include cells that are referenced as many times as they have faces
        // -> all faces in set.
        forAllConstIters(facesPerCell, iter)
        {
            const label celli = iter.key();
            const label count = iter.val();

            if (count == mesh_.cells()[celli].size())
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


void Foam::faceToCell::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    if (isZone_)
    {
        const labelList& faceLabels = mesh_.faceZones()[setName];

        combineImpl(set, add, faceLabels);
    }
    else
    {
        // Load the set
        faceSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
        const labelHashSet& faceLabels = loadedSet;

        combineImpl(set, add, faceLabels);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    const word& setName,
    const faceAction option
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false),
    option_(option)
{}


Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_)),
    option_(faceActionNames_.get("option", dict))
{}


Foam::faceToCell::faceToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells according to face "
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
            Info<< "    Removing cells according to face "
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
