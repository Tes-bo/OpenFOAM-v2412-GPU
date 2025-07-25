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

#include "pointToCell.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, word);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, pointToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, pointToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::pointToCell::usage_
(
    pointToCell::typeName,
    "\n    Usage: pointToCell <pointSet> any|edge\n\n"
    "    Select all cells with any point ('any') or any edge ('edge')"
    " in the pointSet\n\n"
);

const Foam::Enum
<
    Foam::pointToCell::pointAction
>
Foam::pointToCell::pointActionNames_
({
    { pointAction::ANY, "any" },
    { pointAction::EDGE, "edge" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Selector>
void Foam::pointToCell::combineImpl
(
    topoSet& set,
    const bool add,
    const Selector& pointLabels
) const
{
    // Handle any selection
    if (option_ == ANY)
    {
        for (const label pointi : pointLabels)
        {
            const labelList& pCells = mesh_.pointCells()[pointi];

            addOrDelete(set, pCells, add);
        }
    }
    else if (option_ == EDGE)
    {
        const faceList& faces = mesh_.faces();

        forAll(faces, facei)
        {
            const face& f = faces[facei];

            forAll(f, fp)
            {
                if
                (
                    pointLabels.found(f[fp])
                 && pointLabels.found(f.nextLabel(fp))
                )
                {
                    addOrDelete(set, mesh_.faceOwner()[facei], add);
                    if (mesh_.isInternalFace(facei))
                    {
                        addOrDelete(set, mesh_.faceNeighbour()[facei], add);
                    }
                }
            }
        }
    }
}


void Foam::pointToCell::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    if (isZone_)
    {
        const labelList& pointLabels = mesh_.pointZones()[setName];

        combineImpl(set, add, pointLabels);
    }
    else
    {
        // Load the set
        pointSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
        const labelHashSet& pointLabels = loadedSet;

        combineImpl(set, add, pointLabels);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false),
    option_(option)
{}


Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_)),
    option_(pointActionNames_.get("option", dict))
{}


Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells according to point "
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
            Info<< "    Removing cells according to point "
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
