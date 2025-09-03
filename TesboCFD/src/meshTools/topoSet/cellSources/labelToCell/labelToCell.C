/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "labelToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(labelToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, labelToCell, word);
    addToRunTimeSelectionTable(topoSetSource, labelToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, labelToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, labelToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        labelToCell,
        word,
        label
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        labelToCell,
        istream,
        label
    );
}


Foam::topoSetSource::addToUsageTable Foam::labelToCell::usage_
(
    labelToCell::typeName,
    "\n    Usage: labelToCell (i0 i1 .. in)\n\n"
    "    Select cells by cellLabel\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelToCell::labelToCell
(
    const polyMesh& mesh,
    const labelList& labels
)
:
    topoSetCellSource(mesh),
    labels_(labels)
{}


Foam::labelToCell::labelToCell
(
    const polyMesh& mesh,
    labelList&& labels
)
:
    topoSetCellSource(mesh),
    labels_(std::move(labels))
{}


Foam::labelToCell::labelToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    labels_(dict.get<labelList>("value"))
{}


Foam::labelToCell::labelToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    labels_(checkIs(is))
{
    check(labels_, mesh.nCells());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells mentioned in dictionary"
                << " ..." << endl;
        }

        addOrDelete(set, labels_, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells mentioned in dictionary"
                << " ..." << endl;
        }

        addOrDelete(set, labels_, false);
    }
}


// ************************************************************************* //
