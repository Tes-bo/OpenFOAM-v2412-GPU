/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "boundaryToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, boundaryToCell, word);
    addToRunTimeSelectionTable(topoSetSource, boundaryToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, boundaryToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, boundaryToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        boundaryToCell,
        word,
        boundary
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        boundaryToCell,
        istream,
        boundary
    );
}


Foam::topoSetSource::addToUsageTable Foam::boundaryToCell::usage_
(
    boundaryToCell::typeName,
    "\n    Usage: boundaryToCell\n\n"
    "    Select all boundary cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryToCell::combine(topoSet& set, const bool add) const
{
    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        ++facei
    )
    {
        addOrDelete(set, mesh().faceOwner()[facei], add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryToCell::boundaryToCell(const polyMesh& mesh)
:
    topoSetCellSource(mesh)
{}


Foam::boundaryToCell::boundaryToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict)
{}


Foam::boundaryToCell::boundaryToCell
(
    const polyMesh& mesh,
    Istream&
)
:
    topoSetCellSource(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all boundary cells ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all boundary cells ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
