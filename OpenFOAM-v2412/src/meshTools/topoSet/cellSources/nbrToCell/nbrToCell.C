/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "nbrToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nbrToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, nbrToCell, word);
    addToRunTimeSelectionTable(topoSetSource, nbrToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, nbrToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, nbrToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        nbrToCell,
        word,
        nbr
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        nbrToCell,
        istream,
        nbr
    );
}


Foam::topoSetSource::addToUsageTable Foam::nbrToCell::usage_
(
    nbrToCell::typeName,
    "\n    Usage: nbrToCell <nNeighbours>\n\n"
    "    Select all cells with <= nNeighbours neighbouring cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nbrToCell::combine(topoSet& set, const bool add) const
{
    if (minNbrs_ < 1)
    {
        return;  // Nothing to do
    }

    const cellList& cells = mesh().cells();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    boolList isCoupled(mesh_.nBoundaryFaces(), false);

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                isCoupled[facei-mesh_.nInternalFaces()] = true;
                ++facei;
            }
        }
    }

    forAll(cells, celli)
    {
        const cell& cFaces = cells[celli];

        label nNbrCells = 0;

        for (const label facei : cFaces)
        {
            if (mesh_.isInternalFace(facei))
            {
                ++nNbrCells;
            }
            else if (isCoupled[facei-mesh_.nInternalFaces()])
            {
                ++nNbrCells;
            }
        }

        if (nNbrCells <= minNbrs_)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const label minNbrs
)
:
    topoSetCellSource(mesh),
    minNbrs_(minNbrs)
{}


Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    minNbrs_(dict.getCheck<label>("neighbours", labelMinMax::ge(1)))
{}


Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    minNbrs_(readLabel(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nbrToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells with only " << minNbrs_
                << " or fewer neighbouring cells" << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells with only " << minNbrs_
                << " or fewer neighbouring cells" << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
