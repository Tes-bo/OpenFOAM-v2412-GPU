/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "haloToCell.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "topoBitSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(haloToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, haloToCell, word);
    addToRunTimeSelectionTable(topoSetSource, haloToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, haloToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, haloToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        haloToCell,
        word,
        halo
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        haloToCell,
        istream,
        halo
    );
}


Foam::topoSetSource::addToUsageTable Foam::haloToCell::usage_
(
    haloToCell::typeName,
    "\n    Usage: haloToCell\n\n"
    "    Select halo cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::haloToCell::combine(topoSet& set, const bool add) const
{
    if (steps_ < 1)
    {
        return;  // Nothing to do
    }

    const cellList& cells = mesh_.cells();
    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();


    // The starting set of cells
    bitSet current(cells.size());

    const auto* topoBitsPtr = isA<topoBitSet>(set);

    if (topoBitsPtr)
    {
        current |= topoBitsPtr->addressing();
    }
    else
    {
        for (const label celli : set)
        {
            current.set(celli);
        }
    }

    // The perimeter faces of the cell set
    bitSet outsideFaces(mesh_.nFaces());

    bitSet updates(cells.size());

    for (label stepi = 0; stepi < steps_; ++stepi)
    {
        // Mark up perimeter faces. Each mesh face is attached exactly
        // (0,1,2) times to a cell in the set. Using flip() each time means
        // the only 'on' bits correspond to faces that are attached once to
        // the cell set - ie, faces on the perimeter of the set.

        outsideFaces.reset();
        for (const label celli : current)
        {
            for (const label facei : cells[celli])
            {
                outsideFaces.flip(facei);
            }
        }

        // Use xor to eliminate perimeter faces that are actually attached
        // on both sides of the interface.

        syncTools::syncFaceList
        (
            mesh_,
            outsideFaces,
            bitXorEqOp<unsigned int>()
        );

        // Select all cells attached to the perimeter faces.
        updates.reset();
        for (const label facei : outsideFaces)
        {
            updates.set(faceOwn[facei]);
            if (mesh_.isInternalFace(facei))
            {
                updates.set(faceNei[facei]);
            }
        }

        if (add)
        {
            // Restrict to cells not already in the current set
            updates -= current;

            if (verbose_)
            {
                Info<< "    Grow "
                    << returnReduce(current.count(), sumOp<label>()) << " by "
                    << returnReduce(updates.count(), sumOp<label>()) << endl;
            }

            // Add to current set for the next loop
            current |= updates;
        }
        else
        {
            // Restrict to cells already in the current set
            updates &= current;

            if (verbose_)
            {
                Info<< "    Shrink "
                    << returnReduce(current.count(), sumOp<label>()) << " by "
                    << returnReduce(updates.count(), sumOp<label>()) << endl;
            }

            // Remove from current set for the next loop
            current -= updates;
        }

        // Could have early exit, but needs to be parallel-synchronized
        // if (returnReduceAnd(updates.none()))
        // {
        //     break;
        // }

        addOrDelete(set, updates, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::haloToCell::haloToCell
(
    const polyMesh& mesh,
    const label nsteps
)
:
    topoSetCellSource(mesh),
    steps_(nsteps)
{}


Foam::haloToCell::haloToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    steps_(dict.getOrDefault<label>("steps", 1))
{}


Foam::haloToCell::haloToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    haloToCell(mesh, 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::haloToCell::steps() const noexcept
{
    return steps_;
}


Foam::label Foam::haloToCell::steps(const label nsteps) noexcept
{
    label old(steps_);
    steps_ = nsteps;
    return old;
}


void Foam::haloToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    action=new option is not available for haloToCell" << nl
                << "    Cannot create new of halo (needs a starting set)"
                << endl;
        }

        set.clear();
    }
    else if (action == topoSetSource::ADD)
    {
        if (verbose_)
        {
            Info<< "    Adding halo cells to the current set, using "
                << steps_ << " step ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells on the perimeter of current set, using "
                << steps_ << " step ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
