/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "zoneSubSet.H"
#include "cellBitSet.H"
#include "haloToCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{
    defineTypeNameAndDebug(zoneSubSet, 0);
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Detail::zoneSubSet::correct()
{
    subsetter_.clear();
    haloCells_.clearStorage();

    if (zoneMatcher_.empty())
    {
        return false;
    }

    // Select named zones
    cellBitSet selectedCells
    (
        subsetter_.baseMesh(),
        subsetter_.baseMesh().cellZones().selection(zoneMatcher_)
    );

    if (debug)
    {
        Pout<< "Subsetting "
            << selectedCells.addressing().count()
            << " cells based on cellZones "
            << flatOutput(zoneMatcher_) << endl;
    }

    if (nLayers_ > 0)
    {
        // Add halo layer(s)
        haloToCell haloSource(subsetter_.baseMesh(), nLayers_);
        haloSource.verbose(false);

        // Before adding halo cells
        haloCells_ = selectedCells.addressing();

        haloSource.applyToSet(topoSetSource::ADD, selectedCells);

        // Halo cells: anything new, not in the original set
        haloCells_ ^= selectedCells.addressing();
    }

    if (debug)
    {
        const label nHalo = haloCells_.count();
        const label nSubCell = selectedCells.addressing().count();

        Info<< "    overall "
            << returnReduce(nSubCell, sumOp<label>())
            << " cells after adding " << nLayers_ << " layers with "
            << returnReduce(nHalo, sumOp<label>())
            << " halo cells"
            << endl;
    }

    subsetter_.reset(selectedCells.addressing());

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Detail::zoneSubSet::zoneSubSet
(
    const fvMesh& mesh,
    const wordRes& zoneSelector,
    const label nZoneLayers
)
:
    subsetter_(mesh),
    zoneMatcher_(zoneSelector),
    nLayers_(nZoneLayers),
    haloCells_()
{
    correct();
}


Foam::Detail::zoneSubSet::zoneSubSet
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    subsetter_(mesh),
    zoneMatcher_(),
    nLayers_(dict.getOrDefault<label>("nLayers", 0)),
    haloCells_()
{
    dict.readIfPresent("cellZones", zoneMatcher_);

    correct();
}


// ************************************************************************* //
