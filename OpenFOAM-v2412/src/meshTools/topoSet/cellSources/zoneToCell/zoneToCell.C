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

#include "zoneToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToCell, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, zoneToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, zoneToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        zoneToCell,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        zoneToCell,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToCell::usage_
(
    zoneToCell::typeName,
    "\n    Usage: zoneToCell zone\n\n"
    "    Select all cells in the cellZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToCell::combine
(
    topoSet& set,
    const labelUList& zoneIDs,
    const bool add,
    const bool verbosity
) const
{
    const label nZones = mesh_.cellZones().size();

    if (zoneIDs.empty() || !nZones)
    {
        return;  // Nothing to do
    }

    for (const label zonei : zoneIDs)
    {
        if (zonei < 0 || zonei >= nZones)
        {
            continue;
        }

        const auto& zone = mesh_.cellZones()[zonei];

        if (verbosity)
        {
            Info<< "    Using zone " << zone.name() << " with "
                << returnReduce(zone.size(), sumOp<label>()) << " cells" << nl;
        }

        for (const label celli : zone)
        {
            // Only do active cells
            if (celli >= 0 && celli < mesh_.nCells())
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


void Foam::zoneToCell::combine(topoSet& set, const bool add) const
{
    if (!zoneIDs_.empty())
    {
        combine(set, zoneIDs_, add, false);
        return;
    }

    if (zoneMatcher_.empty())
    {
        return;  // Nothing to do
    }

    const labelList matched(mesh_.cellZones().indices(zoneMatcher_));

    if (matched.empty())
    {
        WarningInFunction
            << "Cannot find any cellZone matching "
            << flatOutput(zoneMatcher_) << nl
            << "Valid names: " << flatOutput(mesh_.cellZones().names())
            << endl;

        return;  // Nothing to do
    }

    combine(set, matched, add, verbose_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const wordRes& zoneSelector
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(zoneSelector)
{}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(one{}, zoneName)
{}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const labelUList& zoneIDs
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(),
    zoneIDs_(zoneIDs)
{}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    zoneMatcher_()
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", zoneMatcher_))
    {
        zoneMatcher_.resize(1);
        zoneMatcher_.front() = dict.getCompat<wordRe>("zone", {{"name", 1806}});
    }
}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    zoneToCell(mesh, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordRes& Foam::zoneToCell::zones() const noexcept
{
    return zoneMatcher_;
}


void Foam::zoneToCell::zones(const wordRes& zonesSelector)
{
    zoneMatcher_ = zonesSelector;
    zoneIDs_.clear();
}


void Foam::zoneToCell::zones(const wordRe& zoneName)
{
    zoneMatcher_.resize(1);
    zoneMatcher_.front() = zoneName;
    zoneIDs_.clear();
}


void Foam::zoneToCell::zones(const labelUList& zoneIDs)
{
    zoneMatcher_.clear();
    zoneIDs_ = zoneIDs;
}


void Foam::zoneToCell::zones(const label zoneID)
{
    zoneMatcher_.clear();
    zoneIDs_.resize(1);
    zoneIDs_.front() = zoneID;
}


void Foam::zoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Adding all cells of cell zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Removing all cells of cell zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
