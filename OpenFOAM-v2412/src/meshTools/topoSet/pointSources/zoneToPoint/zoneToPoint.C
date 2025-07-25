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

#include "zoneToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, istream);

    addToRunTimeSelectionTable(topoSetPointSource, zoneToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, zoneToPoint, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        zoneToPoint,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        zoneToPoint,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToPoint::usage_
(
    zoneToPoint::typeName,
    "\n    Usage: zoneToPoint zone\n\n"
    "    Select all points in the pointZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToPoint::combine
(
    topoSet& set,
    const labelUList& zoneIDs,
    const bool add,
    const bool verbosity
) const
{
    const label nZones = mesh_.pointZones().size();

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

        const auto& zone = mesh_.pointZones()[zonei];

        if (verbosity)
        {
            Info<< "    Using zone " << zone.name() << " with "
                << returnReduce(zone.size(), sumOp<label>()) << " points"
                << endl;
        }

        for (const label pointi : zone)
        {
            // Only do active points
            if (pointi >= 0 && pointi < mesh_.nPoints())
            {
                addOrDelete(set, pointi, add);
            }
        }
    }
}


void Foam::zoneToPoint::combine(topoSet& set, const bool add) const
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

    const labelList matched(mesh_.pointZones().indices(zoneMatcher_));

    if (matched.empty())
    {
        WarningInFunction
            << "Cannot find any pointZone matching "
            << flatOutput(zoneMatcher_) << nl
            << "Valid names: " << flatOutput(mesh_.pointZones().names())
            << endl;

        return;  // Nothing to do
    }

    combine(set, matched, add, verbose_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const wordRes& zoneSelector
)
:
    topoSetPointSource(mesh),
    zoneMatcher_(zoneSelector)
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetPointSource(mesh),
    zoneMatcher_(one{}, zoneName)
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const labelUList& zoneIDs
)
:
    topoSetPointSource(mesh),
    zoneMatcher_(),
    zoneIDs_(zoneIDs)
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    zoneMatcher_()
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", zoneMatcher_))
    {
        zoneMatcher_.resize(1);
        zoneMatcher_.front() = dict.getCompat<wordRe>("zone", {{"name", 1806}});
    }
}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    zoneToPoint(mesh, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordRes& Foam::zoneToPoint::zones() const noexcept
{
    return zoneMatcher_;
}


void Foam::zoneToPoint::zones(const wordRes& zonesSelector)
{
    zoneMatcher_ = zonesSelector;
    zoneIDs_.clear();
}


void Foam::zoneToPoint::zones(const wordRe& zoneName)
{
    zoneMatcher_.resize(1);
    zoneMatcher_.front() = zoneName;
    zoneIDs_.clear();
}


void Foam::zoneToPoint::zones(const labelUList& zoneIDs)
{
    zoneMatcher_.clear();
    zoneIDs_ = zoneIDs;
}


void Foam::zoneToPoint::zones(const label zoneID)
{
    zoneMatcher_.clear();
    zoneIDs_.resize(1);
    zoneIDs_.front() = zoneID;
}


void Foam::zoneToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Adding all points of point zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Removing all points of point zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
