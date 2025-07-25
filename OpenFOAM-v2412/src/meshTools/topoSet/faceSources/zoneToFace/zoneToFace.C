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

#include "zoneToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToFace, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToFace, istream);

    addToRunTimeSelectionTable(topoSetFaceSource, zoneToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, zoneToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        zoneToFace,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        zoneToFace,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToFace::usage_
(
    zoneToFace::typeName,
    "\n    Usage: zoneToFace zone\n\n"
    "    Select all faces in the faceZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToFace::combine
(
    topoSet& set,
    const labelUList& zoneIDs,
    const bool add,
    const bool verbosity
) const
{
    const label nZones = mesh_.faceZones().size();

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

        const auto& zone = mesh_.faceZones()[zonei];

        if (verbosity)
        {
            Info<< "    Using zone " << zone.name() << " with "
                << returnReduce(zone.size(), sumOp<label>()) << " faces"
                << endl;
        }

        for (const label facei : zone)
        {
            // Only do active faces
            if (facei >= 0 && facei < mesh_.nFaces())
            {
                addOrDelete(set, facei, add);
            }
        }
    }
}


void Foam::zoneToFace::combine(topoSet& set, const bool add) const
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

    const labelList matched(mesh_.faceZones().indices(zoneMatcher_));

    if (matched.empty())
    {
        WarningInFunction
            << "Cannot find any faceZone matching "
            << flatOutput(zoneMatcher_) << nl
            << "Valid names: " << flatOutput(mesh_.faceZones().names())
            << endl;

        return;  // Nothing to do
    }

    combine(set, matched, add, verbose_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const wordRes& zoneSelector
)
:
    topoSetFaceSource(mesh),
    zoneMatcher_(zoneSelector)
{}


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetFaceSource(mesh),
    zoneMatcher_(one{}, zoneName)
{}


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const labelUList& zoneIDs
)
:
    topoSetFaceSource(mesh),
    zoneMatcher_(),
    zoneIDs_(zoneIDs)
{}


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    zoneMatcher_()
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", zoneMatcher_))
    {
        zoneMatcher_.resize(1);
        zoneMatcher_.front() = dict.getCompat<wordRe>("zone", {{"name", 1806}});
    }
}


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    zoneToFace(mesh, wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordRes& Foam::zoneToFace::zones() const noexcept
{
    return zoneMatcher_;
}


void Foam::zoneToFace::zones(const wordRes& zonesSelector)
{
    zoneMatcher_ = zonesSelector;
    zoneIDs_.clear();
}


void Foam::zoneToFace::zones(const wordRe& zoneName)
{
    zoneMatcher_.resize(1);
    zoneMatcher_.front() = zoneName;
    zoneIDs_.clear();
}


void Foam::zoneToFace::zones(const labelUList& zoneIDs)
{
    zoneMatcher_.clear();
    zoneIDs_ = zoneIDs;
}


void Foam::zoneToFace::zones(const label zoneID)
{
    zoneMatcher_.clear();
    zoneIDs_.resize(1);
    zoneIDs_.front() = zoneID;
}


void Foam::zoneToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Adding all faces of face zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Removing all faces of face zones: "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
