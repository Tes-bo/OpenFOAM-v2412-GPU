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

#include "cellZoneSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(cellZoneSet);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, word);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, size);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, set);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellZoneSet::updateSet()
{
    Foam::sort(addressing_);

    cellSet::clearStorage();
    cellSet::reserve(addressing_.size());
    cellSet::set(addressing_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label initialCapacity,
    IOobjectOption::writeOption wOpt
)
:
    cellSet(mesh, name, initialCapacity, wOpt),  // Construct no-read
    mesh_(mesh)
{}


Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt
)
:
    cellZoneSet(mesh, name, label(0), wOpt)  // Construct no-read
{
    const auto& zones = mesh.cellZones();
    const auto* zonePtr = zones.cfindZone(name);

    if (!zonePtr)
    {
        if (IOobjectOption::isReadRequired(rOpt))
        {
            FatalErrorInFunction
                << "Zone named " << name << " not found.  "
                << "List of available zone names: " << zones.names() << nl
                << exit(FatalError);
        }
    }
    else if (IOobjectOption::isAnyRead(rOpt))
    {
        const auto& zn = *zonePtr;
        addressing_ = zn;
    }

    updateSet();
    check(mesh.nCells());
}


Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    IOobjectOption::writeOption wOpt
)
:
    cellZoneSet(mesh, name, label(0), wOpt)  // Construct no-read
{
    const auto* zonePtr = isA<cellZoneSet>(set);

    if (zonePtr)
    {
        addressing_ = zonePtr->addressing();
    }
    else
    {
        addressing_ = set.sortedToc();
    }

    updateSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label celli = 0; celli < maxLen; ++celli)
    {
        if (!found(celli))
        {
            ++n;
        }
    }

    // Fill
    addressing_.resize_nocopy(n);
    n = 0;

    for (label celli = 0; celli < maxLen; ++celli)
    {
        if (!found(celli))
        {
            addressing_[n] = celli;
            ++n;
        }
    }

    updateSet();
}


void Foam::cellZoneSet::subset(const labelUList& elems)
{
    DynamicList<label> newAddressing(addressing_.size());

    for (const label id : elems)
    {
        if (found(id))
        {
            newAddressing.push_back(id);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::subset(const topoSet& set)
{
    const auto* zonePtr = isA<cellZoneSet>(set);

    if (zonePtr)
    {
        // Is a cellZoneSet
        this->subset(zonePtr->addressing());
    }
    else
    {
        // Assume a cellSet
        this->subset(refCast<const cellSet>(set).sortedToc());
    }
}


void Foam::cellZoneSet::addSet(const labelUList& elems)
{
    DynamicList<label> newAddressing(addressing_);

    for (const label id : elems)
    {
        if (!found(id))
        {
            newAddressing.push_back(id);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::addSet(const topoSet& set)
{
    const auto* zonePtr = isA<cellZoneSet>(set);

    if (zonePtr)
    {
        // Is a cellZoneSet
        this->addSet(zonePtr->addressing());
    }
    else
    {
        // Assume a cellSet
        this->addSet(refCast<const cellSet>(set).sortedToc());
    }
}


void Foam::cellZoneSet::subtractSet(const labelUList& elems)
{
    DynamicList<label> newAddressing(addressing_.size());

    const labelHashSet set(elems);

    for (const label id : addressing_)
    {
        if (!set.found(id))
        {
            // Retain if not in the topoSet (parameter)
            newAddressing.push_back(id);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::subtractSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    for (const label id : addressing_)
    {
        if (!set.found(id))
        {
            // Not found in zoneSet so add
            newAddressing.push_back(id);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::sync(const polyMesh& mesh)
{
    cellSet::sync(mesh);

    // Take over contents of cellSet into addressing.
    addressing_ = sortedToc();
    updateSet();
}


Foam::label Foam::cellZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


bool Foam::cellZoneSet::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Write shadow cellSet
    const word oldTypeName = typeName;
    const_cast<word&>(type()) = cellSet::typeName;
    bool ok = cellSet::writeObject(streamOpt, writeOnProc);
    const_cast<word&>(type()) = oldTypeName;

    // Modify cellZone
    auto& zones = const_cast<polyMesh&>(mesh_).cellZones();
    auto* zonePtr = zones.findZone(name());

    if (zonePtr)
    {
        zonePtr->resetAddressing(addressing_);
    }
    else
    {
        zones.emplace_back
        (
            name(),
            addressing_,
            zones.size(),  // zoneID
            zones
        );
    }
    zones.clearAddressing();

    return ok && zones.write(writeOnProc);
}


void Foam::cellZoneSet::updateMesh(const mapPolyMesh& morphMap)
{
    DynamicList<label> newAddressing(addressing_.size());

    for (const label celli : addressing_)
    {
        label newCelli = morphMap.reverseCellMap()[celli];
        if (newCelli >= 0)
        {
            newAddressing.push_back(newCelli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    cellSet::writeDebug(os, mesh, maxLen);
}


// ************************************************************************* //
