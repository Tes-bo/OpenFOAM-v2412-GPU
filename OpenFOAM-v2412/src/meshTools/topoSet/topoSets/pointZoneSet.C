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

#include "pointZoneSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(pointZoneSet);
    addToRunTimeSelectionTable(topoSet, pointZoneSet, word);
    addToRunTimeSelectionTable(topoSet, pointZoneSet, size);
    addToRunTimeSelectionTable(topoSet, pointZoneSet, set);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointZoneSet::updateSet()
{
    Foam::sort(addressing_);

    pointSet::clearStorage();
    pointSet::reserve(addressing_.size());
    pointSet::set(addressing_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label initialCapacity,
    IOobjectOption::writeOption wOpt
)
:
    pointSet(mesh, name, initialCapacity, wOpt),  // Construct no-read
    mesh_(mesh)
{}


Foam::pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt
)
:
    pointZoneSet(mesh, name, label(0), wOpt)  // Construct no-read
{
    const auto& zones = mesh.pointZones();
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
    check(mesh.nPoints());
}


Foam::pointZoneSet::pointZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    IOobjectOption::writeOption wOpt
)
:
    pointZoneSet(mesh, name, label(0), wOpt)  // Construct no-read
{
    const auto* zonePtr = isA<pointZoneSet>(set);

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

void Foam::pointZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label id = 0; id < maxLen; ++id)
    {
        if (!found(id))
        {
            ++n;
        }
    }

    // Fill
    addressing_.resize_nocopy(n);
    n = 0;

    for (label id = 0; id < maxLen; ++id)
    {
        if (!found(id))
        {
            addressing_[n] = id;
            ++n;
        }
    }
    updateSet();
}


void Foam::pointZoneSet::subset(const labelUList& elems)
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


void Foam::pointZoneSet::subset(const topoSet& set)
{
    const auto* zonePtr = isA<pointZoneSet>(set);

    if (zonePtr)
    {
        // Is a pointZoneSet
        this->subset(zonePtr->addressing());
    }
    else
    {
        // Assume a pointSet
        this->subset(refCast<const pointSet>(set).sortedToc());
    }
}


void Foam::pointZoneSet::addSet(const labelUList& elems)
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


void Foam::pointZoneSet::addSet(const topoSet& set)
{
    const auto* zonePtr = isA<pointZoneSet>(set);

    if (zonePtr)
    {
        // Is a pointZoneSet
        this->addSet(zonePtr->addressing());
    }
    else
    {
        // Assume a pointSet
        this->addSet(refCast<const pointSet>(set).sortedToc());
    }
}


void Foam::pointZoneSet::subtractSet(const labelUList& elems)
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


void Foam::pointZoneSet::subtractSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

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


void Foam::pointZoneSet::sync(const polyMesh& mesh)
{
    pointSet::sync(mesh);

    // Take over contents of pointSet into addressing.
    addressing_ = sortedToc();
    updateSet();
}


Foam::label Foam::pointZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nPoints();
}


bool Foam::pointZoneSet::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Write shadow pointSet
    const word oldTypeName = typeName;
    const_cast<word&>(type()) = pointSet::typeName;
    bool ok = pointSet::writeObject(streamOpt, writeOnProc);
    const_cast<word&>(type()) = oldTypeName;

    // Modify pointZone
    auto& zones = const_cast<polyMesh&>(mesh_).pointZones();
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


void Foam::pointZoneSet::updateMesh(const mapPolyMesh& morphMap)
{
    DynamicList<label> newAddressing(addressing_.size());

    for (const label pointi : addressing_)
    {
        const label newPointi = morphMap.reversePointMap()[pointi];
        if (newPointi >= 0)
        {
            newAddressing.push_back(newPointi);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::pointZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    pointSet::writeDebug(os, mesh, maxLen);
}


// ************************************************************************* //
