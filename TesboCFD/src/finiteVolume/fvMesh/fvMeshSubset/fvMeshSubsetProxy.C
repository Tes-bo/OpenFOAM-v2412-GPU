/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "fvMeshSubsetProxy.H"
#include "cellSet.H"
#include "cellZone.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshSubsetProxy::clearOut()
{
    subsetter_.clear();
    type_ = subsetType::NONE;
    name_.clear();
    names_.clear();
    selectedCells_.clearStorage();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubsetProxy::fvMeshSubsetProxy(fvMesh& baseMesh)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(-1),
    type_(subsetType::NONE),
    name_(),
    names_(),
    selectedCells_()
{}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const subsetType type,
    const word& selectionName,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(selectionName.empty() ? subsetType::NONE : type),
    name_(),
    names_(),
    selectedCells_()
{
    if (type_ == subsetType::ZONES)
    {
        // Populate wordRes for ZONES
        names_.resize(1);
        names_.first() = selectionName;
    }
    else if (type_ == subsetType::SET || type_ == subsetType::ZONE)
    {
        name_ = selectionName;
    }

    correct();
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const wordRes& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(subsetType::ZONES),
    name_(),
    names_(zoneNames),
    selectedCells_()
{
    correct();
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    wordRes&& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(subsetType::ZONES),
    name_(),
    names_(std::move(zoneNames)),
    selectedCells_()
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubsetProxy::resetZones(const wordRes& zoneNames)
{
    fvMeshSubsetProxy::clearOut();

    if (!zoneNames.empty())
    {
        type_ = subsetType::ZONES;
        names_ = zoneNames;
        correct();
    }
}


bool Foam::fvMeshSubsetProxy::correct(bool verbose)
{
    if (type_ == subsetType::NONE)
    {
        subsetter_.clear();
        selectedCells_.clearStorage();
        return false;
    }

    const label nCells = baseMesh_.nCells();

    bitSet selectedCells;

    if (type_ == subsetType::SET)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet cset(baseMesh_, name_, IOobject::NO_REGISTER);

        selectedCells.resize(nCells);
        for (const label idx : cset)
        {
            selectedCells.set(idx);
        }
    }
    else if (type_ == subsetType::ZONE)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        selectedCells.resize(nCells);
        selectedCells.set(baseMesh_.cellZones()[name_]);
    }
    else if (type_ == subsetType::ZONES)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZones "
                << flatOutput(names_) << endl;
        }

        selectedCells = baseMesh_.cellZones().selection(names_);
    }


    const bool changed = (selectedCells_ != selectedCells);

    // Use as a cached value for next time
    selectedCells_.transfer(selectedCells);

    if (changed || selectedCells_.empty())
    {
        subsetter_.reset(selectedCells_, exposedPatchId_);
    }

    return returnReduceOr(changed);
}


Foam::polyMesh::readUpdateState Foam::fvMeshSubsetProxy::readUpdate()
{
    const polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if (meshState == polyMesh::POINTS_MOVED)
    {
        if (correct(true))
        {
            // The cellSet/cellZone changed on POINTS_MOVED,
            // treat like TOPO_CHANGE
            return polyMesh::TOPO_CHANGE;
        }
    }
    else if
    (
        meshState == polyMesh::TOPO_CHANGE
     || meshState == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        correct(true);
    }

    return meshState;
}


// ************************************************************************* //
