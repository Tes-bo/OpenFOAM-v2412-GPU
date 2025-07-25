/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "ensightWrite.H"
#include "dictionary.H"
#include "cellBitSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::updateSubset
(
    fvMeshSubset& subsetter
) const
{
    if (selection_.empty())
    {
        return false;
    }

    bitSet selectedCells
    (
        cellBitSet::select(subsetter.baseMesh(), selection_)
    );

    subsetter.reset(selectedCells);

    return true;
}


bool Foam::functionObjects::ensightWrite::update()
{
    if (meshState_ == polyMesh::UNCHANGED)
    {
        return false;
    }

    // Even if selection doesn't nominally change, a new subMesh is built.
    // Must update the reference for ensightMesh
    if (meshSubset_.hasSubMesh())
    {
        ensMesh_.clear();
        meshSubset_.clear();
    }
    // Probably not needed...
    /// else if (ensMesh_)
    /// {
    ///     ensMesh_->expire();
    /// }

    updateSubset(meshSubset_);

    if (!ensMesh_)
    {
        ensMesh_.reset(new ensightMesh(meshSubset_.mesh(), writeOpts_));
    }
    else if (ensMesh_->needsUpdate())
    {
        ensMesh_->correct();
    }

    meshState_ = polyMesh::UNCHANGED;
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::readSelection(const dictionary& dict)
{
    // Ensure consistency
    ensMesh_.clear();

    meshSubset_.clear();
    meshState_ = polyMesh::TOPO_CHANGE;

    selectFields_.clear();
    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    blockFields_.clear();
    dict.readIfPresent("excludeFields", blockFields_);
    blockFields_.uniq();

    // Actions to define selection
    selection_ = dict.subOrEmptyDict("selection");

    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh&)
{
    meshState_ = polyMesh::TOPO_CHANGE;
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh&)
{
    // Only move to worse states
    if (meshState_ == polyMesh::UNCHANGED)
    {
        meshState_ = polyMesh::POINTS_MOVED;
    }
}


// ************************************************************************* //
