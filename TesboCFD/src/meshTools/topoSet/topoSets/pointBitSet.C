/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "pointBitSet.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(pointBitSet);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointBitSet::pointBitSet(const polyMesh& mesh, const bool val)
:
    topoBitSet(mesh, "pointBitSet", mesh.nPoints(), val)
{}


Foam::pointBitSet::pointBitSet
(
    const polyMesh& mesh,
    const bitSet& bits
)
:
    topoBitSet(mesh, "pointBitSet", mesh.nPoints(), bits)
{}


Foam::pointBitSet::pointBitSet
(
    const polyMesh& mesh,
    bitSet&& bits
)
:
    topoBitSet(mesh, "pointBitSet", mesh.nPoints(), std::move(bits))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointBitSet::sync(const polyMesh& mesh)
{
    // The nullValue = '0u'
    syncTools::syncPointList(mesh, selected_, orEqOp<unsigned int>(), 0u);
}


Foam::label Foam::pointBitSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nPoints();
}


void Foam::pointBitSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reversePointMap());
}


void Foam::pointBitSet::distribute(const mapDistributePolyMesh& map)
{
    bitSet& labels = selected_;

    boolList contents(labels.values());

    map.distributePointData(contents);

    // The new length is contents.size();
    labels.assign(contents);
}


void Foam::pointBitSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.points(), maxLen);
}


// ************************************************************************* //
