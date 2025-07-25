/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "faceBitSet.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(faceBitSet);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceBitSet::faceBitSet(const polyMesh& mesh, const bool val)
:
    topoBitSet(mesh, "faceBitSet", mesh.nFaces(), val)
{}


Foam::faceBitSet::faceBitSet
(
    const polyMesh& mesh,
    const bitSet& bits
)
:
    topoBitSet(mesh, "faceBitSet", mesh.nFaces(), bits)
{}


Foam::faceBitSet::faceBitSet
(
    const polyMesh& mesh,
    bitSet&& bits
)
:
    topoBitSet(mesh, "faceBitSet", mesh.nFaces(), std::move(bits))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceBitSet::sync(const polyMesh& mesh)
{
    syncTools::syncFaceList(mesh, selected_, orEqOp<unsigned int>());
}


Foam::label Foam::faceBitSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


void Foam::faceBitSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseFaceMap());
}


void Foam::faceBitSet::distribute(const mapDistributePolyMesh& map)
{
    bitSet& labels = selected_;

    boolList contents(labels.values());

    map.distributeFaceData(contents);

    // The new length is contents.size();
    labels.assign(contents);
}


void Foam::faceBitSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.faceCentres(), maxLen);
}


// ************************************************************************* //
