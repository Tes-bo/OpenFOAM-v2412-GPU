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

#include "faceBoolSet.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(faceBoolSet);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceBoolSet::faceBoolSet(const polyMesh& mesh)
:
    faceBoolSet(mesh, false)
{}


Foam::faceBoolSet::faceBoolSet(const polyMesh& mesh, const bool val)
:
    topoBoolSet(mesh, "faceBoolSet", mesh.nFaces(), val)
{}


Foam::faceBoolSet::faceBoolSet
(
    const polyMesh& mesh,
    const boolList& bools
)
:
    topoBoolSet(mesh, "faceBoolSet", mesh.nFaces(), bools)
{}


Foam::faceBoolSet::faceBoolSet
(
    const polyMesh& mesh,
    boolList&& bools
)
:
    topoBoolSet(mesh, "faceBoolSet", mesh.nFaces(), std::move(bools))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceBoolSet::sync(const polyMesh& mesh)
{
    syncTools::syncFaceList(mesh, selected_, orEqOp<bool>());
}


Foam::label Foam::faceBoolSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


void Foam::faceBoolSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseFaceMap());
}


void Foam::faceBoolSet::distribute(const mapDistributePolyMesh& map)
{
    map.distributeFaceData(selected_);
}


void Foam::faceBoolSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.faceCentres(), maxLen);
}


// ************************************************************************* //
