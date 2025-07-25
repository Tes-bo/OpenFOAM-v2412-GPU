/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::removeZones()
{
    DebugInFunction << "Removing surface zones." << endl;

    // Remove the surface zones
    surfZones_.clear();

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::clearGeom()
{
    DebugInFunction << "Clearing geometric data" << endl;

    MeshReference::clearGeom();
}


void Foam::surfMesh::clearAddressing()
{
    DebugInFunction << "Clearing topology" << endl;

    MeshReference::clearPatchMeshAddr();
}


void Foam::surfMesh::clearOut()
{
    MeshReference::clearOut();

    clearGeom();
    clearAddressing();
}


void Foam::surfMesh::clearFields()
{
    // Clear the entire registry
    surfaceRegistry::clear();
}


// ************************************************************************* //
