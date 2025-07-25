/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faMeshMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshMapper::faMeshMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    nOldPoints_(mesh.nPoints()),
    nOldEdges_(mesh.nEdges()),
    nOldInternalEdges_(mesh.nInternalEdges()),
    nOldFaces_(mesh.nFaces()),
    oldPatchSizes_(),
    oldPatchStarts_(),
    oldPatchEdgeFaces_(),
    areaMap_(mesh, mpm),
    edgeMap_(mesh, mpm),
    boundaryMap_(mesh, mpm)
{
    // Capture old patch information
    const faBoundaryMesh& patches = mesh.boundary();

    oldPatchSizes_.resize(patches.size());
    oldPatchStarts_.resize(patches.size());
    oldPatchEdgeFaces_.resize(patches.size());

    forAll(patches, patchI)
    {
        oldPatchSizes_[patchI] = patches[patchI].size();
        oldPatchStarts_[patchI] = patches[patchI].start();

        oldPatchEdgeFaces_[patchI] = patches[patchI].edgeFaces();
    }
}


// ************************************************************************* //
