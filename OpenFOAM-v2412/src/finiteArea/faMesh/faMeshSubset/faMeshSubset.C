/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "faMeshSubset.H"
#include "boolList.H"
#include "BitOps.H"
#include "Pstream.H"
#include "emptyFaPatch.H"
#include "cyclicFaPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::faMeshSubset::exposedPatchName("oldInternalEdges");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::faMeshSubset::checkHasSubMesh() const
{
    if (!subMeshPtr_)
    {
        FatalErrorInFunction
            << "Mesh is not subsetted!" << nl
            << abort(FatalError);

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshSubset::faMeshSubset(const faMesh& baseMesh)
:
    baseMesh_(baseMesh),
    subMeshPtr_(nullptr),
    edgeFlipMapPtr_(nullptr),
    pointMap_(),
    faceMap_(),
    cellMap_(),
    patchMap_()
{}


Foam::faMeshSubset::faMeshSubset(const faMesh& baseMesh, const Foam::zero)
:
    faMeshSubset(baseMesh)
{
    reset(Foam::zero{});
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshSubset::clear()
{
    subMeshPtr_.reset(nullptr);
    edgeFlipMapPtr_.reset(nullptr);

    pointMap_.clear();
    faceMap_.clear();
    cellMap_.clear();
    patchMap_.clear();
}


void Foam::faMeshSubset::reset()
{
    clear();
}


void Foam::faMeshSubset::reset(const Foam::zero)
{
    clear();

    // Create zero-sized subMesh
    subMeshPtr_.reset
    (
        new faMesh(baseMesh_, Foam::zero{})
    );
    auto& newSubMesh = subMeshPtr_();


    // Clone non-processor patches
    {
        const faBoundaryMesh& oldBoundary = baseMesh_.boundary();
        const faBoundaryMesh& newBoundary = newSubMesh.boundary();

        faPatchList newPatches(oldBoundary.nNonProcessor());

        patchMap_ = identity(newPatches.size());

        forAll(newPatches, patchi)
        {
            newPatches.set
            (
                patchi,
                oldBoundary[patchi].clone
                (
                    newBoundary,
                    labelList(),   // edgeLabels
                    patchi,
                    oldBoundary[patchi].ngbPolyPatchIndex()
                )
            );
        }

        newSubMesh.addFaPatches(newPatches);
    }
}


// ************************************************************************* //
