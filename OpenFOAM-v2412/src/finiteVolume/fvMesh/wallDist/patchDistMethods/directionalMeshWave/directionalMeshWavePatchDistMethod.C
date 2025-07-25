/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020,2023 OpenCFD Ltd.
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

#include "patchDataWave.H"
#include "directionalMeshWavePatchDistMethod.H"
#include "fvMesh.H"
#include "volFields.H"
#include "directionalWallPointData.H"
#include "emptyFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(directionalMeshWave, 0);
    addToRunTimeSelectionTable
    (
        patchDistMethod,
        directionalMeshWave,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::directionalMeshWave::directionalMeshWave
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    Foam::patchDistMethods::meshWave(dict, mesh, patchIDs),
    n_(dict.get<vector>("normal"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::directionalMeshWave::correct(volScalarField& y)
{
    y = dimensionedScalar(dimLength, GREAT);

    volVectorField n
    (
        IOobject
        (
            "nWall",
            y.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        patchDistMethod::patchTypes<scalar>(mesh_, patchIDs_)
    );

    const fvPatchList& patches = mesh_.boundary();

    volVectorField::Boundary& nbf = n.boundaryFieldRef();

    for (const label patchi : patchIDs_)
    {
        nbf[patchi] == patches[patchi].nf();
    }

    return correct(y, n);
}


bool Foam::patchDistMethods::directionalMeshWave::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    y = dimensionedScalar(dimLength, GREAT);

    // Collect pointers to data on patches
    UPtrList<vectorField> patchData(mesh_.boundaryMesh().size());

    volVectorField::Boundary& nbf = n.boundaryFieldRef();

    forAll(nbf, patchi)
    {
        patchData.set(patchi, &nbf[patchi]);
    }

    // Do mesh wave
    vector testDirection(n_);

    patchDataWave<directionalWallPointData<vector>, vector> wave
    (
        mesh_,
        patchIDs_,
        patchData,
        correctWalls_,
        testDirection
    );

    // Transfer cell values from wave into y and n
    y.transfer(wave.distance());

    n.transfer(wave.cellData());

    // Transfer values on patches into boundaryField of y and n
    volScalarField::Boundary& ybf = y.boundaryFieldRef();

    forAll(ybf, patchi)
    {
        scalarField& waveFld = wave.patchDistance()[patchi];

        if (!isA<emptyFvPatchScalarField>(ybf[patchi]))
        {
            ybf[patchi].transfer(waveFld);

            vectorField& wavePatchData = wave.patchData()[patchi];

            nbf[patchi].transfer(wavePatchData);
        }
    }

    // Make sure boundary values are up-to-date
    y.correctBoundaryConditions();
    n.correctBoundaryConditions();

    // Transfer number of unset values
    this->nUnset_ = wave.nUnset();

    return this->nUnset_ > 0;
}


// ************************************************************************* //
