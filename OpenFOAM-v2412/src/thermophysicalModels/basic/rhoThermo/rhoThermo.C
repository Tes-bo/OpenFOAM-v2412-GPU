/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoThermo, 0);
    defineRunTimeSelectionTable(rhoThermo, fvMesh);
    defineRunTimeSelectionTable(rhoThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoThermo::rhoThermo(const fvMesh& mesh, const word& phaseName)
:
    fluidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
            phaseScopedName("thermo", "rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phaseScopedName("thermo", "psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phaseScopedName("thermo", "mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


Foam::rhoThermo::rhoThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phaseScopedName("thermo", "rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phaseScopedName("thermo", "psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phaseScopedName("thermo", "mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


Foam::rhoThermo::rhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    fluidThermo(mesh, phaseName, dictionaryName),
    rho_
    (
        IOobject
        (
            phaseScopedName("thermo", "rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phaseScopedName("thermo", "psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phaseScopedName("thermo", "mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoThermo> Foam::rhoThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::rhoThermo> Foam::rhoThermo::New
(
     const fvMesh& mesh,
     const word& phaseName,
     const word& dictName
)
{
    return basicThermo::New<rhoThermo>(mesh, phaseName, dictName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoThermo::~rhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rhoThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::rhoThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rhoThermo::rho()
{
    return rho_;
}


void Foam::rhoThermo::correctRho
(
    const Foam::volScalarField& deltaRho,
    const dimensionedScalar& rhoMin,
    const dimensionedScalar& rhoMax
)
{
    rho_ += deltaRho;
    rho_ = max(rho_, rhoMin);
    rho_ = min(rho_, rhoMax);
}

void Foam::rhoThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    rho_ += deltaRho;
}


const Foam::volScalarField& Foam::rhoThermo::psi() const
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::rhoThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::rhoThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


// ************************************************************************* //
