/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "phasePressureModel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::phasePressureModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    alphaMax_(coeffDict_.get<scalar>("alphaMax")),
    preAlphaExp_(coeffDict_.get<scalar>("preAlphaExp")),
    expMax_(coeffDict_.get<scalar>("expMax")),
    g0_("g0", dimPressure, coeffDict_)
{
    nut_ == dimensionedScalar(nut_.dimensions());

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::phasePressureModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        coeffDict().readEntry("alphaMax", alphaMax_);
        coeffDict().readEntry("preAlphaExp", preAlphaExp_);
        coeffDict().readEntry("expMax", expMax_);
        g0_.readIfPresent(coeffDict());

        return true;
    }

    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::omega() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::R() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("R", U_.group()),
        IOobject::NO_REGISTER,
        mesh_,
        dimensioned<symmTensor>(dimensionSet(0, 2, -2, 0, 0), Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(alpha_ - alphaMax_)),
            expMax_
        )
    );

    auto& bpPrime = tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::phasePressureModel::pPrimef() const
{
    tmp<surfaceScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(fvc::interpolate(alpha_) - alphaMax_)),
            expMax_
        )
    );

    auto& bpPrime = tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::devRhoReff() const
{
    return devRhoReff(U_);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::devRhoReff
(
    const volVectorField& U
) const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("devRhoReff", U.group()),
        IOobject::NO_REGISTER,
        mesh_,
        dimensioned<symmTensor>
        (
            rho_.dimensions()*dimensionSet(0, 2, -2, 0, 0),
            Zero
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::phasePressureModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>::New
    (
        U,
        rho_.dimensions()*dimensionSet(0, 4, -2, 0, 0)
    );
}


void Foam::RASModels::phasePressureModel::correct()
{}


// ************************************************************************* //
