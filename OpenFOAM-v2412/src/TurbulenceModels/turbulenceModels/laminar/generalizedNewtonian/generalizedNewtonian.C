/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "generalizedNewtonian.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
generalizedNewtonian<BasicMomentumTransportModel>::generalizedNewtonian
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    linearViscousStress<laminarModel<BasicMomentumTransportModel>>
    (
        typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    viscosityModel_
    (
        generalizedNewtonianViscosityModel::New
        (
            this->coeffDict_
        )
    ),

    nu_
    (
        IOobject
        (
            IOobject::groupName("generalizedNewtonian:nu", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        viscosityModel_->nu(this->nu(), strainRate())
    )
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
generalizedNewtonian<BasicMomentumTransportModel>::strainRate() const
{
    return sqrt(scalar{2})*mag(symm(fvc::grad(this->U())));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool generalizedNewtonian<BasicMomentumTransportModel>::read()
{
    viscosityModel_->read(this->coeffDict_);

    return true;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
generalizedNewtonian<BasicMomentumTransportModel>::nut() const
{
    return volScalarField::New
    (
        IOobject::groupName("nut", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicMomentumTransportModel>
tmp<scalarField>
generalizedNewtonian<BasicMomentumTransportModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>::New(this->mesh_.boundary()[patchi].size(), Zero);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
generalizedNewtonian<BasicMomentumTransportModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        nu_
    );
}


template<class BasicMomentumTransportModel>
tmp<scalarField>
generalizedNewtonian<BasicMomentumTransportModel>::nuEff
(
    const label patchi
) const
{
    return nu_.boundaryField()[patchi];
}


template<class BasicMomentumTransportModel>
void generalizedNewtonian<BasicMomentumTransportModel>::correct()
{
    nu_ = viscosityModel_->nu(this->nu(), strainRate());
    laminarModel<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// ************************************************************************* //
