/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "Stokes.H"
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

template<class BasicTurbulenceModel>
Stokes<BasicTurbulenceModel>::Stokes
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
    linearViscousStress<laminarModel<BasicTurbulenceModel>>
    (
        typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
const dictionary&
Stokes<BasicTurbulenceModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicTurbulenceModel>
bool Stokes<BasicTurbulenceModel>::read()
{
    return true;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
Stokes<BasicTurbulenceModel>::nut() const
{
    return volScalarField::New
    (
        IOobject::groupName("nut", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicTurbulenceModel>
tmp<scalarField>
Stokes<BasicTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>::New(this->mesh_.boundary()[patchi].size(), Zero);
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
Stokes<BasicTurbulenceModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->nu()
    );
}


template<class BasicTurbulenceModel>
tmp<scalarField>
Stokes<BasicTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicTurbulenceModel>
void Stokes<BasicTurbulenceModel>::correct()
{
    laminarModel<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// ************************************************************************* //
