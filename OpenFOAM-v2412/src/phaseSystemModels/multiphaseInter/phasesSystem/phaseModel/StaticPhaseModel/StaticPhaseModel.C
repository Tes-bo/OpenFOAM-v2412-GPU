/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "StaticPhaseModel.H"

#include "multiphaseInterSystem.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StaticPhaseModel<BasePhaseModel>::StaticPhaseModel
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    U_(fluid.mesh().lookupObject<volVectorField>("U")),
    phi_
    (
        IOobject
        (
            IOobject::groupName("phi", multiphaseInter::phaseModel::name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName
            (
                "alphaPhi",
                multiphaseInter::phaseModel::name()
            ),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::StaticPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::phi() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("phi", multiphaseInter::phaseModel::name()),
        IOobject::NO_REGISTER,
        U_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    );
}


template<class BasePhaseModel>
const Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::phi()
{
    phi_ = dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero);
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "alphaPhi",
            multiphaseInter::phaseModel::name()
        ),
        IOobject::NO_REGISTER,
        U_.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    );
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhi()
{
    alphaPhi_ = dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero);
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StaticPhaseModel<BasePhaseModel>::U() const
{
    return volVectorField::New
    (
        IOobject::groupName("U", multiphaseInter::phaseModel::name()),
        IOobject::NO_REGISTER,
        U_.mesh(),
        dimensionedVector(dimVelocity, Zero)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField> Foam::StaticPhaseModel<BasePhaseModel>
::diffNo() const
{
    tmp<surfaceScalarField> tkapparhoCpbyDelta
    (
        sqr(U_.mesh().surfaceInterpolation::deltaCoeffs())
       *fvc::interpolate(this->kappa().ref())
       /fvc::interpolate((this->Cp()*this->rho())())
    );

    return tkapparhoCpbyDelta;
}


// ************************************************************************* //
