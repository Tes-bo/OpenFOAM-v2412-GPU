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

#include "phaseModel.H"
#include "multiphaseInterSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseInter
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, multiphaseInterSystem);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInter::phaseModel::phaseModel
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    fluid_(fluid),
    name_(phaseName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiphaseInter::phaseModel>
Foam::multiphaseInter::phaseModel::New
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
{
    const dictionary& dict = fluid.subDict(phaseName);

    const word modelType(dict.get<word>("type"));

    Info<< "Selecting phaseModel for "
        << phaseName << ": " << modelType << endl;

    auto* ctorPtr = multiphaseInterSystemConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "phaseModel",
            modelType,
            *multiphaseInterSystemConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(fluid, phaseName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::multiphaseInterSystem&
Foam::multiphaseInter::phaseModel::fluid() const
{
    return fluid_;
}


void Foam::multiphaseInter::phaseModel::correct()
{
    thermo().correct();
}


void Foam::multiphaseInter::phaseModel::correctTurbulence()
{
    // do nothing
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::rho() const
{
    return thermo().rho();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInter::phaseModel::rho(const label patchI) const
{
     return thermo().rho(patchI);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::hc() const
{
     return thermo().hc();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::Cp() const
{
     return thermo().Cp();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return (thermo().Cp(p, T, patchI));
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::Cv() const
{
    return thermo().Cv();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().Cv(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::gamma() const
{
    return thermo().gamma();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().gamma(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::Cpv() const
{
    return thermo().Cpv();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().Cpv(p, T, patchI);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::CpByCpv() const
{
     return thermo().CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return thermo().CpByCpv(p, T, patchI);
}


const Foam::volScalarField& Foam::multiphaseInter::phaseModel::alpha() const
{
    return thermo().alpha();
}


const Foam::scalarField&
Foam::multiphaseInter::phaseModel::alpha(const label patchI) const
{
    return thermo().alpha(patchI);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::kappa() const
{
    return thermo().kappa();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInter::phaseModel::kappa(const label patchI) const
{
    return thermo().kappa(patchI);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseInter::phaseModel::alphahe() const
{
    return thermo().alphahe();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInter::phaseModel::alphahe(const label patchI) const
{
    return thermo().alphahe(patchI);
}


Foam::tmp<Foam::volScalarField>Foam::multiphaseInter::phaseModel::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff" + name_);
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::kappaEff
(
    const scalarField& kappat,
    const label patchI
) const
{
    return (kappa(patchI) + kappat);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return (thermo().alpha() + alphat);
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInter::phaseModel::alphaEff
(
    const scalarField& alphat,
    const label patchI
) const
{
    return (thermo().alpha(patchI) + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::mu() const
{
    return thermo().mu();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInter::phaseModel::mu(const label patchi) const
{
    return thermo().mu(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInter::phaseModel::nu() const
{
    return thermo().nu();
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInter::phaseModel::nu(const label patchi) const
{
    return thermo().nu(patchi);
}


bool Foam::multiphaseInter::phaseModel::read()
{
    return true;
}


// ************************************************************************* //
