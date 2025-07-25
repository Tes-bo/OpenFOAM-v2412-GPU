/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "liquidFilmThermo.H"
#include "demandDrivenData.H"
#include "thermoSingleLayer.H"
#include "SLGThermo.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmThermo, 0);

addToRunTimeSelectionTable
(
    filmThermoModel,
    liquidFilmThermo,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const thermoSingleLayer& liquidFilmThermo::thermoFilm() const
{
    if (!isA<thermoSingleLayer>(filmModel_))
    {
        FatalErrorInFunction
            << "Thermo model requires a " << thermoSingleLayer::typeName
            << " film to supply the pressure and temperature, but "
            << filmModel_.type() << " film model selected.  "
            << "Use the 'useReferenceValues' flag to employ reference "
            << "pressure and temperature" << exit(FatalError);
    }

    return refCast<const thermoSingleLayer>(filmModel_);
}


void liquidFilmThermo::initLiquid(const dictionary& dict)
{
    if (liquidPtr_ != nullptr)
    {
        return;
    }

    dict.readEntry("liquid", name_);

    const SLGThermo* thermoPtr =
        filmModel_.primaryMesh().findObject<SLGThermo>("SLGThermo");

    if (thermoPtr)
    {
        // Retrieve from film thermo
        ownLiquid_ = false;

        const SLGThermo& thermo = *thermoPtr;

        const label id = thermo.liquidId(name_);

        liquidPtr_ = &thermo.liquids().properties()[id];
    }
    else
    {
        // New liquid create
        ownLiquid_ = true;

        liquidPtr_ =
            liquidProperties::New(dict.optionalSubDict(name_ + "Coeffs")).ptr();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmThermo::liquidFilmThermo
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmThermoModel(typeName, film, dict),
    name_("unknown_liquid"),
    liquidPtr_(nullptr),
    ownLiquid_(false),
    useReferenceValues_(coeffDict_.get<bool>("useReferenceValues")),
    pRef_(0.0),
    TRef_(0.0)
{
    initLiquid(coeffDict_);

    if (useReferenceValues_)
    {
        coeffDict_.readEntry("pRef", pRef_);
        coeffDict_.readEntry("TRef", TRef_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmThermo::~liquidFilmThermo()
{
    if (ownLiquid_)
    {
        deleteDemandDrivenData(liquidPtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const word& liquidFilmThermo::name() const
{
    return name_;
}


scalar liquidFilmThermo::rho
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->rho(p, T);
}


scalar liquidFilmThermo::mu
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->mu(p, T);
}


scalar liquidFilmThermo::sigma
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->sigma(p, T);
}


scalar liquidFilmThermo::Cp
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->Cp(p, T);
}


scalar liquidFilmThermo::kappa
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->kappa(p, T);
}


scalar liquidFilmThermo::D
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->D(p, T);
}


scalar liquidFilmThermo::hl
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->hl(p, T);
}


scalar liquidFilmThermo::pv
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->pv(p, T);
}


scalar liquidFilmThermo::W() const
{
    return liquidPtr_->W();
}


scalar liquidFilmThermo::Tb(const scalar p) const
{
    return liquidPtr_->pvInvert(p);
}


tmp<volScalarField> liquidFilmThermo::rho() const
{
    auto trho = volScalarField::New
    (
        IOobject::scopedName(type(), "rho"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimDensity, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    scalarField& rho = trho.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        rho = this->rho(pRef_, TRef_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(rho, celli)
        {
            rho[celli] = this->rho(p[celli], T[celli]);
        }
    }

    trho.ref().correctBoundaryConditions();

    return trho;
}


tmp<volScalarField> liquidFilmThermo::mu() const
{
    auto tmu = volScalarField::New
    (
        IOobject::scopedName(type(), "mu"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimPressure*dimTime, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );
    scalarField& mu = tmu.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        mu = this->mu(pRef_, TRef_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(mu, celli)
        {
            mu[celli] = this->mu(p[celli], T[celli]);
        }
    }

    tmu.ref().correctBoundaryConditions();

    return tmu;
}


tmp<volScalarField> liquidFilmThermo::sigma() const
{
    auto tsigma = volScalarField::New
    (
        IOobject::scopedName(type(), "sigma"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );
    scalarField& sigma = tsigma.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        sigma = this->sigma(pRef_, TRef_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(sigma, celli)
        {
            sigma[celli] = this->sigma(p[celli], T[celli]);
        }
    }

    tsigma.ref().correctBoundaryConditions();

    return tsigma;
}


tmp<volScalarField> liquidFilmThermo::Cp() const
{
    auto tCp = volScalarField::New
    (
        IOobject::scopedName(type(), "Cp"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );
    scalarField& Cp = tCp.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        Cp = this->Cp(pRef_, TRef_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(Cp, celli)
        {
            Cp[celli] = this->Cp(p[celli], T[celli]);
        }
    }

    tCp.ref().correctBoundaryConditions();

    return tCp;
}


tmp<volScalarField> liquidFilmThermo::kappa() const
{
    auto tkappa = volScalarField::New
    (
        IOobject::scopedName(type(), "kappa"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );
    scalarField& kappa = tkappa.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        kappa = this->kappa(pRef_, TRef_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(kappa, celli)
        {
            kappa[celli] = this->kappa(p[celli], T[celli]);
        }
    }

    tkappa.ref().correctBoundaryConditions();

    return tkappa;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
