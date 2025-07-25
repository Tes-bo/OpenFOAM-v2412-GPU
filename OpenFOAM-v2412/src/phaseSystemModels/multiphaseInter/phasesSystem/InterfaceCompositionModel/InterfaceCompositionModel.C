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

#include "InterfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "rhoThermo.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::multiphaseInter;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::multiComponentMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& globalThermo
) const
{
    return
        globalThermo.getLocalThermo
        (
            globalThermo.species().find(speciesName)
        );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::pureMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const pureMixture<ThermoType>& globalThermo
) const
{
    return globalThermo.cellMixture(0);
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getSpecieMassFraction
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    auto tY = volScalarField::New
    (
        "tY",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    auto& Ys = tY.ref();

    Ys = mixture.Y(speciesName);

    return tY;
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getSpecieMassFraction
(
    const word& speciesName,
    const pureMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    return volScalarField::New
    (
        "tY",
        IOobject::NO_REGISTER,
        mesh,
        scalar(1),
        dimless,
        fvPatchFieldBase::zeroGradientType()
    );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::MwMixture
(
    const pureMixture<ThermoType>& mixture
) const
{
    const fvMesh& mesh = fromThermo_.p().mesh();

    return volScalarField::New
    (
        "tM",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar
        (
            "Mw",
            dimMass/dimMoles,
            1e-3*mixture.cellMixture(0).W()
        ),
        fvPatchFieldBase::zeroGradientType()
    );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::MwMixture
(
    const multiComponentMixture<ThermoType>& mixture
) const
{
    return refCast<const basicSpecieMixture>(mixture).W();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::InterfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceCompositionModel(dict, pair),
    fromThermo_
    (
        pair.from().mesh().lookupObject<Thermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                pair.from().name()
            )
        )
    ),
    toThermo_
    (
        pair.to().mesh().lookupObject<OtherThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                pair.to().name()
            )
        )
    ),
    Le_("Le", dimless, 1.0, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::Dto
(
    const word& speciesName
) const
{
    const typename OtherThermo::thermoType& toThermoType =
        getLocalThermo
        (
            speciesName,
            toThermo_
        );

    const volScalarField& p = toThermo_.p();

    const volScalarField& T = toThermo_.T();

    auto tmpD = volScalarField::New
    (
        IOobject::groupName("D", pair_.name()),
        IOobject::NO_REGISTER,
        p.mesh(),
        dimensionedScalar(dimArea/dimTime, Zero)
    );
    auto& D = tmpD.ref();

    forAll(p, celli)
    {
        D[celli] =
            toThermoType.alphah(p[celli], T[celli])
           /toThermoType.rho(p[celli], T[celli]);
    }

    D /= Le_;
    D.correctBoundaryConditions();

    return tmpD;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::Dfrom
(
    const word& speciesName
) const
{
    const typename Thermo::thermoType& fromThermoType =
        getLocalThermo
        (
            speciesName,
            fromThermo_
        );

    const volScalarField& p(fromThermo_.p());

    const volScalarField& T(fromThermo_.T());

    auto tmpD = volScalarField::New
    (
        IOobject::groupName("D", pair_.name()),
        IOobject::NO_REGISTER,
        p.mesh(),
        dimensionedScalar(dimArea/dimTime, Zero)
    );
    auto& D = tmpD.ref();

    forAll(p, celli)
    {
        D[celli] =
            fromThermoType.alphah(p[celli], T[celli])
           /fromThermoType.rho(p[celli], T[celli]);
    }

    D /= Le_;
    D.correctBoundaryConditions();

    return tmpD;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::L
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const typename Thermo::thermoType& fromThermo =
        getLocalThermo(speciesName, fromThermo_);
    const typename OtherThermo::thermoType& toThermo =
        getLocalThermo(speciesName, toThermo_);

    const volScalarField& p(fromThermo_.p());

    auto tmpL = volScalarField::New
    (
        IOobject::groupName("L", pair_.name()),
        IOobject::NO_REGISTER,
        p.mesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    auto& L = tmpL.ref();

    // from Thermo (from) to Thermo (to)
    forAll(p, celli)
    {
        L[celli] = fromThermo.Hc() - toThermo.Hc();
    }

    L.correctBoundaryConditions();

    return tmpL;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::dY
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    NotImplemented;
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    NotImplemented;
    return nullptr;
}


// ************************************************************************* //
