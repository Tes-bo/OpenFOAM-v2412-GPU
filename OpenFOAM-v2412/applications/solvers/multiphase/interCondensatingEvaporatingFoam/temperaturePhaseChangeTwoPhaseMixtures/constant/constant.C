/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "twoPhaseMixtureEThermo.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        temperaturePhaseChangeTwoPhaseMixture,
        constant,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::constant
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    temperaturePhaseChangeTwoPhaseMixture(mixture, mesh),
    coeffC_
    (
        "coeffC",
        dimless/dimTime/dimTemperature,
        optionalSubDict(type() + "Coeffs")
    ),
    coeffE_
    (
        "coeffE",
        dimless/dimTime/dimTemperature,
        optionalSubDict(type() + "Coeffs")
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotAlphal() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*max(TSat - T, T0),
       -coeffE_*mixture_.rho1()*max(T - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDot() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField mDotE
    (
        "mDotE",
        coeffE_*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
      * max(T - TSat, T0)
    );
    volScalarField mDotC
    (
        "mDotC",
        coeffC_*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
      * max(TSat - T, T0)
    );

    if (mesh_.time().writeTime())
    {
        mDotC.write();
        mDotE.write();
    }

    return Pair<tmp<volScalarField>>
    (
        tmp<volScalarField>(new volScalarField(mDotC)),
        tmp<volScalarField>(new volScalarField(-mDotE))
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotDeltaT() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    return Pair<tmp<volScalarField>>
    (
        (
            coeffC_*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
          * pos(TSat - T)
        ),
        (
            coeffE_*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
          * pos(T - TSat)
        )
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::TSource() const
{

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    auto tTSource = tmp<fvScalarMatrix>::New(T, dimEnergy/dimTime);
    auto& TSource = tTSource.ref();

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar L = mixture_.Hf2() - mixture_.Hf1();

    const volScalarField Vcoeff
    (
        coeffE_*mixture_.rho1()*clamp(mixture_.alpha1(), zero_one{})
      * L*pos(T - TSat)
    );
    const volScalarField Ccoeff
    (
        coeffC_*mixture_.rho2()*clamp(mixture_.alpha2(), zero_one{})
      * L*pos(TSat - T)
    );

    TSource =
        fvm::Sp(Vcoeff, T) - Vcoeff*TSat
      + fvm::Sp(Ccoeff, T) - Ccoeff*TSat;

    return tTSource;
}


void Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::correct()
{
}


bool Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::read()
{
    if (temperaturePhaseChangeTwoPhaseMixture::read())
    {
        subDict(type() + "Coeffs").readEntry("coeffC", coeffC_);
        subDict(type() + "Coeffs").readEntry("coeffE", coeffE_);

        return true;
    }

    return false;
}


// ************************************************************************* //
