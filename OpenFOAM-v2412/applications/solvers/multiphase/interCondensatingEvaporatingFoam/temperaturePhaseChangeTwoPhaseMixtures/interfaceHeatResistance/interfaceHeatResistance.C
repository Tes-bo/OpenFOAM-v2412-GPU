/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
    Copyright (C) 2020 Henning Scheufler
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

#include "interfaceHeatResistance.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseMixtureEThermo.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "calculatedFvPatchFields.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(interfaceHeatResistance, 0);
    addToRunTimeSelectionTable
    (
        temperaturePhaseChangeTwoPhaseMixture,
        interfaceHeatResistance,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
interfaceHeatResistance
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    temperaturePhaseChangeTwoPhaseMixture(mixture, mesh),
    R_
    (
        "R",
        dimPower/dimArea/dimTemperature, optionalSubDict(type() + "Coeffs")
    ),

    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    mDotc_
    (
        IOobject
        (
            "mDotc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDote_
    (
        IOobject
        (
            "mDote",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDotcSpread_
    (
        IOobject
        (
            "mDotcSpread",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDoteSpread_
    (
        IOobject
        (
            "mDoteSpread",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    spread_
    (
        optionalSubDict(type() + "Coeffs").get<scalar>("spread")
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
vDotAlphal() const
{
    dimensionedScalar alphalCoeff(1.0/mixture_.rho1());

    return Pair<tmp<volScalarField>>
    (
        (alphalCoeff*mDotc_)/(mixture_.alpha2() + SMALL),
       -(alphalCoeff*mDote_)/(mixture_.alpha1() + SMALL)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
mDotAlphal() const
{
    return Pair<tmp<volScalarField>>
    (
        (mDotc_/clamp(mixture_.alpha2(), scalarMinMax(SMALL, 1))),
       -(mDote_/clamp(mixture_.alpha1(), scalarMinMax(SMALL, 1)))
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
mDot() const
{
    return Pair<tmp<volScalarField>>
    (
        tmp<volScalarField>(mDotc_),
        tmp<volScalarField>(mDote_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
mDotDeltaT() const
{
   const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const dimensionedScalar& TSat = thermo.TSat();

    Pair<tmp<volScalarField>> mDotce(mDot());

    return Pair<tmp<volScalarField>>
    (
        mDotc_*pos(TSat - T.oldTime())/(TSat - T.oldTime()),
       -mDote_*pos(T.oldTime() - TSat)/(T.oldTime() - TSat)
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
TSource() const
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

    // interface heat resistance
    volScalarField IHRcoeff(interfaceArea_*R_);

    TSource = fvm::Sp(IHRcoeff, T) - IHRcoeff*TSat;

    return tTSource;
}


void Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
correct()
{
    // Update Interface
    updateInterface();

    // Update mDotc_ and mDote_
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();
    const dimensionedScalar T0(dimTemperature, Zero);

    const dimensionedScalar L(mag(mixture_.Hf2() - mixture_.Hf1()));

    // interface heat resistance
    mDotc_ = interfaceArea_*R_*max(TSat - T, T0)/L;
    mDote_ = interfaceArea_*R_*max(T - TSat, T0)/L;

    // Calculate the spread sources
    dimensionedScalar D
    (
        "D",
        dimArea,
        spread_/sqr(gAverage(mesh_.nonOrthDeltaCoeffs()))
    );


    const volScalarField& alpha1 = mixture_.alpha1();
    const volScalarField& alpha2 = mixture_.alpha2();

    const dimensionedScalar MDotMin("MdotMin", mDotc_.dimensions(), 1e-3);

    if (max(mDotc_) > MDotMin)
    {
        fvc::spreadSource
        (
            mDotcSpread_,
            mDotc_,
            alpha1,
            alpha2,
            D,
            1e-3
        );
    }

    if (max(mDote_) > MDotMin)
    {
        fvc::spreadSource
        (
            mDoteSpread_,
            mDote_,
            alpha1,
            alpha2,
            D,
            1e-3
        );
    }
}


void Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
updateInterface()
{

    // interface heat resistance
    // Interpolating alpha1 cell centre values to mesh points (vertices)
    scalarField ap
    (
        volPointInterpolation::New(mesh_).interpolate(mixture_.alpha1())
    );

    cutCellIso cutCell(mesh_, ap);

    forAll(interfaceArea_, celli)
    {
        label status = cutCell.calcSubCell(celli, 0.5);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] =
                mag(cutCell.faceArea())/mesh_.V()[celli];
        }
    }
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
vDot() const
{

    dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());

    return Pair<tmp<volScalarField>>
    (
        pCoeff*mDotcSpread_,
       -pCoeff*mDoteSpread_
    );
}


bool Foam::temperaturePhaseChangeTwoPhaseMixtures::interfaceHeatResistance::
read()
{
    if (temperaturePhaseChangeTwoPhaseMixture::read())
    {
        optionalSubDict(type() + "Coeffs").readEntry("R", R_);
        optionalSubDict(type() + "Coeffs").readEntry("spread", spread_);

        return true;
    }

    return false;
}


// ************************************************************************* //
