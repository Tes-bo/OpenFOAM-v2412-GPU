/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020,2024 OpenCFD Ltd.
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

#include "plenumPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    gamma_(1.4),
    R_(287.04),
    supplyMassFlowRate_(1.0),
    supplyTotalTemperature_(300.0),
    plenumVolume_(1.0),
    plenumDensity_(1.0),
    plenumDensityOld_(1.0),
    plenumTemperature_(300.0),
    plenumTemperatureOld_(300.0),
    rho_(1.0),
    hasRho_(false),
    inletAreaRatio_(1.0),
    inletDischargeCoefficient_(1.0),
    timeScale_(0.0),
    timeIndex_(-1),
    phiName_("phi"),
    UName_("U")
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    gamma_(dict.get<scalar>("gamma")),
    R_(dict.get<scalar>("R")),
    supplyMassFlowRate_(dict.get<scalar>("supplyMassFlowRate")),
    supplyTotalTemperature_(dict.get<scalar>("supplyTotalTemperature")),
    plenumVolume_(dict.get<scalar>("plenumVolume")),
    plenumDensity_(dict.get<scalar>("plenumDensity")),
    plenumDensityOld_(1.0),
    plenumTemperature_(dict.get<scalar>("plenumTemperature")),
    plenumTemperatureOld_(300.0),
    rho_(1.0),
    hasRho_(false),
    inletAreaRatio_(dict.get<scalar>("inletAreaRatio")),
    inletDischargeCoefficient_(dict.get<scalar>("inletDischargeCoefficient")),
    timeScale_(dict.getOrDefault<scalar>("timeScale", 0)),
    timeIndex_(-1),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    UName_(dict.getOrDefault<word>("U", "U"))
{
    hasRho_ = dict.readIfPresent("rho", rho_);
}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supplyMassFlowRate_(ptf.supplyMassFlowRate_),
    supplyTotalTemperature_(ptf.supplyTotalTemperature_),
    plenumVolume_(ptf.plenumVolume_),
    plenumDensity_(ptf.plenumDensity_),
    plenumDensityOld_(ptf.plenumDensityOld_),
    plenumTemperature_(ptf.plenumTemperature_),
    plenumTemperatureOld_(ptf.plenumTemperatureOld_),
    rho_(ptf.rho_),
    hasRho_(ptf.hasRho_),
    inletAreaRatio_(ptf.inletAreaRatio_),
    inletDischargeCoefficient_(ptf.inletDischargeCoefficient_),
    timeScale_(ptf.timeScale_),
    timeIndex_(ptf.timeIndex_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    gamma_(tppsf.gamma_),
    R_(tppsf.R_),
    supplyMassFlowRate_(tppsf.supplyMassFlowRate_),
    supplyTotalTemperature_(tppsf.supplyTotalTemperature_),
    plenumVolume_(tppsf.plenumVolume_),
    plenumDensity_(tppsf.plenumDensity_),
    plenumDensityOld_(tppsf.plenumDensityOld_),
    plenumTemperature_(tppsf.plenumTemperature_),
    plenumTemperatureOld_(tppsf.plenumTemperatureOld_),
    rho_(tppsf.rho_),
    hasRho_(tppsf.hasRho_),
    inletAreaRatio_(tppsf.inletAreaRatio_),
    inletDischargeCoefficient_(tppsf.inletDischargeCoefficient_),
    timeScale_(tppsf.timeScale_),
    timeIndex_(tppsf.timeIndex_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_)
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    gamma_(tppsf.gamma_),
    R_(tppsf.R_),
    supplyMassFlowRate_(tppsf.supplyMassFlowRate_),
    supplyTotalTemperature_(tppsf.supplyTotalTemperature_),
    plenumVolume_(tppsf.plenumVolume_),
    plenumDensity_(tppsf.plenumDensity_),
    plenumDensityOld_(tppsf.plenumDensityOld_),
    plenumTemperature_(tppsf.plenumTemperature_),
    plenumTemperatureOld_(tppsf.plenumTemperatureOld_),
    rho_(tppsf.rho_),
    hasRho_(tppsf.hasRho_),
    inletAreaRatio_(tppsf.inletAreaRatio_),
    inletDischargeCoefficient_(tppsf.inletDischargeCoefficient_),
    timeScale_(tppsf.timeScale_),
    timeIndex_(tppsf.timeIndex_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plenumPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch properties
    const fvPatchField<scalar>& p = *this;
    const fvPatchField<scalar>& p_old =
        db().lookupObject<volScalarField>
        (
            internalField().name()
        ).oldTime().boundaryField()[patch().index()];

    const auto& U = patch().lookupPatchField<volVectorField>(UName_);
    const auto& phi = patch().lookupPatchField<surfaceScalarField>(phiName_);

    // Get the timestep
    const scalar dt = db().time().deltaTValue();

    // Check if operating at a new time index and update the old-time properties
    // if so
    if (timeIndex_ != db().time().timeIndex())
    {
        timeIndex_ = db().time().timeIndex();
        plenumDensityOld_ = plenumDensity_;
        plenumTemperatureOld_ = plenumTemperature_;
    }

    // Calculate the current mass flow rate
    scalar massFlowRate(1.0);
    if (phi.internalField().dimensions() == dimVolume/dimTime)
    {
        if (hasRho_)
        {
            massFlowRate = - gSum(rho_*phi);
        }
        else
        {
            FatalErrorInFunction
                << "The density must be specified when using a volumetric flux."
                << exit(FatalError);
        }
    }
    else if (phi.internalField().dimensions() == dimMass/dimTime)
    {
        if (hasRho_)
        {
            FatalErrorInFunction
                << "The density must be not specified when using a mass flux."
                << exit(FatalError);
        }
        else
        {
            massFlowRate = - gSum(phi);
        }
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    // Calculate the specific heats
    const scalar cv = R_/(gamma_ - 1), cp = R_*gamma_/(gamma_ - 1);

    // Calculate the new plenum properties
    plenumDensity_ =
        plenumDensityOld_
      + (dt/plenumVolume_)*(supplyMassFlowRate_ - massFlowRate);
    plenumTemperature_ =
        plenumTemperatureOld_
      + (dt/(plenumDensity_*cv*plenumVolume_))
      * (
            supplyMassFlowRate_
           *(cp*supplyTotalTemperature_ - cv*plenumTemperature_)
          - massFlowRate*R_*plenumTemperature_
        );
    const scalar plenumPressure = plenumDensity_*R_*plenumTemperature_;

    // Squared velocity magnitude at exit of channels
    const scalarField U_e(magSqr(U/inletAreaRatio_));

    // Exit temperature to plenum temperature ratio
    const scalarField r
    (
        1.0 - (gamma_ - 1.0)*U_e/(2.0*gamma_*R_*plenumTemperature_)
    );

    // Quadratic coefficient (others not needed as b = +1.0 and c = -1.0)
    const scalarField a
    (
        (1.0 - r)/(r*r*inletDischargeCoefficient_*inletDischargeCoefficient_)
    );

    // Isentropic exit temperature to plenum temperature ratio
    const scalarField s(2.0/(1.0 + sqrt(1.0 + 4.0*a)));

    // Exit pressure to plenum pressure ratio
    const scalarField t(pow(s, gamma_/(gamma_ - 1.0)));

    // Limit to prevent outflow
    const scalarField p_new
    (
        lerp
        (
            t*plenumPressure,           // Negative phi
            max(p, plenumPressure),     // Positive phi
            pos0(phi)                   // 0-1 selector
        )
    );

    // Relaxation fraction
    const scalar oneByFraction = timeScale_/dt;
    const scalar fraction = oneByFraction < 1.0 ? 1.0 : 1.0/oneByFraction;

    // Set the new value
    operator==(lerp(p_old, p_new, fraction));
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::plenumPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("gamma", gamma_);
    os.writeEntry("R", R_);
    os.writeEntry("supplyMassFlowRate", supplyMassFlowRate_);
    os.writeEntry("supplyTotalTemperature", supplyTotalTemperature_);
    os.writeEntry("plenumVolume", plenumVolume_);
    os.writeEntry("plenumDensity", plenumDensity_);
    os.writeEntry("plenumTemperature", plenumTemperature_);
    if (hasRho_)
    {
        os.writeEntry("rho", rho_);
    }
    os.writeEntry("inletAreaRatio", inletAreaRatio_);
    os.writeEntry("inletDischargeCoefficient", inletDischargeCoefficient_);
    os.writeEntryIfDifferent<scalar>("timeScale", 0.0, timeScale_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        plenumPressureFvPatchScalarField
    );
}

// ************************************************************************* //
