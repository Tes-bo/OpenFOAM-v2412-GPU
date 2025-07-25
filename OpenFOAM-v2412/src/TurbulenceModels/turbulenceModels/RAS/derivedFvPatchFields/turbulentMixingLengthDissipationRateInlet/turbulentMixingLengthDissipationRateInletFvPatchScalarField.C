/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "turbulentMixingLengthDissipationRateInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    kName_("k"),
    mixingLength_(0.0),
    Cmu_(0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    mixingLength_(ptf.mixingLength_),
    Cmu_(ptf.Cmu_)
{}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    kName_(dict.getOrDefault<word>("k", "k")),
    mixingLength_
    (
        dict.getCheck<scalar>("mixingLength", scalarMinMax::ge(SMALL))
    ),
    Cmu_(dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL)))
{
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    kName_(ptf.kName_),
    mixingLength_(ptf.mixingLength_),
    Cmu_(ptf.Cmu_)
{}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    kName_(ptf.kName_),
    mixingLength_(ptf.mixingLength_),
    Cmu_(ptf.Cmu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMixingLengthDissipationRateInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    Cmu_ = turbModel.coeffDict().getOrDefault<scalar>("Cmu", Cmu_);

    const scalar Cmu75 = pow(Cmu_, 0.75);

    const auto& kp =
        patch().lookupPatchField<volScalarField>(kName_);

    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    this->refValue() = (Cmu75/mixingLength_)*pow(kp, 1.5);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void turbulentMixingLengthDissipationRateInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("mixingLength", mixingLength_);
    os.writeEntry("phi", this->phiName_);
    os.writeEntry("k", kName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMixingLengthDissipationRateInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
