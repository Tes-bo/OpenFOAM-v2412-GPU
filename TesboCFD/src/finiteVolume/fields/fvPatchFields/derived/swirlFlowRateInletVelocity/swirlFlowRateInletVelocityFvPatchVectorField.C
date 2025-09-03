/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "swirlFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    origin_(),
    axis_(Zero),
    flowRate_(),
    rpm_()
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    origin_
    (
        dict.getOrDefault
        (
            "origin",
            returnReduceOr(patch().size())
          ? gSum(patch().Cf()*patch().magSf())/gSum(patch().magSf())
          : Zero
        )
    ),
    axis_
    (
        dict.getOrDefault
        (
            "axis",
            returnReduceOr(patch().size())
          ? -gSum(patch().Sf())/gSum(patch().magSf())
          : Zero
        )
    ),
    flowRate_(Function1<scalar>::New("flowRate", dict, &db())),
    rpm_(Function1<scalar>::New("rpm", dict, &db()))
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    flowRate_(ptf.flowRate_.clone()),
    rpm_(ptf.rpm_.clone())
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    flowRate_(ptf.flowRate_.clone()),
    rpm_(ptf.rpm_.clone())
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    flowRate_(ptf.flowRate_.clone()),
    rpm_(ptf.rpm_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const scalar totArea = gSum(patch().magSf());

    if (totArea > ROOTVSMALL && axis_ != vector(Zero))
    {
        const scalar t = this->db().time().timeOutputValue();
        const scalar flowRate = flowRate_->value(t);
        const scalar omega = rpmToRads(rpm_->value(t));

        const scalar avgU = -flowRate/totArea;

        const vector axisHat = axis_/mag(axis_);

        // Update angular velocity
        tmp<vectorField> tangentialVelocity
        (
            axisHat ^ omega*(patch().Cf() - origin_)
        );

        tmp<vectorField> n = patch().nf();

        const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);

        if (phi.dimensions() == dimVolume/dimTime)
        {
            // volumetric flow-rate
            operator==(tangentialVelocity + n*avgU);
        }
        else if (phi.dimensions() == dimMass/dimTime)
        {
            const auto& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            // mass flow-rate
            operator==(tangentialVelocity + n*avgU/rhop);
        }
        else
        {
            FatalErrorInFunction
                << "dimensions of " << phiName_ << " are incorrect" << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::swirlFlowRateInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    flowRate_->writeData(os);
    rpm_->writeData(os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
