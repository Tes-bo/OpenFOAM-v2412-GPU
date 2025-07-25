/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(nullptr),
    rhoName_("rho"),
    rhoInlet_(0),
    volumetric_(false),
    extrapolateProfile_(false)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    flowRate_(nullptr),
    rhoName_("rho"),
    rhoInlet_(dict.getOrDefault<scalar>("rhoInlet", -VGREAT)),
    volumetric_(false),
    extrapolateProfile_
    (
        dict.getOrDefault<Switch>("extrapolateProfile", false)
    )
{
    flowRate_ =
        Function1<scalar>::NewIfPresent("volumetricFlowRate", dict, &db());

    if (flowRate_)
    {
        volumetric_ = true;
    }
    else
    {
        dict.readIfPresent("rho", rhoName_);
        flowRate_ =
            Function1<scalar>::NewIfPresent("massFlowRate", dict, &db());
    }

    if (!flowRate_)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' (optional: with 'rho')" << nl
            << exit(FatalIOError);
    }

    // Value field required if mass based
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::buffered);
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::flowRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    const vectorField n(patch().nf());

    if (extrapolateProfile_)
    {
        vectorField Up(this->patchInternalField());

        // Patch normal extrapolated velocity
        scalarField nUp(n & Up);

        // Remove the normal component of the extrapolate patch velocity
        Up -= nUp*n;

        // Remove any reverse flow
        nUp = min(nUp, scalar(0));

        const scalar flowRate = flowRate_->value(t);
        const scalar estimatedFlowRate = -gSum(rho*(this->patch().magSf()*nUp));

        if (estimatedFlowRate > 0.5*flowRate)
        {
            nUp *= (mag(flowRate)/mag(estimatedFlowRate));
        }
        else
        {
            nUp -= ((flowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
        }

        // Add the corrected normal component of velocity to the patch velocity
        Up += nUp*n;

        // Correct the patch velocity
        this->operator==(Up);
    }
    else
    {
        const scalar avgU = -flowRate_->value(t)/gSum(rho*patch().magSf());
        operator==(avgU*n);
    }
}


void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one{});
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const auto& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    }
    if (extrapolateProfile_)
    {
        os.writeEntry("extrapolateProfile", extrapolateProfile_);
    }
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
