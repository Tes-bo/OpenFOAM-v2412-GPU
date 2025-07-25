/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "variableHeightFlowRateFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    lowerBound_(0.0),
    upperBound_(1.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


Foam::variableHeightFlowRateFvPatchScalarField
::variableHeightFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    lowerBound_(dict.get<scalar>("lowerBound")),
    upperBound_(dict.get<scalar>("upperBound"))
{
    fvPatchFieldBase::readDict(dict);

    if (!this->readValueEntry(dict))
    {
        // Fallback: set to the internal field
        this->extrapolateInternal();
    }

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::variableHeightFlowRateFvPatchScalarField
    ::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


Foam::variableHeightFlowRateFvPatchScalarField
    ::variableHeightFlowRateFvPatchScalarField
(
    const variableHeightFlowRateFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    lowerBound_(ptf.lowerBound_),
    upperBound_(ptf.upperBound_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::variableHeightFlowRateFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);

    scalarField alphap(this->patchInternalField());


    forAll(phip, i)
    {
        if (phip[i] < -SMALL)
        {
            if (alphap[i] < lowerBound_)
            {
                this->refValue()[i] = 0.0;
            }
            else if (alphap[i] > upperBound_)
            {
                this->refValue()[i] = 1.0;
            }
            else
            {
                this->refValue()[i] = alphap[i];
            }

            this->valueFraction()[i] = 1.0;
        }
        else
        {
            this->refValue()[i] = 0.0;
            this->valueFraction()[i] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::variableHeightFlowRateFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("lowerBound", lowerBound_);
    os.writeEntry("upperBound", upperBound_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        variableHeightFlowRateFvPatchScalarField
    );
}

// ************************************************************************* //
