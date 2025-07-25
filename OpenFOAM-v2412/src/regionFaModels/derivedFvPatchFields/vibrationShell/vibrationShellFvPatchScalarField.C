/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "vibrationShellFvPatchScalarField.H"
#include "dictionaryContent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    baffle_(nullptr),
    dict_()
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const vibrationShellFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>
    (
        ptf,
        p,
        iF,
        mapper
    ),
    baffle_(nullptr),
    dict_(ptf.dict_)
{}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    baffle_(nullptr),
    dict_
    (
        // Copy dictionary, but without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value", "refValue", "refGradient", "valueFraction"
            })
        )
    )
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1;
    }

    if (!baffle_)
    {
        baffle_.reset(baffleType::New(p.boundaryMesh().mesh(), dict_));
    }
}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const vibrationShellFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    baffle_(nullptr),
    dict_(ptf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vibrationShellFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const auto& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar rho("rho", dimDensity, transportProperties);

    baffle_->evolve();

    // rho * acceleration

    refGrad() = Zero;  // safety (for any unmapped values)

    baffle_->vsm().mapToVolumePatch(baffle_->a(), refGrad(), patch().index());

    refGrad() *= rho.value();

    refValue() = Zero;
    valueFraction() = Zero;

    mixedFvPatchField<scalar>::updateCoeffs();
}


void vibrationShellFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    vibrationShellFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
