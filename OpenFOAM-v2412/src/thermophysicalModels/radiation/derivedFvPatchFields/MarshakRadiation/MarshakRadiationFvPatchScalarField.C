/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016,2020 OpenCFD Ltd.
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

#include "MarshakRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationModel.H"
#include "physicoChemicalConstants.H"
#include "boundaryRadiationProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::MarshakRadiationFvPatchScalarField::
MarshakRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::radiation::MarshakRadiationFvPatchScalarField::
MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_)
{}


Foam::radiation::MarshakRadiationFvPatchScalarField::
MarshakRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.getOrDefault<word>("T", "T"))
{
    if (!this->readValueEntry(dict))
    {
        fvPatchScalarField::operator=(Zero);
    }

    refValue() = *this;
    refGrad() = Zero;  // zero gradient
    valueFraction() = 1.0;
    fvPatchScalarField::operator=(refValue());
}


Foam::radiation::MarshakRadiationFvPatchScalarField::
MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_)
{}


Foam::radiation::MarshakRadiationFvPatchScalarField::
MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::MarshakRadiationFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Temperature field
    const auto& Tp = patch().lookupPatchField<volScalarField>(TName_);

    // Re-calc reference value
    refValue() = 4.0*constant::physicoChemical::sigma.value()*pow4(Tp);

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const auto& gamma = patch().lookupPatchField<volScalarField>("gammaRad");

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());

    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index(), 0, nullptr, &Tp)
    );

    const scalarField& emissivity = temissivity();

    const scalarField Ep(emissivity/(2.0*(scalar(2) - emissivity)));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    UPstream::msgType(oldTag);  // Restore tag

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::MarshakRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFvPatchScalarField
    );
}
}

// ************************************************************************* //
