/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "prghTotalHydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prghTotalHydrostaticPressureFvPatchScalarField::
prghTotalHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    ph_rghName_("ph_rgh")
{}


Foam::prghTotalHydrostaticPressureFvPatchScalarField::
prghTotalHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.getOrDefault<word>("U", "U")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    ph_rghName_(dict.getOrDefault<word>("ph_rgh", "ph_rgh"))
{}


Foam::prghTotalHydrostaticPressureFvPatchScalarField::
prghTotalHydrostaticPressureFvPatchScalarField
(
    const prghTotalHydrostaticPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    ph_rghName_(ptf.ph_rghName_)
{}


Foam::prghTotalHydrostaticPressureFvPatchScalarField::
prghTotalHydrostaticPressureFvPatchScalarField
(
    const prghTotalHydrostaticPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    ph_rghName_(ptf.ph_rghName_)
{}


Foam::prghTotalHydrostaticPressureFvPatchScalarField::
prghTotalHydrostaticPressureFvPatchScalarField
(
    const prghTotalHydrostaticPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    ph_rghName_(ptf.ph_rghName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prghTotalHydrostaticPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& rhop =
        patch().lookupPatchField<volScalarField>(rhoName_);

    const scalarField& ph_rghp =
        patch().lookupPatchField<volScalarField>(ph_rghName_);

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField>(phiName_);

    const vectorField& Up =
        patch().lookupPatchField<volVectorField>(UName_);

    operator==
    (
        ph_rghp
      - 0.5*rhop*(neg(phip))*magSqr(Up)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::prghTotalHydrostaticPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("ph_rgh", "ph_rgh", ph_rghName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        prghTotalHydrostaticPressureFvPatchScalarField
    );
}

// ************************************************************************* //
