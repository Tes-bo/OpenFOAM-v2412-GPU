/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "freestreamPressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U")
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.getOrDefault<word>("U", "U"))
{
    // freestreamValue() and refValue() are identical
    freestreamValue().assign("freestreamValue", dict, p.size());
    refGrad() = Zero;
    valueFraction() = 0;

    if (!this->readValueEntry(dict))
    {
        fvPatchScalarField::operator=(freestreamValue());
    }
}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& wbppsf
)
:
    mixedFvPatchScalarField(wbppsf),
    UName_(wbppsf.UName_)
{}


Foam::freestreamPressureFvPatchScalarField::
freestreamPressureFvPatchScalarField
(
    const freestreamPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wbppsf, iF),
    UName_(wbppsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const Field<vector>& Up =
        patch().template lookupPatchField<volVectorField>(UName_);

    valueFraction() = 0.5 + 0.5*(Up & patch().nf())/mag(Up);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::freestreamPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    freestreamValue().writeEntry("freestreamValue", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freestreamPressureFvPatchScalarField
    );
}

// ************************************************************************* //
