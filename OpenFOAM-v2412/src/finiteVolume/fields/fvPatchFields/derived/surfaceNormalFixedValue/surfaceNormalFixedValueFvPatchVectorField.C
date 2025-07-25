/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceNormalFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(p.size()),
    ramp_(nullptr)
{}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    refValue_("refValue", dict, p.size()),
    ramp_(Function1<scalar>::NewIfPresent("ramp", dict, word::null, &db()))
{
    tmp<vectorField> tvalues(refValue_*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }

    fvPatchVectorField::operator=(tvalues);
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(ptf.refValue_, mapper, pTraits<scalar>::zero),
    ramp_(ptf.ramp_.clone())
{
    // Note: refValue_ will have default value of 0 for unmapped faces. This
    // can temporarily happen during e.g. redistributePar. We should not
    // access ptf.patch() instead since redistributePar has destroyed this
    // at the time of mapping.

    tmp<vectorField> tvalues(refValue_*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }

    fvPatchVectorField::operator=(tvalues);
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    refValue_(ptf.refValue_),
    ramp_(ptf.ramp_.clone())
{}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    refValue_(ptf.refValue_),
    ramp_(ptf.ramp_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceNormalFixedValueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
    refValue_.autoMap(mapper);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const surfaceNormalFixedValueFvPatchVectorField& tiptf =
        refCast<const surfaceNormalFixedValueFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> tvalues(refValue_*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }

    fvPatchVectorField::operator=(tvalues);
    fvPatchVectorField::updateCoeffs();
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    refValue_.writeEntry("refValue", os);
    if (ramp_)
    {
        ramp_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        surfaceNormalFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
