/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "uniformNormalFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformNormalFixedValueFvPatchVectorField::
uniformNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValueFunc_(nullptr),
    ramp_(nullptr)
{}


Foam::uniformNormalFixedValueFvPatchVectorField::
uniformNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    refValueFunc_(PatchFunction1<scalar>::New(p.patch(), "uniformValue", dict)),
    ramp_(Function1<scalar>::NewIfPresent("ramp", dict, word::null, &db()))
{
    if (!this->readValueEntry(dict))
    {
        this->evaluate();
    }
}


Foam::uniformNormalFixedValueFvPatchVectorField::
uniformNormalFixedValueFvPatchVectorField
(
    const uniformNormalFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),   // Don't map
    refValueFunc_(ptf.refValueFunc_.clone(p.patch())),
    ramp_(ptf.ramp_.clone())
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


Foam::uniformNormalFixedValueFvPatchVectorField::
uniformNormalFixedValueFvPatchVectorField
(
    const uniformNormalFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    refValueFunc_(ptf.refValueFunc_.clone(this->patch().patch())),
    ramp_(ptf.ramp_.clone())
{}


Foam::uniformNormalFixedValueFvPatchVectorField::
uniformNormalFixedValueFvPatchVectorField
(
    const uniformNormalFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    refValueFunc_(ptf.refValueFunc_.clone(this->patch().patch())),
    ramp_(ptf.ramp_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformNormalFixedValueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
    refValueFunc_().autoMap(mapper);

    if (refValueFunc_().constant())
    {
        // If mapper is not dependent on time we're ok to evaluate
        this->evaluate();
    }
}


void Foam::uniformNormalFixedValueFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const uniformNormalFixedValueFvPatchVectorField& tiptf =
        refCast<const uniformNormalFixedValueFvPatchVectorField>(ptf);

    refValueFunc_().rmap(tiptf.refValueFunc_(), addr);
}


void Foam::uniformNormalFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    tmp<vectorField> tvalues(refValueFunc_->value(t)*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(t);
    }

    fvPatchVectorField::operator=(tvalues);
    fvPatchVectorField::updateCoeffs();
}


void Foam::uniformNormalFixedValueFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    refValueFunc_->writeData(os);
    if (ramp_)
    {
        ramp_->writeData(os);
    }
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        uniformNormalFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
