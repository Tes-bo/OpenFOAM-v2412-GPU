/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "uniformFixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    refGradFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedGradientFvPatchField<Type>(p, iF, fld),
    refGradFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),  // Bypass dictionary constructor
    refGradFunc_
    (
        PatchFunction1<Type>::New(p.patch(), "uniformGradient", dict)
    )
{
    fvPatchFieldBase::readDict(dict);

    if (!this->readValueEntry(dict))
    {
        // Ensure field has initialised values
        this->extrapolateInternal();

        // Evaluate to assign a value
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    refGradFunc_(ptf.refGradFunc_.clone())
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    refGradFunc_(ptf.refGradFunc_.clone())
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    refGradFunc_(ptf.refGradFunc_.clone())
{
    // Evaluate the profile if defined
    if (refGradFunc_)
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    // Extra safety on the evaluation
    if (refGradFunc_)
    {
        this->gradient() = refGradFunc_->value(t);
    }
    else
    {
        this->gradient() = Zero;
    }

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    if (refGradFunc_)
    {
        refGradFunc_->writeData(os);
    }
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
