/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "uniformInletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    uniformInletValue_
    (
        Function1<Type>::New("uniformInletValue", dict, &this->db())
    )
{
    fvPatchFieldBase::readDict(dict);

    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());

    if (!this->readValueEntry(dict))
    {
        fvPatchField<Type>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(p, iF),  // Don't map
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_.clone())
{
    this->patchType() = ptf.patchType();

    // Evaluate refValue since not mapped
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    // Initialize the patch value to the refValue
    fvPatchField<Type>::operator=(this->refValue());

    this->map(ptf, mapper);
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_.clone())
{}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update the uniform value as a function of time
    const scalar t = this->db().time().timeOutputValue();
    this->refValue() = uniformInletValue_->value(t);

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField>
        (
            phiName_
        );

    this->valueFraction() = neg(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    this->uniformInletValue_->writeData(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().timeOutputValue());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
