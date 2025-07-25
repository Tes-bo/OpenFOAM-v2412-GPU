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

#include "mappedMixedFieldFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedMixedFieldFvPatchField<Type>::mappedMixedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedPatchFieldBase<Type>(*this, *this),
    weightFieldName_()
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::mappedMixedFieldFvPatchField<Type>::mappedMixedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedPatchFieldBase<Type>(*this, *this, dict),
    weightFieldName_(dict.getOrDefault<word>("weightField", word::null))
{}


template<class Type>
Foam::mappedMixedFieldFvPatchField<Type>::mappedMixedFieldFvPatchField
(
    const mappedMixedFieldFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


template<class Type>
Foam::mappedMixedFieldFvPatchField<Type>::mappedMixedFieldFvPatchField
(
    const mappedMixedFieldFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


template<class Type>
Foam::mappedMixedFieldFvPatchField<Type>::mappedMixedFieldFvPatchField
(
    const mappedMixedFieldFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedMixedFieldFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedMixedFieldFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedMixedFieldFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const tmp<Field<Type>> nbrIntFld(this->mappedInternalField());

    //- Unweighted
    //const tmp<scalarField> nbrKDelta(this->mappedWeightField());

    //- Weighted
    tmp<scalarField> myKDelta;
    tmp<scalarField> nbrKDelta;
    this->mappedWeightField(weightFieldName_, myKDelta, nbrKDelta);

    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld;
    this->refGrad() = Zero;
    this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

    mixedFvPatchField<Type>::updateCoeffs();

    if (debug)
    {
        Info<< this->patch().boundaryMesh().mesh().name() << ':'
            << this->patch().name() << ':'
            << this->internalField().name() << " <- "
            << this->mapper_.sampleRegion() << ':'
            << this->mapper_.samplePatch() << ':'
            << this->fieldName_ << " :"
            << " value "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


template<class Type>
void Foam::mappedMixedFieldFvPatchField<Type>::write(Ostream& os) const
{
    mappedPatchBase::write(os);
    mappedPatchFieldBase<Type>::write(os);
    os.writeEntryIfDifferent<word>("weightField", word::null, weightFieldName_);
    mixedFvPatchField<Type>::write(os);
}


// ************************************************************************* //
