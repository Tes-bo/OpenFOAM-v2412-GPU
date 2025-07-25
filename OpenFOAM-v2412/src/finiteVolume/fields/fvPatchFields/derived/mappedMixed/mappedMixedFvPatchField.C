/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "mappedMixedFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"
#include "mappedFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>
    (
        mappedFixedValueFvPatchField<Type>::mapper(p, iF),
        *this
    ),
    weightFieldName_()
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    // Bypass dictionary constructor (all reading handled later)
    // but cannot use NO_READ since will still trigger an evaluate()
    mixedFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>
    (
        mappedFixedValueFvPatchField<Type>::mapper(p, iF),
        *this,
        dict
    ),
    weightFieldName_(dict.getOrDefault<word>("weightField", word::null))
{
    fvPatchFieldBase::readDict(dict);  // Consistent with a dict constructor

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        this->refValue() = *this;
        this->refGrad() = Zero;
        this->valueFraction() = 1.0;
    }

// This blocks (crashes) with more than two worlds!
//
///    // Store patch value as initial guess when running in database mode
///    mappedPatchFieldBase<Type>::initRetrieveField
///    (
///        this->internalField().name(),
///        *this
///    );
///    mappedPatchFieldBase<Type>::initRetrieveField
///    (
///        this->internalField().name() + "_weights",
///        this->patch().deltaCoeffs()
///    );
}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchFieldBase<Type>
    (
        mappedFixedValueFvPatchField<Type>::mapper(p, iF),
        *this,
        ptf
    ),
    weightFieldName_(ptf.weightFieldName_)
{}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    mappedPatchFieldBase<Type>(ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    mappedPatchFieldBase<Type>
    (
        mappedFixedValueFvPatchField<Type>::mapper(ptf.patch(), iF),
        *this,
        ptf
    ),
    weightFieldName_(ptf.weightFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedMixedFvPatchField<Type>::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::mappedMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::mappedMixedFvPatchField<Type>::updateCoeffs()
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
void Foam::mappedMixedFvPatchField<Type>::write(Ostream& os) const
{
    mappedPatchFieldBase<Type>::write(os);
    os.writeEntryIfDifferent<word>("weightField", word::null, weightFieldName_);
    mixedFvPatchField<Type>::write(os);
}


// ************************************************************************* //
