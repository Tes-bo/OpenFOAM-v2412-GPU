/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "scaledFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::scaledFixedValueFvPatchField<Type>::scaledFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    scalePtr_(),
    refValuePtr_(fvPatchField<Type>::New("refValue", p, iF))
{}


template<class Type>
Foam::scaledFixedValueFvPatchField<Type>::scaledFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    scalePtr_(PatchFunction1<scalar>::New(p.patch(), "scale", dict)),
    refValuePtr_(fvPatchField<Type>::New(p, iF, dict.subDict("refValue")))
{
    if (!isA<fixedValueFvPatchField<Type>>(refValuePtr_()))
    {
        FatalIOErrorInFunction(dict)
            << typeName << " condition can only be applied to fixed value "
            << "conditions"
            << exit(FatalIOError);
    }

    const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));
    this->operator==(s*refValuePtr_());
}


template<class Type>
Foam::scaledFixedValueFvPatchField<Type>::scaledFixedValueFvPatchField
(
    const scaledFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    scalePtr_(ptf.scalePtr_.clone(p.patch())),
    refValuePtr_(fvPatchField<Type>::New(ptf.refValue(), p, iF, mapper))
{}


template<class Type>
Foam::scaledFixedValueFvPatchField<Type>::scaledFixedValueFvPatchField
(
    const scaledFixedValueFvPatchField& spf
)
:
    fixedValueFvPatchField<Type>(spf),
    scalePtr_(spf.scalePtr_.clone(spf.patch().patch())),
    refValuePtr_(spf.refValue().clone())
{}


template<class Type>
Foam::scaledFixedValueFvPatchField<Type>::scaledFixedValueFvPatchField
(
    const scaledFixedValueFvPatchField& spf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(spf, iF),
    scalePtr_(spf.scalePtr_.clone(spf.patch().patch())),
    refValuePtr_(spf.refValue().clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    refValuePtr_->autoMap(m);

    scalePtr_().autoMap(m);

    if (scalePtr_().constant())
    {
        // If mapper is not dependent on time we're ok to evaluate
        this->evaluate();
    }
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const scaledFixedValueFvPatchField& sptf =
        refCast<const scaledFixedValueFvPatchField>(ptf);

    refValuePtr_->rmap(sptf.refValue(), addr);

    scalePtr_().rmap(sptf.scalePtr_(), addr);
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    refValuePtr_->evaluate();

    const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));

    // Note: setting this field value using = operator (not ==)
    Field<Type>::operator=(s*refValuePtr_());

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    scalePtr_->writeData(os);

    os.beginBlock("refValue");
    refValuePtr_->write(os);
    os.endBlock();
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));

    forAll(s, facei)
    {
        const scalar si = s[facei];
        if (mag(si) > ROOTVSMALL)
        {
            refValuePtr_->operator[](facei) = ptf[facei]/si;
        }
    }

    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::operator==(const Field<Type>& tf)
{
    const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));

    forAll(s, facei)
    {
        const scalar si = s[facei];
        if (mag(si) > ROOTVSMALL)
        {
            refValuePtr_->operator[](facei) = tf[facei]/si;
        }
    }

    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::scaledFixedValueFvPatchField<Type>::operator==(const Type& t)
{
    const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));

    forAll(s, facei)
    {
        const scalar si = s[facei];
        if (mag(si) > ROOTVSMALL)
        {
            refValuePtr_->operator[](facei) = t/si;
        }
    }

    Field<Type>::operator=(t);
}


// ************************************************************************* //
