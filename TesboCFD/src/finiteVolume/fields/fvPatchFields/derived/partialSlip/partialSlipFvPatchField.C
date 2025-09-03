/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "partialSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Foam::zero{}),
    valueFraction_(p.size(), 1.0),
    writeValue_(false)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    valueFraction_(ptf.valueFraction_, mapper),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Foam::zero{}),
    valueFraction_("valueFraction", dict, p.size()),
    writeValue_(dict.getOrDefault("writeValue", false))
{
    fvPatchFieldBase::readDict(dict);

    // Backwards compatibility - leave refValue as zero unless specified
    refValue_.assign("refValue", dict, p.size(), IOobjectOption::LAZY_READ);

    evaluate();
}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf
)
:
    partialSlipFvPatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::partialSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    parent_bctype::autoMap(m);
    refValue_.autoMap(m);
    valueFraction_.autoMap(m);
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    parent_bctype::rmap(ptf, addr);

    const auto& dmptf =
        refCast<const partialSlipFvPatchField<Type>>(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    valueFraction_.rmap(dmptf.valueFraction_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::partialSlipFvPatchField<Type>::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    const Field<Type> pif(this->patchInternalField());

    return
    (
        lerp
        (
            transform(I - sqr(nHat), pif),
            refValue_,
            valueFraction_
        ) - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    tmp<vectorField> nHat = this->patch().nf();

    Field<Type>::operator=
    (
        lerp
        (
            transform(I - sqr(nHat), this->patchInternalField()),
            refValue_,
            valueFraction_
        )
    );

    parent_bctype::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::partialSlipFvPatchField<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().nf()));

    return
        valueFraction_*pTraits<Type>::one
      + (1.0 - valueFraction_)
       *transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    if (writeValue_)
    {
        os.writeEntry("writeValue", "true");
    }
    refValue_.writeEntry("refValue", os);
    valueFraction_.writeEntry("valueFraction", os);
    if (writeValue_)
    {
        fvPatchField<Type>::writeValueEntry(os);
    }
}


// ************************************************************************* //
