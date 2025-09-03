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

#include "fixedNormalSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    fixedValue_(p.size(), Foam::zero{}),
    writeValue_(false)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    fixedValue_(ptf.fixedValue_, mapper),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    fixedValue_("fixedValue", dict, p.size()),
    writeValue_(dict.getOrDefault("writeValue", false))
{
    fvPatchFieldBase::readDict(dict);

    evaluate();
}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf
)
:
    parent_bctype(ptf),
    fixedValue_(ptf.fixedValue_),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    fixedValue_(ptf.fixedValue_),
    writeValue_(ptf.writeValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    parent_bctype::autoMap(m);
    fixedValue_.autoMap(m);
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    parent_bctype::rmap(ptf, addr);

    const auto& dmptf =
        refCast<const fixedNormalSlipFvPatchField<Type>>(ptf);

    fixedValue_.rmap(dmptf.fixedValue_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    const Field<Type> pif(this->patchInternalField());

    return
    (
        (nHat*(nHat & fixedValue_) + transform(I - sqr(nHat), pif)) - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().nf());

    Field<Type>::operator=
    (
        nHat*(nHat & fixedValue_)
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    this->parent_bctype::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().nf()));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    if (writeValue_)
    {
        os.writeEntry("writeValue", "true");
    }
    fixedValue_.writeEntry("fixedValue", os);
    if (writeValue_)
    {
        fvPatchField<Type>::writeValueEntry(os);
    }
}


// ************************************************************************* //
