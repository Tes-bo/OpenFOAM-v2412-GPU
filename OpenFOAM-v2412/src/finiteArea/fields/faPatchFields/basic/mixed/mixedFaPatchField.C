/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "mixedFaPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::mixedFaPatchField<Type>::readMixedEntries
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = faPatchFieldBase::patch();


    // If there is a 'refValue', also require all others
    const auto* hasValue = dict.findEntry("refValue", keyType::LITERAL);

    if (!hasValue && IOobjectOption::isReadOptional(readOpt))
    {
        return false;
    }

    const auto* hasGrad = dict.findEntry("refGradient", keyType::LITERAL);
    const auto* hasFrac = dict.findEntry("valueFraction", keyType::LITERAL);

    // Combined error message on failure
    if (!hasValue || !hasGrad || !hasFrac)
    {
        FatalIOErrorInFunction(dict)
            << "Required entries:";

        if (!hasValue) FatalIOError << " 'refValue'";
        if (!hasGrad)  FatalIOError << " 'refGradient'";
        if (!hasFrac)  FatalIOError << " 'valueFraction'";

        FatalIOError
            << " : missing for patch " << p.name()
            << " : in dictionary " << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    // Everything verified - can assign
    refValue_.assign(*hasValue, p.size());
    refGrad_.assign(*hasGrad, p.size());
    valueFraction_.assign(*hasFrac, p.size());

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{}


template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Foam::zero
)
:
    faPatchField<Type>(p, iF),
    refValue_(p.size(), Zero),
    refGrad_(p.size(), Zero),
    valueFraction_(p.size(), Zero)
{}


template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireMixed
)
:
    // The "value" entry is not required
    faPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{
    if (!readMixedEntries(dict, requireMixed))
    {
        // Not read (eg, optional and missing): no evaluate possible/need
        return;
    }

    // Could also check/clamp fraction to 0-1 range
    evaluate();
}


template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


template<class Type>
Foam::mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mixedFaPatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
}


template<class Type>
void Foam::mixedFaPatchField<Type>::rmap
(
    const faPatchField<Type>& ptf,
    const labelList& addr
)
{
    faPatchField<Type>::rmap(ptf, addr);

    const auto& mptf = refCast<const mixedFaPatchField<Type>>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
}


template<class Type>
void Foam::mixedFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        lerp
        (
            this->patchInternalField() + refGrad_/this->patch().deltaCoeffs(),
            refValue_,
            valueFraction_
        )
    );

    faPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mixedFaPatchField<Type>::snGrad() const
{
    return lerp
    (
        refGrad_,
        (refValue_ - this->patchInternalField())*this->patch().deltaCoeffs(),
        valueFraction_
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mixedFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mixedFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return lerp
    (
        refGrad_/this->patch().deltaCoeffs(),
        refValue_,
        valueFraction_
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}


template<class Type>
void Foam::mixedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
