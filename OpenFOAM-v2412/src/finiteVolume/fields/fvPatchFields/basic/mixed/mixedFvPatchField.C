/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::mixedFvPatchField<Type>::readMixedEntries
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = fvPatchFieldBase::patch();


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
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size()),
    source_(p.size(), Zero)
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Foam::zero
)
:
    fvPatchField<Type>(p, iF),
    refValue_(p.size(), Zero),
    refGrad_(p.size(), Zero),
    valueFraction_(p.size(), Zero),
    source_(p.size(), Zero)
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireMixed
)
:
    // The "value" entry is not required
    fvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size()),
    source_(p.size(), Zero)
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
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper),
    source_(ptf.source_, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
    source_.autoMap(m);
}


template<class Type>
void Foam::mixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const mixedFvPatchField<Type>& mptf =
        refCast<const mixedFvPatchField<Type>>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
    source_.rmap(mptf.source_, addr);
}


template<class Type>
void Foam::mixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
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

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::snGrad() const
{
    return lerp
    (
        refGrad_,
        (refValue_ - this->patchInternalField())*this->patch().deltaCoeffs(),
        valueFraction_
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueBoundaryCoeffs
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
Foam::mixedFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}


template<class Type>
void Foam::mixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    source_.writeEntry("source", os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
