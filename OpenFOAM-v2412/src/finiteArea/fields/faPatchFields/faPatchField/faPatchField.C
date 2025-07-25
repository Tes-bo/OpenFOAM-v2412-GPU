/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "dictionary.H"
#include "faPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::faPatchField<Type>::readValueEntry
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = faPatchFieldBase::patch();


    const auto* eptr = dict.findEntry("value", keyType::LITERAL);

    if (eptr)
    {
        Field<Type>::assign(*eptr, p.size());
        return true;
    }

    if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Required entry 'value' : missing for patch " << p.name()
            << " in dictionary " << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class Type>
void Foam::faPatchField<Type>::extrapolateInternal()
{
    faPatchFieldBase::patch().patchInternalField(internalField_, *this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchFieldBase(p),
    Field<Type>(p.size()),
    internalField_(iF)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Type& value
)
:
    faPatchFieldBase(p),
    Field<Type>(p.size(), value),
    internalField_(iF)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& pfld
)
:
    faPatchFieldBase(p),
    Field<Type>(pfld),
    internalField_(iF)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    Field<Type>&& pfld
)
:
    faPatchFieldBase(p),
    Field<Type>(std::move(pfld)),
    internalField_(iF)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    faPatchFieldBase(p, dict),
    Field<Type>(p.size()),
    internalField_(iF)
{
    if (!readValueEntry(dict, requireValue))
    {
        // Not read (eg, optional and missing): define zero
        Field<Type>::operator=(Zero);
    }
}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchFieldBase(ptf, p),
    Field<Type>(ptf, mapper),
    internalField_(iF)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf
)
:
    faPatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(ptf.internalField_)
{}


template<class Type>
Foam::faPatchField<Type>::faPatchField
(
    const faPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::faPatchField<Type>::check(const faPatchField<Type>& rhs) const
{
    faPatchFieldBase::checkPatch(rhs);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faPatchField<Type>::snGrad() const
{
    return (*this - patchInternalField())*patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::faPatchField<Type>::patchInternalField() const
{
    return patch().patchInternalField(internalField_);
}


template<class Type>
void Foam::faPatchField<Type>::patchInternalField(Field<Type>& pfld) const
{
    patch().patchInternalField(internalField_, pfld);
}


template<class Type>
void Foam::faPatchField<Type>::autoMap(const faPatchFieldMapper& m)
{
    Field<Type>::autoMap(m);
}


template<class Type>
void Foam::faPatchField<Type>::rmap
(
    const faPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::faPatchField<Type>::updateCoeffs()
{
    faPatchFieldBase::setUpdated(true);
}


template<class Type>
void Foam::faPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    faPatchFieldBase::setUpdated(false);
}


template<class Type>
void Foam::faPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (!patchType().empty())
    {
        os.writeEntry("patchType", patchType());
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::faPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::faPatchField<Type>::operator=
(
    const faPatchField<Type>& ptf
)
{
    faPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator+=
(
    const faPatchField<Type>& ptf
)
{
    faPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator-=
(
    const faPatchField<Type>& ptf
)
{
    faPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator*=
(
    const faPatchField<scalar>& ptf
)
{
    faPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator/=
(
    const faPatchField<scalar>& ptf
)
{
    faPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::faPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::faPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::faPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::faPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::faPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::faPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::faPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::faPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::faPatchField<Type>::operator==
(
    const faPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::faPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::faPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const faPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
