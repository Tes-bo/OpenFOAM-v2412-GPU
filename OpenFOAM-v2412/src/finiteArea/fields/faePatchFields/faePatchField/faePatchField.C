/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "faePatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::faePatchField<Type>::readValueEntry
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = faePatchFieldBase::patch();


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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchFieldBase(p),
    Field<Type>(p.size()),
    internalField_(iF)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const Type& value
)
:
    faePatchFieldBase(p),
    Field<Type>(p.size(), value),
    internalField_(iF)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const Field<Type>& pfld
)
:
    faePatchFieldBase(p),
    Field<Type>(pfld),
    internalField_(iF)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    Field<Type>&& pfld
)
:
    faePatchFieldBase(p),
    Field<Type>(std::move(pfld)),
    internalField_(iF)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    faePatchFieldBase(p, dict),
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
Foam::faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchFieldBase(ptf, p),
    Field<Type>(ptf, mapper),
    internalField_(iF)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf
)
:
    faePatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(ptf.internalField_)
{}


template<class Type>
Foam::faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchFieldBase(ptf),
    Field<Type>(ptf),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::faePatchField<Type>::check(const faePatchField<Type>& rhs) const
{
    faePatchFieldBase::checkPatch(rhs);
}


template<class Type>
void Foam::faePatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


template<class Type>
void Foam::faePatchField<Type>::rmap
(
    const faePatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::faePatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    // if (!patchType().empty())
    // {
    //     os.writeEntry("patchType", patchType());
    // }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::faePatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::faePatchField<Type>::operator=
(
    const faePatchField<Type>& ptf
)
{
    faePatchFieldBase::checkPatch(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator+=
(
    const faePatchField<Type>& ptf
)
{
    faePatchFieldBase::checkPatch(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator-=
(
    const faePatchField<Type>& ptf
)
{
    faePatchFieldBase::checkPatch(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator*=
(
    const faePatchField<scalar>& ptf
)
{
    faePatchFieldBase::checkPatch(ptf);
    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator/=
(
    const faePatchField<scalar>& ptf
)
{
    faePatchFieldBase::checkPatch(ptf);
    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::faePatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::faePatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::faePatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::faePatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::faePatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::faePatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::faePatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::faePatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::faePatchField<Type>::operator==
(
    const faePatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::faePatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::faePatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const faePatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
