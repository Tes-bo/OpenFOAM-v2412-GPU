/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef FULLDEBUG

template<template<class> class Field, class Type1, class Type2>
void checkFields
(
    const FieldField<Field, Type1>& f1,
    const FieldField<Field, Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " FieldField<" << pTraits<Type1>::typeName
            << "> f1(" << f1.size() << ')'
            << " and FieldField<" << pTraits<Type2>::typeName
            << "> f2(" << f2.size() << ')'
            << endl << " for operation " << op
            << abort(FatalError);
    }
}

template<template<class> class Field, class Type1, class Type2, class Type3>
void checkFields
(
    const FieldField<Field, Type1>& f1,
    const FieldField<Field, Type2>& f2,
    const FieldField<Field, Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorInFunction
            << " FieldField<" << pTraits<Type1>::typeName
            << "> f1(" << f1.size() << ')'
            << ", FieldField<" <<pTraits<Type2>::typeName
            << "> f2(" << f2.size() << ')'
            << " and FieldField<"<<pTraits<Type3>::typeName
            << "> f3("<<f3.size() << ')'
            << endl << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<template<class> class Field, class Type1, class Type2>
void checkFields
(
    const FieldField<Field, Type1>&,
    const FieldField<Field, Type2>&,
    const char* op
)
{}

template<template<class> class Field, class Type1, class Type2, class Type3>
void checkFields
(
    const FieldField<Field, Type1>&,
    const FieldField<Field, Type2>&,
    const FieldField<Field, Type3>&,
    const char* op
)
{}

#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
constexpr FieldField<Field, Type>::FieldField() noexcept
:
    PtrList<Field<Type>>()
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(const label size)
:
    PtrList<Field<Type>>(size)
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField
(
    const word& type,
    const FieldField<Field, Type>& ff
)
:
    PtrList<Field<Type>>(ff.size())
{
    forAll(*this, i)
    {
        set(i, Field<Type>::New(type, ff[i]));
    }
}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(const FieldField<Field, Type>& ff)
:
    PtrList<Field<Type>>(ff)
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(FieldField<Field, Type>&& ff)
:
    PtrList<Field<Type>>(std::move(ff))
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(FieldField<Field, Type>& ff, bool reuse)
:
    PtrList<Field<Type>>(ff, reuse)
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(const PtrList<Field<Type>>& list)
:
    PtrList<Field<Type>>(list)
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(PtrList<Field<Type>>&& list)
:
    PtrList<Field<Type>>(std::move(list))
{}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(const tmp<FieldField<Field, Type>>& tf)
:
    PtrList<Field<Type>>(tf.constCast(), tf.movable())
{
    tf.clear();
}


template<template<class> class Field, class Type>
FieldField<Field, Type>::FieldField(Istream& is)
:
    PtrList<Field<Type>>(is)
{}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> FieldField<Field, Type>::clone() const
{
    return tmp<FieldField<Field, Type>>::New(*this);
}


template<template<class> class Field, class Type>
template<class Type2>
tmp<FieldField<Field, Type>> FieldField<Field, Type>::NewCalculatedType
(
    const FieldField<Field, Type2>& ff
)
{
    const label len = ff.size();

    auto tresult = tmp<FieldField<Field, Type>>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result.set(i, Field<Type>::NewCalculatedType(ff[i]).ptr());
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void FieldField<Field, Type>::negate()
{
    forAll(*this, i)
    {
        this->operator[](i).negate();
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::normalise()
{
    forAll(*this, i)
    {
        this->operator[](i).normalise();
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType>>
FieldField<Field, Type>::component
(
    const direction d
) const
{
    auto tres =
        FieldField
        <
            Field, typename FieldField<Field, Type>::cmptType
        >::NewCalculatedType(*this);

    ::Foam::component(tres.ref(), *this, d);

    return tres;
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::replace
(
    const direction d,
    const FieldField<Field, cmptType>& sf
)
{
    forAll(*this, i)
    {
        this->operator[](i).replace(d, sf[i]);
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::replace
(
    const direction d,
    const cmptType& s
)
{
    forAll(*this, i)
    {
        this->operator[](i).replace(d, s);
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::clamp_min
(
    const Type& lower
)
{
    for (auto& ff : *this)
    {
        ff.clamp_min(lower);
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::clamp_max
(
    const Type& upper
)
{
    for (auto& ff : *this)
    {
        ff.clamp_max(upper);
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::clamp_range
(
    const Type& lower,
    const Type& upper
)
{
    // Note: no checks for bad/invalid clamping ranges

    for (auto& ff : *this)
    {
        ff.clamp_range(lower, upper);
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::clamp_range
(
    const MinMax<Type>& range
)
{
    // Note: no checks for bad/invalid clamping ranges

    for (auto& ff : *this)
    {
        ff.clamp_range(range.min(), range.max());
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> FieldField<Field, Type>::T() const
{
    auto tres
    (
        FieldField<Field, Type>::NewCalculatedType(*this)
    );

    ::Foam::T(tres.ref(), *this);
    return tres;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
const Type& FieldField<Field, Type>::operator[](const labelPair& index) const
{
    return this->operator[](index.first())[index.second()];
}


template<template<class> class Field, class Type>
Type& FieldField<Field, Type>::operator[](const labelPair& index)
{
    return this->operator[](index.first())[index.second()];
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::operator=(const FieldField<Field, Type>& ff)
{
    if (this == &ff)
    {
        return;  // Self-assignment is a no-op
    }

    // No size checking done

    forAll(*this, i)
    {
        this->operator[](i) = ff[i];
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::operator=(FieldField<Field, Type>&& ff)
{
    if (this == &ff)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<Field<Type>>::transfer(ff);
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::operator=(const tmp<FieldField>& tf)
{
    // The cref() method also checks that tmp is not nullptr.
    if (this == &(tf.cref()))
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<Field<Type>>::clear();

    // Release the tmp pointer, or clone const reference for a new pointer.
    // Error potential when tmp is non-unique.

    auto* tptr = tf.ptr();
    PtrList<Field<Type>>::transfer(*tptr);
    delete tptr;
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::operator=(const Type& val)
{
    forAll(*this, i)
    {
        this->operator[](i) = val;
    }
}


template<template<class> class Field, class Type>
void FieldField<Field, Type>::operator=(const Foam::zero)
{
    forAll(*this, i)
    {
        this->operator[](i) = Foam::zero{};
    }
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<template<class> class Field, class Type>                              \
void FieldField<Field, Type>::operator op(const FieldField<Field, TYPE>& f)    \
{                                                                              \
    forAll(*this, i)                                                           \
    {                                                                          \
        this->operator[](i) op f[i];                                           \
    }                                                                          \
}                                                                              \
                                                                               \
template<template<class> class Field, class Type>                              \
void FieldField<Field, Type>::operator op                                      \
(                                                                              \
    const tmp<FieldField<Field, TYPE>>& tf                                     \
)                                                                              \
{                                                                              \
    operator op(tf());                                                         \
    tf.clear();                                                                \
}                                                                              \
                                                                               \
template<template<class> class Field, class Type>                              \
void FieldField<Field, Type>::operator op(const TYPE& t)                       \
{                                                                              \
    forAll(*this, i)                                                           \
    {                                                                          \
        this->operator[](i) op t;                                              \
    }                                                                          \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
Ostream& operator<<(Ostream& os, const FieldField<Field, Type>& f)
{
    os << static_cast<const PtrList<Field<Type>>&>(f);
    return os;
}


template<template<class> class Field, class Type>
Ostream& operator<<(Ostream& os, const tmp<FieldField<Field, Type>>& tf)
{
    os << tf();
    tf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFieldFunctions.C"

// ************************************************************************* //
