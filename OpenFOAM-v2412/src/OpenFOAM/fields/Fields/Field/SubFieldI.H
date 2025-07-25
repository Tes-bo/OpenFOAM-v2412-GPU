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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubField<Type>& fld
)
:
    SubList<Type>(fld)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubList<Type>& list
)
:
    SubList<Type>(list)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list
) noexcept
:
    SubList<Type>(list)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label len
)
:
    SubList<Type>(list, len)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label len,
    const label start
)
:
    SubList<Type>(list, len, start)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const labelRange& range
)
:
    SubList<Type>(list, range)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const labelRange& range,
    const UList<Type>& list
)
:
    SubList<Type>(range, list)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::Field<typename Foam::SubField<Type>::cmptType>>
Foam::SubField<Type>::component
(
    const direction d
) const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).component(d);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::SubField<Type>::T() const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).T();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline Foam::SubField<Type>::operator const Foam::Field<Type>&() const
{
    return *reinterpret_cast<const Field<Type>*>(this);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const SubField<Type>& rhs)
{
    SubList<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const Field<Type>& rhs)
{
    SubList<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const Type& val)
{
    SubList<Type>::operator=(val);
}


template<class Type>
inline void Foam::SubField<Type>::operator=(const Foam::zero)
{
    SubList<Type>::operator=(Foam::zero{});
}


template<class Type>
template<class Form, Foam::direction Ncmpts>
inline void Foam::SubField<Type>::operator=
(
    const VectorSpace<Form, Type, Ncmpts>& rhs
)
{
    forAll(rhs, i)
    {
        this->operator[](i) = rhs[i];
    }
}


template<class Type>
inline void Foam::SubField<Type>::operator+=(const Type& val)
{
    for (Type& lhs : *this)
    {
        lhs += val;
    }
}


template<class Type>
inline void Foam::SubField<Type>::operator-=(const Type& val)
{
    for (Type& lhs : *this)
    {
        lhs -= val;
    }
}


template<class Type>
inline void Foam::SubField<Type>::operator*=(const scalar& s)
{
    for (Type& lhs : *this)
    {
        lhs *= s;
    }
}


template<class Type>
inline void Foam::SubField<Type>::operator/=(const scalar& s)
{
    for (Type& lhs : *this)
    {
        lhs /= s;
    }
}


// ************************************************************************* //
