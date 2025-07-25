/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
inline constexpr Foam::Field<Type>::Field() noexcept
:
    List<Type>()
{}


template<class Type>
inline Foam::Field<Type>::Field(const label len)
:
    List<Type>(len)
{}


template<class Type>
inline Foam::Field<Type>::Field(const label len, const Type& val)
:
    List<Type>(len, val)
{}


template<class Type>
inline Foam::Field<Type>::Field(const label len, const Foam::zero)
:
    List<Type>(len, Foam::zero{})
{}


template<class Type>
inline Foam::Field<Type>::Field(const Foam::one, const Type& val)
:
    List<Type>(Foam::one{}, val)
{}


template<class Type>
inline Foam::Field<Type>::Field(const Foam::one, Type&& val)
:
    List<Type>(Foam::one{}, std::move(val))
{}


template<class Type>
inline Foam::Field<Type>::Field(const Foam::one, const Foam::zero)
:
    List<Type>(Foam::one{}, Foam::zero{})
{}


template<class Type>
inline Foam::Field<Type>::Field(const Field<Type>& fld)
:
    List<Type>(fld)
{}


template<class Type>
inline Foam::Field<Type>::Field(const UList<Type>& list)
:
    List<Type>(list)
{}


template<class Type>
template<class Addr>
inline Foam::Field<Type>::Field(const IndirectListBase<Type, Addr>& list)
:
    List<Type>(list)
{}


template<class Type>
inline Foam::Field<Type>::Field(Field<Type>&& fld) noexcept
:
    List<Type>(std::move(static_cast<List<Type>&>(fld)))
{}


template<class Type>
inline Foam::Field<Type>::Field(List<Type>&& list) noexcept
:
    List<Type>(std::move(list))
{}


template<class Type>
template<int SizeMin>
inline Foam::Field<Type>::Field(DynamicList<Type, SizeMin>&& list)
:
    List<Type>()
{
    List<Type>::transfer(list);
}


template<class Type>
inline Foam::Field<Type>::Field(Field<Type>& fld, bool reuse)
:
    List<Type>(fld, reuse)
{}


template<class Type>
inline Foam::Field<Type>::Field(const tmp<Field<Type>>& tfld)
:
    List<Type>(tfld.constCast(), tfld.movable())
{
    tfld.clear();
}


template<class Type>
inline Foam::Field<Type>::Field(Istream& is)
:
    List<Type>(is)
{}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Field<Type>::clone() const
{
    return tmp<Field<Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::Field<Type>::operator=(const UList<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::Field<Type>::operator=(const SubField<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


template<class Type>
template<class Addr>
inline void Foam::Field<Type>::operator=
(
    const IndirectListBase<Type, Addr>& rhs
)
{
    List<Type>::operator=(rhs);
}


template<class Type>
inline void Foam::Field<Type>::operator=(Field<Type>&& rhs)
{
    List<Type>::transfer(rhs);
}


template<class Type>
inline void Foam::Field<Type>::operator=(List<Type>&& rhs)
{
    List<Type>::transfer(rhs);
}


template<class Type>
template<int SizeMin>
inline void Foam::Field<Type>::operator=(DynamicList<Type, SizeMin>&& rhs)
{
    List<Type>::transfer(rhs);
}


template<class Type>
inline void Foam::Field<Type>::operator=(const Type& val)
{
    List<Type>::operator=(val);
}


template<class Type>
inline void Foam::Field<Type>::operator=(const Foam::zero)
{
    List<Type>::operator=(Foam::zero{});
}


// ************************************************************************* //
