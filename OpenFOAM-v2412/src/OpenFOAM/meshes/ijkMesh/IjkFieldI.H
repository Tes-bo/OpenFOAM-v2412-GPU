/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
inline Foam::IjkField<Type>::IjkField()
:
    Field<Type>(),
    ijk_()
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField(const IjkField<Type>& field)
:
    Field<Type>(field),
    ijk_(field.ijk())
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField(IjkField<Type>&& field)
:
    Field<Type>(std::move(field)),
    ijk_(field.ijk())
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField(const labelVector& ijk)
:
    Field<Type>(cmptProduct(ijk)),
    ijk_(ijk)
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField
(
    const labelVector& ijk,
    const Type& val
)
:
    Field<Type>(cmptProduct(ijk), val),
    ijk_(ijk)
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField
(
    const labelVector& ijk,
    const Foam::zero
)
:
    Field<Type>(cmptProduct(ijk), Foam::zero{}),
    ijk_(ijk)
{}


template<class Type>
inline Foam::IjkField<Type>::IjkField
(
    const labelVector& ijk,
    const UList<Type>& list
)
:
    Field<Type>(list),
    ijk_(ijk)
{
    if (ijk_.size() != Field<Type>::size())
    {
        #ifdef FULLDEBUG
        WarningInFunction
            << "Resizing field to match i-j-k sizing " << sizes()
            << nl << nl;
        #endif

        Field<Type>::resize(ijk_.size(), Zero);
    }
}


template<class Type>
inline Foam::IjkField<Type>::IjkField
(
    const labelVector& ijk,
    Field<Type>&& field
)
:
    Field<Type>(std::move(field)),
    ijk_(ijk)
{
    if (ijk_.size() != Field<Type>::size())
    {
        #ifdef FULLDEBUG
        WarningInFunction
            << "Resizing field to match i-j-k sizing " << sizes()
            << nl << nl;
        #endif

        Field<Type>::resize(ijk_.size(), Zero);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Foam::ijkAddressing& Foam::IjkField<Type>::ijk() const
{
    return ijk_;
}


template<class Type>
inline Foam::ijkAddressing& Foam::IjkField<Type>::ijk()
{
    return ijk_;
}


template<class Type>
inline const Foam::labelVector& Foam::IjkField<Type>::sizes() const
{
    return ijk_.sizes();
}


template<class Type>
inline Foam::labelVector& Foam::IjkField<Type>::sizes()
{
    return ijk_.sizes();
}


template<class Type>
inline const Foam::label& Foam::IjkField<Type>::size
(
    const vector::components cmpt
) const
{
    return ijk_.size(cmpt);
}


template<class Type>
inline void Foam::IjkField<Type>::clear()
{
    ijk_.clear();
    Field<Type>::clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::IjkField<Type>::operator()
(
    const label i,
    const label j,
    const label k
) const
{
    return Field<Type>::operator[](ijk_.index(i, j, k));
}


template<class Type>
inline Type& Foam::IjkField<Type>::operator()
(
    const label i,
    const label j,
    const label k
)
{
    return Field<Type>::operator[](ijk_.index(i, j, k));
}


template<class Type>
inline const Type& Foam::IjkField<Type>::operator()
(
    const labelVector& ijk
) const
{
    return Field<Type>::operator[](ijk_.index(ijk));
}


template<class Type>
inline Type& Foam::IjkField<Type>::operator()
(
    const labelVector& ijk
)
{
    return Field<Type>::operator[](ijk_.index(ijk));
}


template<class Type>
inline void Foam::IjkField<Type>::operator=(IjkField<Type>&& rhs)
{
    if (this != &rhs)
    {
        sizes() = rhs.sizes();

        List<Type>::transfer(rhs);

        rhs.clear();
    }
}


template<class Type>
inline void Foam::IjkField<Type>::operator=(const Type& val)
{
    Field<Type>::operator=(val);
}


template<class Type>
inline void Foam::IjkField<Type>::operator=(const Foam::zero)
{
    Field<Type>::operator=(Foam::zero{});
}


// ************************************************************************* //
