/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "exprValue.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::expressions::exprValue::good
(
    const expressions::valueTypeCode valTypeCode
) noexcept
{
    return
    (
        valTypeCode != expressions::valueTypeCode::NONE
     && valTypeCode != expressions::valueTypeCode::INVALID
    );
}


inline bool Foam::expressions::exprValue::is_integral
(
    const expressions::valueTypeCode valTypeCode
) noexcept
{
    return
    (
        valTypeCode == expressions::valueTypeCode::type_bool
     || valTypeCode == expressions::valueTypeCode::type_label
    );
}


inline Foam::direction Foam::expressions::exprValue::nComponents
(
    const expressions::valueTypeCode valTypeCode
) noexcept
{
    return expressions::Detail::nComponents(valTypeCode);
}


inline Foam::direction Foam::expressions::exprValue::rank
(
    const expressions::valueTypeCode valTypeCode
) noexcept
{
    return expressions::Detail::rank(valTypeCode);
}


template<class Type>
inline bool Foam::expressions::exprValue::supportedType()
{
    return
    (
        expressions::valueTypeCode::INVALID != exprTypeTraits<Type>::value
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::expressions::exprValue::fill_zero()
{
    // Preserve the type
    const auto oldType(typeCode_);
    clear();
    typeCode_ = oldType;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::expressions::exprValue::exprValue()
{
    clear();  // zero-initialized (type is NONE)
}


inline Foam::expressions::exprValue::exprValue(const exprValue& rhs)
{
    deepCopy(rhs);
}


template<class Type>
inline Foam::expressions::exprValue::exprValue(const Type& val)
{
    clear();  // type: none, value: zero

    set<Type>(val);

    if (!good())
    {
        FatalErrorInFunction
            << "Cannot construct for unsupported type: "
            << typeid(Type).name() << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::exprValue::good() const noexcept
{
    return exprValue::good(typeCode_);
}


inline bool Foam::expressions::exprValue::is_integral() const noexcept
{
    return exprValue::is_integral(typeCode_);
}


inline Foam::direction
Foam::expressions::exprValue::nComponents() const noexcept
{
    return expressions::Detail::nComponents(typeCode_);
}


inline Foam::direction
Foam::expressions::exprValue::rank() const noexcept
{
    return expressions::Detail::rank(typeCode_);
}


inline Foam::word
Foam::expressions::exprValue::valueTypeName() const
{
    return Foam::name(typeCode_);
}


inline Foam::word Foam::expressions::exprValue::listCompoundName() const
{
    return
    (
        good()
        // List<scalar>, List<vector> ...
      ? word("List<" + Foam::name(typeCode_) + '>')
      : word()
    );
}


template<class Type>
inline bool Foam::expressions::exprValue::set(const Type& val)
{
    const auto whichCode(exprTypeTraits<Type>::value);

    if
    (
        whichCode == expressions::valueTypeCode::NONE
     || whichCode == expressions::valueTypeCode::INVALID
    )
    {
        clear();
        typeCode_ = expressions::valueTypeCode::NONE;
        return false;
    }

    data_.set<Type>(val);
    typeCode_ = whichCode;
    return true;
}


template<class Type>
inline const Type& Foam::expressions::exprValue::get() const
{
    const Type* dataPtr = this->isA<Type>();

    if (dataPtr)
    {
        return *dataPtr;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline const Type* Foam::expressions::exprValue::isA() const noexcept
{
    // NOTE: typeCode_ is *never* set to INVALID,
    // so no false positives for unsupported types
    return
    (
        (typeCode_ == exprTypeTraits<Type>::value)
      ? data_.get<Type>()
      : nullptr
    );
}


// ************************************************************************* //
