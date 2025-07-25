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

#include "exprTraits.H"

//TBD: handle complex?

#undef  FOR_ALL_EXPR_TYPE_CODES
#define FOR_ALL_EXPR_TYPE_CODES(Macro, ...)                                   \
    Macro(bool, __VA_ARGS__)                                                  \
    Macro(label, __VA_ARGS__)                                                 \
    Macro(scalar, __VA_ARGS__)                                                \
    Macro(vector, __VA_ARGS__)                                                \
    Macro(sphericalTensor, __VA_ARGS__)                                       \
    Macro(symmTensor, __VA_ARGS__)                                            \
    Macro(tensor, __VA_ARGS__)


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::direction
Foam::expressions::Detail::nComponents
(
    const expressions::valueTypeCode typeCode
) noexcept
{
    switch (typeCode)
    {
        case expressions::valueTypeCode::NONE :
        case expressions::valueTypeCode::INVALID :
        {
            break;
        }

        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            return exprTypeTraits<Type>::nComponents;                 \
        }

        FOR_ALL_EXPR_TYPE_CODES(doLocalCode);
        #undef doLocalCode
    }

    return 0;
}


Foam::direction
Foam::expressions::Detail::rank
(
    const expressions::valueTypeCode typeCode
) noexcept
{
    switch (typeCode)
    {
        case expressions::valueTypeCode::NONE :
        case expressions::valueTypeCode::INVALID :
        {
            break;
        }

        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            return exprTypeTraits<Type>::rank;                        \
        }

        FOR_ALL_EXPR_TYPE_CODES(doLocalCode);
        #undef doLocalCode
    }

    return 0;
}


Foam::expressions::valueTypeCode
Foam::expressions::valueTypeCodeOf
(
    const word& dataTypeName,
    const expressions::valueTypeCode deflt
)
{
    if (!dataTypeName.empty())
    {
        // Could compare with pTraits<Type>::typeName instead of
        // exprTypeTraits<Type>::name, but then we might miss
        // possible typos.

        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        if (dataTypeName == exprTypeTraits<Type>::name)               \
        {                                                             \
            return expressions::valueTypeCode::type_##Type;           \
        }

        FOR_ALL_EXPR_TYPE_CODES(doLocalCode);
        #undef doLocalCode
    }

    return deflt;
}


Foam::word Foam::name(const expressions::valueTypeCode typeCode)
{
    switch (typeCode)
    {
        case expressions::valueTypeCode::NONE :
        {
            return "none";
        }

        case expressions::valueTypeCode::INVALID :
        {
            // returns "", could also return "bad"
            break;
        }

        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            return exprTypeTraits<Type>::name;                        \
        }

        FOR_ALL_EXPR_TYPE_CODES(doLocalCode);
        #undef doLocalCode
    }

    return word();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  FOR_ALL_EXPR_TYPE_CODES

// ************************************************************************* //
