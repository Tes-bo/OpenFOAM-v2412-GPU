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

#include "DimensionedScalarField.H"

#define TEMPLATE template<class Type, class GeoMesh>
#include "DimensionedFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, direction r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh>>
pow
(
    const DimensionedField<Type, GeoMesh>& f1,
    typename powProduct<Type, r>::type
);

template<class Type, class GeoMesh, direction r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh>>
pow
(
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,
    typename powProduct<Type, r>::type
);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const DimensionedField<Type, GeoMesh>& f1);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const tmp<DimensionedField<Type, GeoMesh>>& tf1);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const DimensionedField<Type, GeoMesh>& f1);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const tmp<DimensionedField<Type, GeoMesh>>& tf1);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const DimensionedField<Type, GeoMesh>& f1);

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const tmp<DimensionedField<Type, GeoMesh>>& tf1);

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType,
        GeoMesh
    >
>
cmptAv(const DimensionedField<Type, GeoMesh>& f1);

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType,
        GeoMesh
    >
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh>>& tf1);


// Forward to FieldFunction via gFunc()
#define UNARY_REDUCTION_FUNCTION(ReturnType, Func, gFunc)                      \
                                                                               \
template<class Type, class GeoMesh>                                            \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1                                  \
);                                                                             \
template<class Type, class GeoMesh>                                            \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1                            \
);

UNARY_REDUCTION_FUNCTION(Type, max, gMax)
UNARY_REDUCTION_FUNCTION(Type, min, gMin)
UNARY_REDUCTION_FUNCTION(Type, sum, gSum)
UNARY_REDUCTION_FUNCTION(Type, average, gAverage)
UNARY_REDUCTION_FUNCTION(MinMax<Type>, minMax, gMinMax)
UNARY_REDUCTION_FUNCTION(scalarMinMax, minMaxMag, gMinMaxMag)

UNARY_REDUCTION_FUNCTION(typename typeOfMag<Type>::type, sumMag, gSumMag)

#undef UNARY_REDUCTION_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// ------------------------------------------------------------------------- //

// Clamp Methods

template<class Type, class GeoMesh>
void clamp
(
    DimensionedField<Type, GeoMesh>& result,
    const DimensionedField<Type, GeoMesh>& f1,
    const Foam::zero_one
);

template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
clamp
(
    const DimensionedField<Type, GeoMesh>& f1,
    const Foam::zero_one
);

template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
clamp
(
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,
    const Foam::zero_one
);

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clamp)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TERNARY_FUNCTION(Type, Type, Type, scalar, lerp)
TERNARY_TYPE_FUNCTION_FFS(Type, Type, Type, scalar, lerp)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate, transform)

BINARY_OPERATOR(Type, Type, scalar, *, '*', multiply)
BINARY_OPERATOR(Type, scalar, Type, *, '*', multiply)
BINARY_OPERATOR(Type, Type, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, '*', multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, '*', multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& f1,                                \
    const DimensionedField<Type2, GeoMesh>& f2                                 \
);                                                                             \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& f1,                                \
    const tmp<DimensionedField<Type2, GeoMesh>>& tf2                           \
);                                                                             \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tf1,                          \
    const DimensionedField<Type2, GeoMesh>& f2                                 \
);                                                                             \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tf1,                          \
    const tmp<DimensionedField<Type2, GeoMesh>>& tf2                           \
);                                                                             \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1,                                 \
    const dimensioned<Form>& dvs                                               \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1,                                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,                           \
    const dimensioned<Form>& dvs                                               \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,                           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const DimensionedField<Type, GeoMesh>& f2                                  \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const DimensionedField<Type, GeoMesh>& f2                                  \
);                                                                             \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf2                            \
);                                                                             \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<DimensionedField<Type, GeoMesh>>& tf2                            \
);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
