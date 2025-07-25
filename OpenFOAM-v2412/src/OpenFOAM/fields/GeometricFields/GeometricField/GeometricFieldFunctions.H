/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "GeometricScalarField.H"

#define TEMPLATE \
    template<class Type, template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void component
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    const direction d
);

template<class Type, template<class> class PatchField, class GeoMesh>
void T
(
     GeometricField<Type, PatchField, GeoMesh>& result,
     const GeometricField<Type, PatchField, GeoMesh>& f1
);

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
void pow
(
    GeometricField
    <typename powProduct<Type, r>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
);

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
tmp
<
    GeometricField
    <typename powProduct<Type, r>::type, PatchField, GeoMesh>
>
pow
(
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    typename powProduct<Type, r>::type
);

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
tmp
<
    GeometricField
    <typename powProduct<Type, r>::type, PatchField, GeoMesh>
>
pow
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,
    typename powProduct<Type, r>::type
);

template<class Type, template<class> class PatchField, class GeoMesh>
void sqr
(
    GeometricField
    <typename outerProduct<Type, Type>::type, PatchField, GeoMesh>& f1,
    const GeometricField<Type, PatchField, GeoMesh>& gf1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const GeometricField<Type, PatchField, GeoMesh>& f1);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1);

template<class Type, template<class> class PatchField, class GeoMesh>
void magSqr
(
    GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
magSqr
(
    const GeometricField<Type, PatchField, GeoMesh>& f1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
magSqr
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1
);

template<class Type, template<class> class PatchField, class GeoMesh>
void mag
(
    GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& tf1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
mag
(
    const GeometricField<Type, PatchField, GeoMesh>& f1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
mag
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1
);

template<class Type, template<class> class PatchField, class GeoMesh>
void cmptAv
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const GeometricField<Type, PatchField, GeoMesh>& f1);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1);


#define UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(ReturnType, Func, BinaryOp)     \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1                        \
);                                                                             \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1                  \
);

UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, max, maxOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, min, minOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(MinMax<Type>, minMax, minMaxOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(scalarMinMax, minMaxMag, minMaxMagOp)

#undef UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY


#define UNARY_REDUCTION_FUNCTION(ReturnType, Func, gFunc)                      \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1                        \
);                                                                             \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1                  \
);

UNARY_REDUCTION_FUNCTION(Type, sum, gSum)
UNARY_REDUCTION_FUNCTION(Type, average, gAverage)
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

template<class Type, template<class> class PatchField, class GeoMesh>
void clamp
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    const Foam::zero_one
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
clamp
(
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    const Foam::zero_one
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
clamp
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,
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
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>& result,       \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
        <typename product<Type1, Type2>::type, PatchField, GeoMesh>            \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>                \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>                \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>                \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Type, Form>::type, PatchField, GeoMesh>& result,         \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const dimensioned<Form>& dvs                                               \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Type, Form>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const dimensioned<Form>& dvs                                               \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Type, Form>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const dimensioned<Form>& dvs                                               \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>& result,         \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
);                                                                             \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf2                  \
);                                                                             \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp                                                                            \
<                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>                  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf2                  \
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
