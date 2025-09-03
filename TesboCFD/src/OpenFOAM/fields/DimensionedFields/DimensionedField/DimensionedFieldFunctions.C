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

#include "DimensionedFieldReuseFunctions.H"

#define TEMPLATE template<class Type, class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

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
)
{
    typedef typename powProduct<Type, r>::type resultType;

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            f1,
            "pow(" + f1.name() + ',' + Foam::name(r) + ')',
            pow(f1.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tres.ref().field(), f1.field());

    return tres;
}


template<class Type, class GeoMesh, direction r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh>>
pow
(
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            tf1,
            "pow(" + f1.name() + ',' + Foam::name(r) + ')',
            pow(f1.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tres.ref().field(), f1.field());

    tf1.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const DimensionedField<Type, GeoMesh>& f1)
{
    typedef typename outerProduct<Type, Type>::type resultType;

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            f1,
            "sqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    sqr(tres.ref().field(), f1.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const tmp<DimensionedField<Type, GeoMesh>>& tf1)
{
    typedef typename outerProduct<Type, Type>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            tf1,
            "sqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    sqr(tres.ref().field(), f1.field());

    tf1.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const DimensionedField<Type, GeoMesh>& f1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            f1,
            "magSqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    magSqr(tres.ref().field(), f1.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const tmp<DimensionedField<Type, GeoMesh>>& tf1)
{
    typedef typename typeOfMag<Type>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            tf1,
            "magSqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    magSqr(tres.ref().field(), f1.field());

    tf1.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const DimensionedField<Type, GeoMesh>& f1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            f1,
            "mag(" + f1.name() + ')',
            f1.dimensions()
        );

    mag(tres.ref().field(), f1.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const tmp<DimensionedField<Type, GeoMesh>>& tf1)
{
    typedef typename typeOfMag<Type>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            tf1,
            "mag(" + f1.name() + ')',
            f1.dimensions()
        );

    mag(tres.ref().field(), f1.field());

    tf1.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >
>
cmptAv(const DimensionedField<Type, GeoMesh>& f1)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType resultType;

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            f1,
            "cmptAv(" + f1.name() + ')',
            f1.dimensions()
        );

    cmptAv(tres.ref().field(), f1.field());

    return tres;
}


template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh>>& tf1)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New
        (
            tf1,
            "cmptAv(" + f1.name() + ')',
            f1.dimensions()
        );

    cmptAv(tres.ref().field(), f1.field());

    tf1.clear();
    return tres;
}


#define UNARY_REDUCTION_FUNCTION(ReturnType, Func, gFunc)                      \
                                                                               \
template<class Type, class GeoMesh>                                            \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1                                  \
)                                                                              \
{                                                                              \
    return dimensioned<ReturnType>                                             \
    (                                                                          \
        #Func "(" + f1.name() + ')',                                           \
        f1.dimensions(),                                                       \
        gFunc(f1.field())                                                      \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1                            \
)                                                                              \
{                                                                              \
    dimensioned<ReturnType> res = Func(tf1());                                 \
    tf1.clear();                                                               \
    return res;                                                                \
}

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
)
{
    const MinMax<Type> range(Foam::zero_one{});

    clamp(result.field(), f1.field(), range);
    result.oriented() = f1.oriented();
}

template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
clamp
(
    const DimensionedField<Type, GeoMesh>& f1,
    const Foam::zero_one
)
{
    auto tres =
        reuseTmpDimensionedField<Type, Type, GeoMesh>::New
        (
            f1,
            "clamp01(" + f1.name() + ')',
            f1.dimensions()
        );

    clamp(tres.ref(), f1, Foam::zero_one{});

    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
clamp
(
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,
    const Foam::zero_one
)
{
    const auto& f1 = tf1();

    auto tres =
        reuseTmpDimensionedField<Type, Type, GeoMesh>::New
        (
            tf1,
            "clamp01(" + f1.name() + ')',
            f1.dimensions()
        );

    clamp(tres.field(), f1, Foam::zero_one{});

    tf1.clear();
    return tres;
}

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
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type1, GeoMesh>::New              \
        (                                                                      \
            f1,                                                                \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), f2.field());                  \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& f1,                                \
    const tmp<DimensionedField<Type2, GeoMesh>>& tf2                           \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type2, GeoMesh>::New              \
        (                                                                      \
            tf2,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), f2.field());                  \
                                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tf1,                          \
    const DimensionedField<Type2, GeoMesh>& f2                                 \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type1, GeoMesh>::New              \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), f2.field());                  \
                                                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tf1,                          \
    const tmp<DimensionedField<Type2, GeoMesh>>& tf2                           \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpDimensionedField                                            \
        <resultType, Type1, Type1, Type2, GeoMesh>::New                        \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), f2.field());                  \
                                                                               \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1,                                 \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New               \
        (                                                                      \
            f1,                                                                \
            '(' + f1.name() + #Op + dvs.name() + ')',                          \
            (f1.dimensions() Op dvs.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), dvs.value());                 \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& f1,                                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return f1 Op dimensioned<Form>(static_cast<const Form&>(vs));              \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,                           \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
                                                                               \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New               \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + #Op + dvs.name() + ')',                          \
            (f1.dimensions() Op dvs.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), f1.field(), dvs.value());                 \
                                                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf1,                           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return tf1 Op dimensioned<Form>(static_cast<const Form&>(vs));             \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const DimensionedField<Type, GeoMesh>& f2                                  \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New               \
        (                                                                      \
            f2,                                                                \
            '(' + dvs.name() + #Op + f2.name() + ')',                          \
            (dvs.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), dvs.value(), f2.field());                 \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const DimensionedField<Type, GeoMesh>& f2                                  \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) Op f2;              \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tf2                            \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
                                                                               \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<resultType, Type, GeoMesh>::New               \
        (                                                                      \
            tf2,                                                               \
            '(' + dvs.name() + #Op + f2.name() + ')',                          \
            (dvs.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref().field(), dvs.value(), f2.field());                 \
                                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<DimensionedField<Type, GeoMesh>>& tf2                            \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) Op tf2;             \
}


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
