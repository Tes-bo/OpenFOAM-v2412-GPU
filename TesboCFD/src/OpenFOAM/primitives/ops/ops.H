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

InNamespace
    Foam

Description
    Various functors for unary and binary operations.
    Can be used for parallel combine-reduce operations or other places
    requiring a functor.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ops_H
#define Foam_ops_H

#include "stdFoam.H"    // For FOAM_NODISCARD

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Assignment operation taking two parameters, returning void.
// Alters the value of the first parameter.
//     Eg, plusEqOp for (x += y)
#define EqOp(opName, op)                                                       \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        void operator()(T1& x, const T2& y) const                              \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        void operator()(T& x, const T& y) const                                \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };


EqOp(eq, x = y)
EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(eqMag, x = mag(y))
EqOp(eqSqr, x = sqr(y))
EqOp(eqMagSqr, x = magSqr(y))
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(minEq, x = min(x, y))
EqOp(maxEq, x = max(x, y))
EqOp(minMagSqrEq, x = (magSqr(x) <= magSqr(y) ? x : y))
EqOp(maxMagSqrEq, x = (magSqr(x) >= magSqr(y) ? x : y))

EqOp(andEq, x = (x && y))
EqOp(orEq,  x = (x || y))
EqOp(xorEq, x = (x != y))
EqOp(bitAndEq, x &= y)
EqOp(bitOrEq,  x |= y)
EqOp(bitXorEq, x ^= y)

EqOp(eqMinus, x = -y)

EqOp(nopEq, (void)x)

#undef EqOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Operation taking two parameters, returning the first type.
// Neither parameter is altered.
//     Eg, plusOp for (x + y)

#define Op(opName, op)                                                         \
                                                                               \
    template<class T, class T1, class T2>                                      \
    struct opName##Op3                                                         \
    {                                                                          \
        FOAM_NODISCARD T operator()(const T1& x, const T2& y) const            \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        FOAM_NODISCARD T1 operator()(const T1& x, const T2& y) const           \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        FOAM_NODISCARD T operator()(const T& x, const T& y) const              \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Operations taking two parameters (unaltered), returning bool

#define BoolOp(opName, op)                                                     \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        FOAM_NODISCARD bool operator()(const T1& x, const T2& y) const         \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        FOAM_NODISCARD bool operator()(const T& x, const T& y) const           \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Operations taking one parameter, returning bool.
// The comparison value is defined during construction

#define Bool1Op(opName, op)                                                    \
                                                                               \
    template<class T>                                                          \
    struct opName##Op1                                                         \
    {                                                                          \
        const T& value;                                                        \
                                                                               \
        opName##Op1(const T& v) : value(v) {}                                  \
                                                                               \
        FOAM_NODISCARD bool operator()(const T& x) const                       \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Weighting operations

#define WeightedOp(opName, op)                                                 \
                                                                               \
    template<class T, class CombineOp>                                         \
    class opName##WeightedOp                                                   \
    {                                                                          \
        const CombineOp& cop_;                                                 \
                                                                               \
    public:                                                                    \
                                                                               \
        opName##WeightedOp(const CombineOp& cop)                               \
        :                                                                      \
            cop_(cop)                                                          \
        {}                                                                     \
                                                                               \
        void operator()                                                        \
        (                                                                      \
            T& x,                                                              \
            const label index,                                                 \
            const T& y,                                                        \
            const scalar weight                                                \
        ) const                                                                \
        {                                                                      \
            cop_(x, op);                                                       \
        }                                                                      \
    };                                                                         \


Op(sum, x + y)

Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptPow, cmptPow(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(min, min(x, y))
Op(max, max(x, y))
Op(minMagSqr, (magSqr(x) <= magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x) >= magSqr(y) ? x : y))
Op(minMod, minMod(x, y))

Op(bitAnd, (x & y))
Op(bitOr,  (x | y))
Op(bitXor, (x ^ y))

BoolOp(and, x && y)
BoolOp(or,  x || y)
BoolOp(xor, (!x) != (!y))   // With forced bool context
BoolOp(equal, x == y)
BoolOp(notEqual, x != y)
BoolOp(less, x < y)
BoolOp(lessEqual, x <= y)
BoolOp(greater, x > y)
BoolOp(greaterEqual, x >= y)

BoolOp(lessEq, x <= y)          // OpenFOAM-v2112 and earlier
BoolOp(greaterEq, x >= y)       // OpenFOAM-v2112 and earlier

Bool1Op(equal, x == value)
Bool1Op(notEqual, x != value)
Bool1Op(less, x < value)
Bool1Op(lessEqual, x <= value)
Bool1Op(greater, x > value)
Bool1Op(greaterEqual, x >= value)

WeightedOp(multiply, (weight*y))

#undef Op
#undef BoolOp
#undef Bool1Op
#undef WeightedOp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Three-way comparison operation of two parameters,
//  similar to the \c <=> operator in C++20.
//
//  \return a negative value for less, a positive value for greater,
//      and zero for equal value.
template<class T>
struct compareOp
{
    FOAM_NODISCARD int operator()(const T& a, const T& b) const
    {
        return (a < b) ? -1 : (b < a) ? 1 : 0;
    }
};


//- Linear interpolation (cf. std::lerp)
template<class T>
struct lerpOp
{
    FOAM_NODISCARD T operator()(const T& a, const T& b, const scalar t) const
    {
        return lerp(a, b, t);
    }
};


//- Linear interpolation (lerp) with interpolation value defined on construct
template<class T>
struct lerpOp1
{
    const scalar value;

    lerpOp1(scalar v) : value(v) {}

    FOAM_NODISCARD T operator()(const T& a, const T& b) const
    {
        return lerp(a, b, value);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Deprecated(2020-11) use nameOp (word.H)
//  \deprecated(2020-11) use nameOp
template<class T> struct getNameOp : nameOp<T> {};

//- Deprecated(2020-11) use typeOp (word.H)
//  \deprecated(2020-11) use typeOp
template<class T> struct getTypeOp : typeOp<T> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
