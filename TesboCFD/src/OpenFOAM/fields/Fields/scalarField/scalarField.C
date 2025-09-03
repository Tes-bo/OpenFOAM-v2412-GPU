/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Specialisation of Field\<T\> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "unitConversion.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const
{
    return *this;
}

void component(scalarField& sf, const UList<scalar>& f, const direction)
{
    sf = f;
}

template<>
void scalarField::replace(const direction, const UList<scalar>& sf)
{
    *this = sf;
}

template<>
void scalarField::replace(const direction, const scalar& s)
{
    *this = s;
}


void stabilise(scalarField& res, const UList<scalar>& sf, const scalar s)
{
    // std::transform
    TFOR_ALL_F_OP_FUNC_S_F
    (
        scalar, res, =, ::Foam::stabilise, scalar, s, scalar, sf
    )
}

tmp<scalarField> stabilise(const UList<scalar>& sf, const scalar s)
{
    auto tresult = tmp<scalarField>::New(sf.size());
    stabilise(tresult.ref(), sf, s);
    return tresult;
}

tmp<scalarField> stabilise(const tmp<scalarField>& tsf, const scalar s)
{
    tmp<scalarField> tresult = New(tsf);
    stabilise(tresult.ref(), tsf(), s);
    tsf.clear();
    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
float sumProd(const UList<float>& f1, const UList<float>& f2)
{
    float result = 0.0;
    if (f1.size() && (f1.size() == f2.size()))
    {
        // std::inner_product
        TFOR_ALL_S_OP_F_OP_F(float, result, +=, float, f1, *, float, f2)
    }
    return result;
}


template<>
double sumProd(const UList<double>& f1, const UList<double>& f2)
{
    double result = 0.0;
    if (f1.size() && (f1.size() == f2.size()))
    {
        // std::inner_product
        TFOR_ALL_S_OP_F_OP_F(double, result, +=, double, f1, *, double, f2)
    }
    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, divide)

BINARY_FUNCTION(scalar, scalar, scalar, pow)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, pow)

BINARY_FUNCTION(scalar, scalar, scalar, atan2)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, atan2)

BINARY_FUNCTION(scalar, scalar, scalar, hypot)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, hypot)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3)
UNARY_FUNCTION(scalar, scalar, pow4)
UNARY_FUNCTION(scalar, scalar, pow5)
UNARY_FUNCTION(scalar, scalar, pow6)
UNARY_FUNCTION(scalar, scalar, pow025)
UNARY_FUNCTION(scalar, scalar, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt)
UNARY_FUNCTION(scalar, scalar, sign)
UNARY_FUNCTION(scalar, scalar, pos)
UNARY_FUNCTION(scalar, scalar, pos0)
UNARY_FUNCTION(scalar, scalar, neg)
UNARY_FUNCTION(scalar, scalar, neg0)
UNARY_FUNCTION(scalar, scalar, posPart)
UNARY_FUNCTION(scalar, scalar, negPart)
UNARY_FUNCTION(scalar, scalar, exp)
UNARY_FUNCTION(scalar, scalar, log)
UNARY_FUNCTION(scalar, scalar, log10)
UNARY_FUNCTION(scalar, scalar, sin)
UNARY_FUNCTION(scalar, scalar, cos)
UNARY_FUNCTION(scalar, scalar, tan)
UNARY_FUNCTION(scalar, scalar, asin)
UNARY_FUNCTION(scalar, scalar, acos)
UNARY_FUNCTION(scalar, scalar, atan)
UNARY_FUNCTION(scalar, scalar, sinh)
UNARY_FUNCTION(scalar, scalar, cosh)
UNARY_FUNCTION(scalar, scalar, tanh)
UNARY_FUNCTION(scalar, scalar, asinh)
UNARY_FUNCTION(scalar, scalar, acosh)
UNARY_FUNCTION(scalar, scalar, atanh)
UNARY_FUNCTION(scalar, scalar, erf)
UNARY_FUNCTION(scalar, scalar, erfc)
UNARY_FUNCTION(scalar, scalar, lgamma)
UNARY_FUNCTION(scalar, scalar, j0)
UNARY_FUNCTION(scalar, scalar, j1)
UNARY_FUNCTION(scalar, scalar, y0)
UNARY_FUNCTION(scalar, scalar, y1)

UNARY_FUNCTION(scalar, scalar, degToRad)
UNARY_FUNCTION(scalar, scalar, radToDeg)
UNARY_FUNCTION(scalar, scalar, atmToPa)
UNARY_FUNCTION(scalar, scalar, barToPa)
UNARY_FUNCTION(scalar, scalar, paToAtm)
UNARY_FUNCTION(scalar, scalar, paToBar)


#define BesselFunc(func)                                                       \
void func(scalarField& res, const int n, const UList<scalar>& sf)              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_S_F(scalar, res, =, ::Foam::func, int, n, scalar, sf)   \
}                                                                              \
                                                                               \
tmp<scalarField> func(const int n, const UList<scalar>& sf)                    \
{                                                                              \
    auto tresult = tmp<scalarField>::New(sf.size());                           \
    func(tresult.ref(), n, sf);                                                \
    return tresult;                                                            \
}                                                                              \
                                                                               \
tmp<scalarField> func(const int n, const tmp<scalarField>& tsf)                \
{                                                                              \
    tmp<scalarField> tresult = New(tsf);                                       \
    func(tresult.ref(), n, tsf());                                             \
    tsf.clear();                                                               \
    return tresult;                                                            \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
