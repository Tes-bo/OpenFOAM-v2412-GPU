/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "complexField.H"
#include "addToRunTimeSelectionTable.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<complex>, complexList);
    addCompoundToRunTimeSelectionTable(List<complex>, complexList);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::zip
(
    complexField& result,
    const UList<scalar>& realValues,
    const UList<scalar>& imagValues
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != realValues.size() || len != imagValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << realValues.size() << ' ' << imagValues.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i].real(realValues[i]);
        result[i].imag(imagValues[i]);
    }
}


void Foam::zip
(
    complexField& result,
    const UList<scalar>& realValues,
    const scalar imagValue
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != realValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len
            << " != " << realValues.size() << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i].real(realValues[i]);
        result[i].imag(imagValue);
    }
}


void Foam::zip
(
    complexField& result,
    const scalar realValue,
    const UList<scalar>& imagValues
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != imagValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len
            << " != " << imagValues.size() << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i].real(realValue);
        result[i].imag(imagValues[i]);
    }
}


void Foam::unzip
(
    const UList<complex>& input,
    scalarField& realValues,
    scalarField& imagValues
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != realValues.size() || len != imagValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << realValues.size() << ' ' << imagValues.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        realValues[i] = input[i].real();
        imagValues[i] = input[i].imag();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::complexField Foam::ComplexField
(
    const UList<scalar>& realValues,
    const UList<scalar>& imagValues
)
{
    complexField result(realValues.size());

    Foam::zip(result, realValues, imagValues);

    return result;
}


Foam::complexField Foam::ComplexField
(
    const UList<scalar>& realValues,
    const scalar imagValue
)
{
    complexField result(realValues.size());

    Foam::zip(result, realValues, imagValue);

    return result;
}


Foam::complexField Foam::ComplexField
(
    const scalar realValue,
    const UList<scalar>& imagValues
)
{
    complexField result(imagValues.size());

    Foam::zip(result, realValue, imagValues);

    return result;
}


Foam::scalarField Foam::ReImSum(const UList<complex>& cmplx)
{
    scalarField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complex& c) { return c.cmptSum(); }
    );

    return result;
}


Foam::scalarField Foam::Re(const UList<complex>& cmplx)
{
    scalarField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complex& c) { return c.real(); }
    );

    return result;
}


Foam::scalarField Foam::Im(const UList<complex>& cmplx)
{
    scalarField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complex& c) { return c.imag(); }
    );

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
complex sumProd(const UList<complex>& f1, const UList<complex>& f2)
{
    complex result(0);
    if (f1.size() && (f1.size() == f2.size()))
    {
        // std::inner_product
        TFOR_ALL_S_OP_F_OP_F(complex, result, +=, complex, f1, *, complex, f2)
    }
    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(complex, complex, complex, +, add)
BINARY_TYPE_OPERATOR(complex, complex, complex, -, subtract)

BINARY_OPERATOR(complex, complex, complex, *, multiply)
BINARY_OPERATOR(complex, complex, complex, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(complex, complex, pow3)
UNARY_FUNCTION(complex, complex, pow4)
UNARY_FUNCTION(complex, complex, pow5)
UNARY_FUNCTION(complex, complex, pow6)
UNARY_FUNCTION(complex, complex, pow025)
UNARY_FUNCTION(complex, complex, sqrt)
UNARY_FUNCTION(complex, complex, exp)
UNARY_FUNCTION(complex, complex, log)
UNARY_FUNCTION(complex, complex, log10)
UNARY_FUNCTION(complex, complex, sin)
UNARY_FUNCTION(complex, complex, cos)
UNARY_FUNCTION(complex, complex, tan)
UNARY_FUNCTION(complex, complex, asin)
UNARY_FUNCTION(complex, complex, acos)
UNARY_FUNCTION(complex, complex, atan)
UNARY_FUNCTION(complex, complex, sinh)
UNARY_FUNCTION(complex, complex, cosh)
UNARY_FUNCTION(complex, complex, tanh)
UNARY_FUNCTION(complex, complex, asinh)
UNARY_FUNCTION(complex, complex, acosh)
UNARY_FUNCTION(complex, complex, atanh)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
