/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application

Description
    Tests for complex numbers

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "complex.H"
#include "complexFields.H"
#include "scalarField.H"
#include "diagTensor.H"
#include "symmTensor.H"
#include "symmTensor2D.H"
#include "ListOps.H"
#include "ops.H"

using namespace Foam;


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< "complex()     : " << complex() << nl
        << "complex(zero) : " << complex(Zero) << nl
        << "pTraits<complex>::zero : " << pTraits<complex>::zero << nl
        << "pTraits<complex>::one  : " << pTraits<complex>::one << nl
        << "complex(scalar) : " << complex(3.14519) << nl
        << nl;

    std::complex<scalar> c1(10, -3);
    Info<< "std::complex : " << c1 << nl;
    Info<< "sin: " << std::sin(c1) << nl;

    Info<< "complexVector::zero : " << complexVector::zero << nl
        << "complexVector::one  : " << complexVector::one << nl
        << nl;

    {
        const complex a(0, 1);
        const complex b(20, 100);

        Info<< "lerp of " << a << " : " << b << nl;

        for (const double t : { 0.0, 0.5, 1.0, -0.5, 1.5 })
        {
            Info<< "    " << t << " = " << lerp(a, b, t) << nl;
        }
    }

    for (complex c : { complex{1, 0}, complex{1, 2}} )
    {
        Info<< nl << "complex : " << c
            << " mag = " << c.magnitude()
            << " norm = " << c.magSqr()
            << nl;

        Info<< "sin: " << sin(c) << nl
            << "pow(3.0f): " << pow(c, 3.0f) << nl
            << "pow(3): " << pow(c, 3) << nl
            << "pow3: " << pow3(c) << nl
            << "log: " << log(c) << nl
            << "pow025: " << pow025(c) << nl
            ;

        // TDB: allow implicit construct from scalar?
        //
        // if (c == 1.0)
        // {
        //     Info<< c << " == " << 1 << nl;
        // }
    }

    // Test powers of zero
    #if 1
    {
        const complex complex0{0, 0};
        const std::complex<scalar> std0{0, 0};

        const label label0{0};
        const scalar scalar0{0};

        Info<< nl
            << "# std::pow(0, 0)" << nl
            << "    (label, label) = " << std::pow(label0, label0) << nl
            << "    (scalar, scalar) = " << std::pow(scalar0, scalar0) << nl
            << "    (label, scalar) = " << std::pow(label0, scalar0) << nl
            << "    (scalar, label) = " << std::pow(scalar0, label0) << nl
            << "    (std::complex, label) = " << std::pow(std0, label0) << nl
            << "    (std::complex, scalar) = " << std::pow(std0, scalar0) << nl
            << "    (label, std::complex) = " << std::pow(label0, std0) << nl
            << "    (scalar, std::complex) = " << std::pow(scalar0, std0) << nl
            ;

        Info<< nl
            << "# Foam::pow(0, 0)" << nl
            << "    (label, label) = " << Foam::pow(label0, label0) << nl
            << "    (scalar, scalar) = " << Foam::pow(scalar0, scalar0) << nl
            << "    (label, scalar) = " << Foam::pow(label0, scalar0) << nl
            << "    (scalar, label) = " << Foam::pow(scalar0, label0) << nl
            << "    (complex, label) = " << Foam::pow(complex0, label0) << nl
            << "    (complex, scalar) = " << Foam::pow(complex0, scalar0) << nl
            << "    (label, complex) = " << Foam::pow(label0, complex0) << nl
            << "    (scalar, complex) = " << Foam::pow(scalar0, complex0) << nl
            ;
    }
    #endif

    // Test zip/unzip
    {
        scalarField reals(4);
        scalarField imags(4);

        forAll(reals, i)
        {
            reals[i] = i;
        }
        forAll(imags, i)
        {
            imags[i] = (i % 2) ? -i : i;
        }

        complexField cmplx(4);

        zip(cmplx, reals, zero{});
        zip(cmplx, 1, imags);
        zip(cmplx, reals, imags);
        Info<< nl
            << "zip " << reals << nl
            << "    " << imags << nl
            << " => " << cmplx << nl;

        reverse(cmplx);

        Info<< "reverse order: " << cmplx << nl;

        unzip(cmplx, reals, imags);

        Info<< "unzip " << cmplx << nl
            << " => " << reals << nl
            << " => " << imags << nl;
    }

    {
        SymmTensor<complex> st1(SymmTensor<complex>::uniform({3, 4}));

        Info<< "symmTensor: " << st1 << nl
            << "  tr: " << tr(st1) << nl
            << "  diagSqr: " << st1.diagSqr() << nl
            << "  magSqr: " << magSqr(st1) << nl
            << "  mag: " << mag(st1) << nl;

        SymmTensor<scalar> st2(SymmTensor<scalar>::uniform(5));

        Info<< "symmTensor: " << st2 << nl
            << "  tr: " << tr(st2) << nl
            << "  diagSqr: " << st2.diagSqr() << nl
            << "  magSqr: " << magSqr(st2) << nl
            << "  mag: " << mag(st2) << nl;

        st2 = Zero;

        DiagTensor<complex> dt1(SphericalTensor<complex>({3, 4}));

        Info<< "diagTensor: " << dt1 << nl
            << "  tr: " << tr(dt1) << nl
            << "  diagSqr: " << dt1.diagSqr() << nl
            << "  magSqr: " << magSqr(dt1) << nl
            << "  mag: " << mag(dt1) << nl;


        // A bit ugly...
        st1 = SphericalTensor<complex>({3, 4});

        Info<< "symmTensor: " << st1 << nl
            << "  tr: " << tr(st1) << nl
            << "  diagSqr: " << st1.diagSqr() << nl
            << "  magSqr: " << magSqr(st1) << nl
            << "  mag: " << mag(st1) << nl;
    }

    {
        SymmTensor2D<complex> st1(SymmTensor2D<complex>::uniform({3, 4}));

        Info<< "symmTensor: " << st1 << nl
            << "  tr: " << tr(st1) << nl
            << "  diagSqr: " << st1.diagSqr() << nl
            << "  magSqr: " << magSqr(st1) << nl
            << "  mag: " << mag(st1) << nl;
    }

    {
        Tensor<complex> st1(Tensor<complex>::uniform({3, 4}));

        Info<< "tensor: " << st1 << nl
            << "  tr: " << tr(st1) << nl
            << "  diagSqr: " << st1.diagSqr() << nl
            << "  magSqr: " << magSqr(st1) << nl
            << "  mag: " << mag(st1) << endl;

        Tensor<scalar> st2(Tensor<scalar>::uniform(5));

        Info<< "Tensor: " << st2 << nl
            << "  tr: " << tr(st2) << nl
            << "  diagSqr: " << st2.diagSqr() << nl
            << "  magSqr: " << magSqr(st2) << nl
            << "  mag: " << mag(st2) << endl;
    }

    complexField fld1(3, complex(2.0, 1.0));
    complexField fld2(fld1);

    for (complex& c : fld2)
    {
        c = ~c;
    }

    Info<< nl
        << "Field " << flatOutput(fld1) << nl
        << "Conjugate: " << flatOutput(fld2) << nl;

    // Some arbitrary change
    for (complex& c : fld2)
    {
        c.Im() *= 5;
    }

    Info<< "sumProd: " << sumProd(fld1, fld2) << nl;

    fld1 *= 10;
    Info<< "scalar multiply: " << flatOutput(fld1) << nl;

    fld1 /= 10;
    Info<< "scalar divide: " << flatOutput(fld1) << nl;

    Info<< "sin: " << sin(fld1) << nl;

    Info<< "operator + : " << (fld1 + fld2) << nl;
    Info<< "operator + : " << (fld1 + fld2 + complex(1,0)) << nl;

    // Some operators are still incomplete

    // Info<< "operator * : " << (fld1 * fld2) << nl;
    // Info<< "operator / : " << (fld1 / fld2) << nl;
    // Info<< "operator / : " << (fld1 / 2) << nl;
    // Info<< "operator / : " << (fld1 / fld2) << nl;
    // Info<< "sqrt   : " << sqrt(fld1) << nl;
    // Info<< "pow(2) : " << pow(fld1, 2) << nl;

    #if 1
    Info<< nl << "## Elementary complex-complex arithmetic operations:" << nl;
    {
        const complex a(6, 1);
        complex b = a;

        Info<< "# Compound assignment operations:" << nl;

        Info<< "a = " << a << ", b = " << b << nl;

        // Addition
        b += a;
        Info<< "b += a:" << tab << "b =" << b << nl;

        // Subtraction
        b -= a;
        Info<< "b -= a:" << tab << "b =" << b << nl;

        // Multiplication
        b *= a;
        Info<< "b *= a:" << tab << "b =" << b << nl;

        // Division
        b /= a;
        Info<< "b /= a:" << tab << "b =" << b << nl;
    }
    #endif


    #if 1
    Info<< nl << "## Elementary complex-scalar arithmetic operations:" << nl;
    {
        const scalar a = 5;
        complex b(6, 1);

        Info<< "# Non-assignment operations:" << nl;

        Info<< "(scalar) a = " << a << ", b = " << b << nl;

        // Addition
        b = a + b;
        Info<< "b = a + b: " << tab << b << nl;

        b = b + a;
        Info<< "b = b + a: " << tab << b << nl;

        // Subtraction
        b = a - b;
        Info<< "b = a - b: " << tab << b << nl;

        b = b - a;
        Info<< "b = b - a: " << tab << b << nl;

        // Multiplication
        b = a*b;
        Info<< "b = a*b: " << tab << b << nl;

        b = b*a;
        Info<< "b = b*a: " << tab << b << nl;

        // Division
        b = a/b;
        Info<< "b = a/b = scalar(a)/b = complex(a)/b:" << tab << b << nl;

        b = b/a;
        Info<< "b = b/a: " << tab << b << nl;


        Info<< "# Compound assignment operations:" << nl;

        Info<< "(scalar) a = " << a << ", b = " << b << nl;

        // Addition: complex+scalar
        b += a;
        Info<< "b += a (only real part):" << tab << b << nl;

        // Subtraction: complex-scalar
        b -= a;
        Info<< "b -= a (only real part):" << tab << b << nl;

        // Multiplication: complex*scalar
        b *= a;
        Info<< "b *= a (real and imag parts):" << tab << b << nl;

        // Division: complex/scalar
        b /= a;
        Info<< "b /= a (real and imag parts):" << tab << b << nl;
    }
    #endif


    #if 1
    Info<< nl << "## Other mathematical expressions:" << nl;
    {
        const complex a(4.3, -3.14);
        const complex b(0, -4.3);
        const complex c(-4.3, 0);

        Info<< "a = " << a << ", b = " << b << ", c = " << c << nl;

        // Square-root
        Info<< "sqrt(a) = " << Foam::sqrt(a) << ", "
            << "sqrt(b) = " << Foam::sqrt(b) << ", "
            << "sqrt(c) = " << Foam::sqrt(c) << nl;

        // Square
        Info<< "sqr(a) = " << sqr(a) << ", "
            << "sqr(b) = " << sqr(b) << ", "
            << "sqr(c) = " << sqr(c) << nl;

        // n^th power
        Info<< "pow(a, -1) = " << pow(a, -1) << ", "
            << "pow(b, -1) = " << pow(b, -1) << ", "
            << "pow(c, -1) = " << pow(c, -1) << nl;

        // Exponential
        Info<< "exp(a) = " << exp(a) << ", "
            << "exp(b) = " << exp(b) << ", "
            << "exp(c) = " << exp(c) << nl;

        // Natural logarithm
        Info<< "log(a) = " << log(a) << ", "
            << "log(b) = " << log(b) << ", "
            << "log(c) = " << log(c) << nl;
    }
    #endif


    // Make some changes
    {
        label i = 1;
        for (complex& c : fld1)
        {
            c.Re() += i;
            c.Im() -= 10 - i;
            ++i;
        }
    }

    Info<< nl
        << "field = " << fld1 << nl;

    Info<< "magSqr = "
        << ListOps::create<scalar>
           (
               fld1,
               [](const complex& c) { return magSqr(c); }
           )
        << nl;

    Info
        << "sum = " << sum(fld1) << nl
        << "min = " << min(fld1) << nl
        << "max = " << max(fld1) << nl;

    // MinMax fails since there is no less comparison operator
    // Info<< "min/max = " << MinMax<complex>(fld1) << nl;


    // Cross-product
    {
        const vector vec(1, 2, 3);
        const vector realValue(4, 5, 6);
        const vector imagValue(7, 8, 9);

        complexVector cmplxVec(zip(realValue, imagValue));

        Info<< "complexVector: " << cmplxVec << nl;
        Info<< "cross: " << (vec ^ cmplxVec) << nl;

        Info<< "cross real: " << (vec ^ realValue) << nl
            << "cross imag: " << (vec ^ imagValue) << nl
            << "cross     : "
            << zip((vec ^ realValue), (vec ^ imagValue)) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
