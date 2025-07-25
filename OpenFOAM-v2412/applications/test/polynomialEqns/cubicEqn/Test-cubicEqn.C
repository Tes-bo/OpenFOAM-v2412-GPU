/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Test-cubicEqn

Description
    Tests for \c cubicEqn constructors, and member functions.

\*---------------------------------------------------------------------------*/

#include <ctime>
#include <random>

#include "cubicEqn.H"
#include "vector.H"
#include "scalarList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Compare two floating point types, and print output.
// Do ++nFail_ if values of two objects are not equal within a given tolerance.
// The function is converted from PEP-485.
void cmp
(
    const word& msg,
    const scalar& x,
    const scalar& y,
    const scalar relTol = 1e-8,  //<! are values the same within 8 decimals
    const scalar absTol = 0      //<! useful for cmps near zero
)
{
    Info<< msg << x << endl;

    unsigned nFail = 0;

    if (max(absTol, relTol*max(mag(x), mag(y))) < mag(x - y))
    {
        ++nFail;
    }

    if (nFail)
    {
        Info<< nl
            << "        #### Fail in " << nFail << " comps ####" << nl << endl;
        ++nFail_;
    }
    ++nTest_;
}


// Create each constructor of cubicEqn, and print output
void test_constructors()
{
    #if 0
    {
        Info<< "# Construct null:" << nl;
        const cubicEqn T();
        Info<< T << endl;
    }
    #endif

    {
        Info<< "# Construct initialized to zero:" << nl;
        const cubicEqn T(Zero);
        Info<< T << endl;
    }
    {
        Info<< "# Construct from components:" << nl;
        const cubicEqn T(1, 2, 3, 4);
        Info<< T << endl;
    }
}


// Execute each member function of cubicEqn, and print output
void test_member_funcs()
{
    const cubicEqn T(1, -6, -2, 12);
    const scalar x = 9.7;

    Info<< "# Operands: " << nl
        << "  cubicEqn = " << T << nl
        << "  x = " << x << endl;


    {
        Info<< "# Access:" << nl;
        cmp("  a = ",  1, T.a());
        cmp("  b = ", -6, T.b());
        cmp("  c = ", -2, T.c());
        cmp("  d = ", 12, T.d());
    }
    {
        Info<< "# Evaluate:" << nl;
        cmp("  T.value(x) = ", T.value(x), 340.73299999999983);
        cmp("  T.derivative(x) = ", T.derivative(x), 163.87);
        cmp
        (
            "  T.error(x) = ",
            T.error(x),
            2.1854789999999994e-13,
            1e-8,
            1e-8
        );

        const vector rts(6, 1.41421356, -1.41421356);
        forAll(rts, i)
        {
            cmp("  T.roots() = ", T.roots()[i], rts[i]);
        }
    }

    {
        Info<< "# Special cases for the root evaluation: " << nl;

        {
            Info<< "  Three distinct real roots:" << nl;
            const cubicEqn Tc(1.0, -4.0, 1.0, 6.0);
            const vector rts(3, 2, -1);
            forAll(rts, i)
            {
                cmp("    root = ", Tc.roots()[i], rts[i]);
            }
        }
        {
            Info<< "  One real root and one complex conjugate-pair root:" << nl;
            const cubicEqn Tc(1.0, 2.0, 3.0, 4.0);
            const vector rts(-1.65062919, -0.1746854, 1.54686889);
            forAll(rts, i)
            {
                cmp("    root = ", Tc.roots()[i], rts[i], 1e-8, 1e-8);
            }
        }
        {
            Info<< "  Two identical real roots and a distinct real root:" << nl;
            const cubicEqn Tc(1.0, -5.0, 8.0, -4.0);
            const vector rts(2, 2, 1);
            forAll(rts, i)
            {
                cmp("    root = ", Tc.roots()[i], rts[i]);
            }
        }
        {
            Info<< "  Three identical real roots:" << nl;
            const cubicEqn Tc(1.0, -6.0, 12.0, -8.0);
            const vector rts(2, 2, 2);
            forAll(rts, i)
            {
                cmp("    root = ", Tc.roots()[i], rts[i]);
            }
        }
        {
            Info<< "  Arbitrary roots with non-unity leading term:" << nl;
            const cubicEqn Tc(-8.59, 6.77, -0.2, 8.0);
            const vector rts(1.31167947, -0.26177687, 0.80093099);
            forAll(rts, i)
            {
                cmp("    root = ", Tc.roots()[i], rts[i]);
            }
        }
    }
}


scalar randomScalar(const scalar min, const scalar max)
{
    static_assert
    (
        sizeof(long) == sizeof(scalar),
        "Scalar and long are not the same size"
    );
    static std::default_random_engine generator(std::time(0));
    static std::uniform_int_distribution<long>
        distribution
        (
            std::numeric_limits<long>::min(),
            std::numeric_limits<long>::max()
        );

    scalar x = 0;

    do
    {
        x = scalar(distribution(generator));
    }
    while (min > mag(x) || mag(x) > max || !std::isfinite(x));

    return x;
};


template <class Type>
void test(const Type& polynomialEqn, const scalar tol)
{
    Roots<Type::nComponents - 1> r = polynomialEqn.roots();

    const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
    const scalar inf = std::numeric_limits<scalar>::infinity();

    FixedList<label, Type::nComponents - 1> t;
    FixedList<scalar, Type::nComponents - 1> v(nan);
    FixedList<scalar, Type::nComponents - 1> e(nan);
    bool ok = true;
    forAll(r, i)
    {
        t[i] = r.type(i);
        switch (t[i])
        {
            case roots::real:
                v[i] = polynomialEqn.value(r[i]);
                e[i] = polynomialEqn.error(r[i]);
                ok = ok && mag(v[i]) <= tol*mag(e[i]);
                break;
            case roots::posInf:
                v[i] = + inf;
                e[i] = nan;
                break;
            case roots::negInf:
                v[i] = - inf;
                e[i] = nan;
                break;
            default:
                v[i] = e[i] = nan;
                break;
        }
    }

    if (!ok)
    {
        Info<< "Coeffs: " << polynomialEqn << endl
            << " Types: " << t << endl
            << " Roots: " << r << endl
            << "Values: " << v << endl
            << "Errors: " << e << endl << endl;
        ++nFail_;
    }
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    Info<< nl << "    ## Test constructors: ##" << nl;
    test_constructors();

    Info<< nl << "    ## Test member functions: ##" << nl;
    test_member_funcs();


    // Pre-v2006 tests
    const scalar tol = 5;

    const label nTests = 1000000;
    for (label t = 0; t < nTests; ++ t)
    {
        test
        (
            cubicEqn
            (
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50)
            ),
            tol
        );
        ++nTest_;
    }
    Info << nTests << " scalar cubic equations were tested." << endl;

    const label coeffMin = -9, coeffMax = 10, nCoeff = coeffMax - coeffMin;
    for (label a = coeffMin; a < coeffMax; ++ a)
    {
        for (label b = coeffMin; b < coeffMax; ++ b)
        {
            for (label c = coeffMin; c < coeffMax; ++ c)
            {
                for (label d = coeffMin; d < coeffMax; ++ d)
                {
                    test(cubicEqn(a, b, c, d), tol);
                    ++nTest_;
                }
            }
        }
    }
    Info<< nCoeff*nCoeff*nCoeff*nCoeff
        << " label cubic equations were tested." << endl;


    if (nFail_)
    {
        Info<< nl << "        #### "
            << "Failed in " << nFail_ << " tests "
            << "out of total " << nTest_ << " tests "
            << "####\n" << endl;
        return 1;
    }

    Info<< nl << "        #### Passed all " << nTest_ <<" tests ####\n" << endl;
    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
