/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-linearEqn

Description
    Tests for \c linearEqn constructors, and member functions.

\*---------------------------------------------------------------------------*/

#include <ctime>
#include <random>

#include "linearEqn.H"

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


// Create each constructor of linearEqn, and print output
void test_constructors()
{
    #if 0
    {
        Info<< "# Construct null:" << nl;
        const linearEqn T();
        Info<< T << endl;
    }
    #endif

    {
        Info<< "# Construct initialized to zero:" << nl;
        const linearEqn T(Zero);
        Info<< T << endl;
    }
    {
        Info<< "# Construct from components:" << nl;
        const linearEqn T(1, 2);
        Info<< T << endl;
    }
}


// Execute each member function of linearEqn, and print output
void test_member_funcs()
{
    const linearEqn T(1, -6);
    const scalar x = 9.7;

    Info<< "# Operands: " << nl
        << "  linearEqn = " << T << nl
        << "  x = " << x << endl;


    {
        Info<< "# Access:" << nl;
        cmp("  a = ",  1, T.a());
        cmp("  b = ",  -6, T.b());
    }
    {
        Info<< "# Evaluate:" << nl;
        cmp("  T.value(x) = ", T.value(x), 3.6999999999999993);
        cmp("  T.derivative(x) = ", T.derivative(x), 1);
        cmp
        (
            "  T.error(x) = ",
            T.error(x),
            1.5699999999999999e-15,
            1e-8,
            1e-8
        );
        cmp("  T.roots() = ", T.roots()[0], 6);
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


    // Sanity checks
    const scalar tol = 5;

    const label nTests = 1000000;
    for (label t = 0; t < nTests; ++t)
    {
        test
        (
            linearEqn
            (
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50)
            ),
            tol
        );
        ++nTest_;
    }
    Info << nTests << " scalar linear equations were tested." << endl;

    const label coeffMin = -9, coeffMax = 10, nCoeff = coeffMax - coeffMin;
    for (label a = coeffMin; a < coeffMax; ++a)
    {
        for (label b = coeffMin; b < coeffMax; ++b)
        {
            test(linearEqn(a, b), tol);
            ++nTest_;
        }
    }
    Info<< nCoeff*nCoeff << " label linear equations were tested." << endl;


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
