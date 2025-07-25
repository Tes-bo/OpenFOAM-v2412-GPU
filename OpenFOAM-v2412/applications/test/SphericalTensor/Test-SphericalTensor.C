/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Test-SphericalTensor

Description
    Tests for \c SphericalTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "Tensor.H"
#include "SymmTensor.H"
#include "SphericalTensor.H"
#include "DiagTensor.H"
#include "scalar.H"
#include "complex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Compare two floating point types, and print output.
// Do ++nFail_ if values of two objects are not equal within a given tolerance.
// The function is converted from PEP-485.
template<class Type>
typename std::enable_if<pTraits<Type>::rank == 0, void>::type
cmp
(
    const word& msg,
    const Type& x,
    const Type& y,
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


// Compare two containers elementwise, and print output.
// Do ++nFail_ if two components are not equal within a given tolerance.
// The function is converted from PEP-485
template<class Type>
typename std::enable_if<pTraits<Type>::rank != 0, void>::type
cmp
(
    const word& msg,
    const Type& x,
    const Type& y,
    const scalar relTol = 1e-8,
    const scalar absTol = 0
)
{
    Info<< msg << x << endl;

    unsigned nFail = 0;

    for (direction i = 0; i < pTraits<Type>::nComponents; ++i)
    {
        if (max(absTol, relTol*max(mag(x[i]), mag(y[i]))) < mag(x[i] - y[i]))
        {
            ++nFail;
        }
    }

    if (nFail)
    {
        Info<< nl
            << "        #### Fail in " << nFail << " comps ####" << nl << endl;
        ++nFail_;
    }
    ++nTest_;
}


// Create each constructor of SphericalTensor<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const SphericalTensor<Type> spT(Zero);
        Info<< spT << endl;
    }
    {
        Info<< "# Construct given VectorSpace of the same rank:" << nl;
        const VectorSpace<SphericalTensor<Type>, Type, 1> V(Zero);
        const SphericalTensor<Type> spT(V);
        Info<< spT << endl;
    }
    {
        Info<< "# Construct given the component:" << nl;
        const SphericalTensor<Type> spT(Type(1));
        Info<< spT << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const SphericalTensor<Type> spT(Zero);
        const SphericalTensor<Type> copyspT(spT);
        Info<< spT << tab << copyspT << endl;
    }
}


// Execute each member function of SphericalTensor<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SphericalTensor<Type> spT(Type(1));
    const SphericalTensor<Type> cspT(Type(-9));

    Info<< "# Operand: " << nl
        << "  SphericalTensor = " << spT << endl;


    {
        Info<< "# Component access:" << nl;

        SphericalTensor<Type> cpspT(spT.ii());
        cmp("  'SphericalTensor' access:", spT, cpspT);

        const SphericalTensor<Type> cpcspT(cspT.ii());
        cmp("  'const SphericalTensor' access:", cspT, cpcspT);
    }
    {
        Info<< "# SphericalTensor operations:" << nl;

        Info<< "  Transpose:" << nl;
        cmp("  'SphericalTensor'.T():", spT.T(), spT);
    }
}


// Execute each global function of SphericalTensor<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const SphericalTensor<Type> spT(Type(5));

    Info<< "# Operand: " << nl
        << "  SphericalTensor = " << spT << endl;


    cmp("  Trace = ", tr(spT), Type(15));
    cmp("  Spherical part = ", sph(spT), spT);
    cmp("  Determinant = ", det(spT), Type(124.99999999999994));
    cmp
    (
        "  Inverse = ",
        inv(spT),
        SphericalTensor<Type>(Type(0.2))
    );
    cmp("  Square of Frobenius norm = ", magSqr(spT), Type(75));
    cmp("  Max component = ", cmptMax(spT), Type(5));
    cmp("  Min component = ", cmptMax(spT), Type(5));
    cmp("  Sum of components = ", cmptSum(spT), Type(15));
    cmp("  Arithmetic average of components = ", cmptAv(spT), Type(5));
}


// Execute each global operator of SphericalTensor<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor<Type> T
    (
        Type(-1), Type(2), Type(-3),
        Type(4),  Type(5), Type(-6),
        Type(7),  Type(8), Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(-1), Type(2), Type(-3),
                  Type(5), Type(-6),
                           Type(-9)
    );
    const DiagTensor<Type> dT(Type(1), Type(5), Type(-9));
    const SphericalTensor<Type> spT(Type(-2));
    const Vector<Type> v(Type(3), Type(2), Type(1));
    const Type x(4);

    Info<< "# Operands:" << nl
        << "  Tensor = " << T << nl
        << "  SymmTensor = " << sT << nl
        << "  DiagTensor = " << dT << nl
        << "  SphericalTensor = " << spT << nl
        << "  Vector = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Division of Type by SpTensor = ",
        (x/spT),
        SphericalTensor<Type>(Type(-2))
    );
    cmp
    (
        "  Division of SpTensor by Type = ",
        (spT/x),
        SphericalTensor<Type>(Type(-0.5))
    );
    cmp
    (
        "  Inner-product of SpTensor-SpTensor = ",
        (spT & spT),
        SphericalTensor<Type>(Type(4))
    );
    cmp
    (
        "  Inner-product of SpTensor-Vector = ",
        (spT & v),
        Vector<Type>(Type(-6), Type(-4), Type(-2)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector-SpTensor = ",
        (v & spT),
        Vector<Type>(Type(-6), Type(-4), Type(-2)) // Row-vector
    );
    cmp("  D-inner-product of SpTensor-SpTensor = ", (spT && spT), Type(12));
}


// Do compile-time recursion over the given types
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
run_tests(const std::tuple<Tp...>& types, const List<word>& typeID){}


template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
run_tests(const std::tuple<Tp...>& types, const List<word>& typeID)
{
    Info<< nl << "    ## Test constructors: "<< typeID[I] <<" ##" << nl;
    test_constructors(std::get<I>(types));

    Info<< nl << "    ## Test member functions: "<< typeID[I] <<" ##" << nl;
    test_member_funcs(std::get<I>(types));

    Info<< nl << "    ## Test global functions: "<< typeID[I] << " ##" << nl;
    test_global_funcs(std::get<I>(types));

    Info<< nl << "    ## Test global operators: "<< typeID[I] <<" ##" << nl;
    test_global_opers(std::get<I>(types));

    run_tests<I + 1, Tp...>(types, typeID);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    const std::tuple<floatScalar, doubleScalar, complex> types
    (
        std::make_tuple(Zero, Zero, Zero)
    );

    const List<word> typeID
    ({
        "SphericalTensor<floatScalar>",
        "SphericalTensor<doubleScalar>",
        "SphericalTensor<complex>"
    });

    run_tests(types, typeID);


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
