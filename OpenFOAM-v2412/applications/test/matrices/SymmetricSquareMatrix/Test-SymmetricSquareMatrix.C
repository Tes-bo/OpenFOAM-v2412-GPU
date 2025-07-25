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
    Test-SymmetricSquareMatrix

Description
    Tests for \c SymmetricSquareMatrix constructors, member functions, member
    operators, global functions, global operators, and friend functions using
    \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' if no theoretical
    cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

Note
    Pending tests were tagged as "## Pending ##".

\*---------------------------------------------------------------------------*/

#include "scalarMatrices.H"
#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "SymmetricSquareMatrix.H"
#include "scalar.H"
#include "complex.H"
#include "IOmanip.H"
#include "Random.H"
#include "TestTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create each constructor of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct for given size (rows == cols):" << nl;
        const SymmetricSquareMatrix<Type> A(2);
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing all elements to zero:" << nl;
        const SymmetricSquareMatrix<Type> A(2, Zero);
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing all elements to a value:" << nl;
        const SymmetricSquareMatrix<Type> A(2, Type(3));
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing to the identity matrix:" << nl;
        const SymmetricSquareMatrix<Type> A(2, I);
        Info<< A << endl;
    }

    //  ## Pending ##
    //  Construct from Istream
    //  Clone
    //  #############
}


// Execute each member function of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    //  ## Pending ##
    //  #############
}


// Execute each member operators of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_member_opers(Type)
{
    //  ## Pending ##
    //  #############
}


// Execute each friend function of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_friend_funcs(Type)
{
    //  ## Pending ##
    //  #############
}


// Execute each global function of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    //  ## Pending ##
    //  #############
}


// Execute each global operators of SymmetricSquareMatrix<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    //  ## Pending ##
    //  #############
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

    Info<< nl << "    ## Test member opers: "<< typeID[I] <<" ##" << nl;
    test_member_opers(std::get<I>(types));

    Info<< nl << "    ## Test global functions: "<< typeID[I] << " ##" << nl;
    test_global_funcs(std::get<I>(types));

    Info<< nl << "    ## Test global operators: "<< typeID[I] << " ##" << nl;
    test_global_opers(std::get<I>(types));

    Info<< nl << "    ## Test friend funcs: "<< typeID[I] <<" ##" << nl;
    test_friend_funcs(std::get<I>(types));

    run_tests<I + 1, Tp...>(types, typeID);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    Info<< setprecision(15);

    const std::tuple<floatScalar, doubleScalar, complex> types
    (
        std::make_tuple(Zero, Zero, Zero)
    );

    const List<word> typeID
    ({
        "SymmetricSquareMatrix<floatScalar>",
        "SymmetricSquareMatrix<doubleScalar>",
        "SymmetricSquareMatrix<complex>"
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


// ************************************************************************* //
