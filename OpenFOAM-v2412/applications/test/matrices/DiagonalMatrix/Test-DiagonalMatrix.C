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
    Test-DiagonalMatrix

Description
    Tests for \c DiagonalMatrix constructors, member functions and global
    functions using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' if no theoretical
    cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

\*---------------------------------------------------------------------------*/

#include "DiagonalMatrix.H"
#include "RectangularMatrix.H"
#include "scalar.H"
#include "complex.H"
#include "TestTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create each constructor of DiagonalMatrix<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct empty from size:" << nl;
        const DiagonalMatrix<Type> A(5);
        Info<< A << endl;
    }
    {
        Info<< "# Construct from size and initialise all elems to zero:" << nl;
        const DiagonalMatrix<Type> A(5, Zero);
        Info<< A << endl;
    }
    {
        Info<< "# Construct from size and initialise all elems to value" << nl;
        const DiagonalMatrix<Type> A(5, Type(8));
        Info<< A << endl;
    }
    {
        Info<< "# Construct from the diagonal of a Matrix" << nl;
        const RectangularMatrix<Type> M(3, 5, Zero);
        const DiagonalMatrix<Type> A(M);
        Info<< A << endl;
    }
}


// Execute each member function of DiagonalMatrix<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    DiagonalMatrix<Type> A(3, Zero);
    assignMatrix(A, {Type(1), Type(2), Type(-3)});

    Info<< "# Operand: " << nl
        << "  DiagonalMatrix = " << A << endl;


    {
        Info<< "# Return the matrix inverse into itself:" << nl;
        A.invert();
        cmp
        (
            "  DiagonalMatrix<Type>.invert() = ",
            A,
            List<Type>({Type(1), Type(0.5), Type(-0.333333)}),
            1e-6
        );
    }
    {
        Info<< "# Sort:" << nl;

        DiagonalMatrix<Type> B(5, Zero);
        assignMatrix(B, {Type(1), Type(2), Type(-3), Type(5), Type(1.01)});

        auto descend = [&](Type a, Type b){ return mag(a) > mag(b); };
        const List<label> sortPermutation(B.sortPermutation(descend));
        cmp
        (
            "  Return a sort permutation labelList according to "
            "a given comparison on the diagonal entries",
            sortPermutation,
            List<label>({3, 2, 1, 4, 0})
        );

        DiagonalMatrix<Type> sortedB0(5, Zero);
        assignMatrix
        (
            sortedB0,
            {
                Type(5),
                        Type(-3),
                                 Type(2),
                                         Type(1.01),
                                                    Type(1)
            }
        );
        const DiagonalMatrix<Type> sortedB1
        (
            applyPermutation(B, sortPermutation)
        );
        cmp
        (
            "  Return Matrix column-reordered according to "
            "a given permutation labelList",
            sortedB0,
            sortedB1
        );

        DiagonalMatrix<Type> cpB(B);
        cpB.applyPermutation(sortPermutation);
        cmp
        (
            "  Column-reorder this Matrix according to "
            "a given permutation labelList",
            sortedB0,
            cpB
        );
    }

}


// Execute each global function of DiagonalMatrix<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    DiagonalMatrix<Type> A(3, Zero);
    assignMatrix(A, {Type(1), Type(2), Type(-3)});

    Info<< "# Operand: " << nl
        << "  DiagonalMatrix = " << A << nl << endl;


    cmp
    (
        "  Inverse = ",
        inv(A),
        List<Type>({Type(1), Type(0.5), Type(-0.333333)}),
        1e-6
    );
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
        "DiagonalMatrix<floatScalar>",
        "DiagonalMatrix<doubleScalar>",
        "DiagonalMatrix<complex>"
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
