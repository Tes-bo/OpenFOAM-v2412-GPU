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
    This file is derivative work of OpenFOAM.

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
    Test-EigenMatrix

Description
    Tests for \c EigenMatrix constructors, and member functions
    using \c floatScalar, and \c doubleScalar base types.

    Cross-checks were obtained from 'NumPy 1.15.1' if no theoretical
    cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

\*---------------------------------------------------------------------------*/

#include "scalarMatrices.H"
#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "complex.H"
#include "IOmanip.H"
#include "EigenMatrix.H"
#include "TestTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create each constructor of EigenMatrix<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct from a SquareMatrix<Type>" << nl;
        const SquareMatrix<Type> A(5, Zero);
        const EigenMatrix<Type> EM(A);
    }
    {
        Info<< "# Construct from a SquareMatrix<Type> and symmetry flag" << nl;
        const SquareMatrix<Type> A(5, Zero);
        const EigenMatrix<Type> EM(A, true);
    }
}


// Execute each member function of EigenMatrix<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SquareMatrix<Type> A(3, Zero);
    assignMatrix
    (
        A,
        {
            Type(1), Type(2), Type(3),
            Type(4), Type(5), Type(6),
            Type(7), Type(8), Type(9)
        }
    );

    Info<< "# Operand: " << nl
        << "  SquareMatrix = " << A << endl;


    // Since eigenvalues are unique and eigenvectors are not unique,
    // the bitwise comparisons are limited to eigenvalue computations.
    // Here, only the execution of the functions is tested, rather than
    // the verification of the eigendecomposition through theoretical relations.
    {
        const EigenMatrix<Type> EM(A);
        const DiagonalMatrix<Type> EValsRe = EM.EValsRe();
        const DiagonalMatrix<Type>& EValsIm = EM.EValsIm();
        const SquareMatrix<Type>& EVecs = EM.EVecs();

        cmp
        (
            "  Return real eigenvalues or real part of complex eigenvalues = ",
            EValsRe,
            List<Type>
            ({
                Type(16.116844),
                Type(-1.116844),
                Type(0)
            }),
            getTol(Type(0)),
            1e-6
        );
        cmp
        (
            "  Return zero-matrix for real eigenvalues "
            "or imaginary part of complex eigenvalues = ",
            EValsIm,
            List<Type>(3, Zero)
        );
        Info<< "  Return eigenvectors matrix = " << EVecs << endl;
    }
}


// Test the relation: "sum(eigenvalues) = trace(A)"
// w.wiki/4zs (Retrieved: 16-06-19) # Item-1
template<class Type>
void test_eigenvalues_sum
(
    const SquareMatrix<Type>& A,
    const DiagonalMatrix<Type>& EValsRe
)
{
    const Type trace = A.trace();
    // Imaginary part of complex conjugates cancel each other
    const Type EValsSum = sum(EValsRe);
    Info<< "  # A.mRows = " << A.m() << nl;
    cmp
    (
        "  # sum(eigenvalues) = trace(A) = ",
        EValsSum,
        trace,
        getTol(Type(0))
    );
}


// Test the relation: "prod(eigenvalues) = det(A)"
// w.wiki/4zs (Retrieved: 16-06-19) # Item-2
// Note that the determinant computation may fail
// which is not a suggestion that eigendecomposition fails
template<class Type>
void test_eigenvalues_prod
(
    const SquareMatrix<Type>& A,
    const DiagonalMatrix<Type>& EValsRe,
    const DiagonalMatrix<Type>& EValsIm
)
{
    const Type determinant = mag(det(A));
    Type EValsProd = Type(1);

    if (EValsIm.empty())
    {
        for (label i = 0; i < EValsRe.size(); ++i)
        {
            EValsProd *= Foam::sqrt(sqr(EValsRe[i]));
        }
    }
    else
    {
        for (label i = 0; i < EValsRe.size(); ++i)
        {
            EValsProd *= Foam::sqrt(sqr(EValsRe[i]) + sqr(EValsIm[i]));
        }
    }
    cmp
    (
        "  # prod(eigenvalues) = det(A) = ",
        EValsProd,
        determinant,
        getTol(Type(0))
    );
}


// Test eigenvalues in eigendecomposition relations
// Relations: (Beauregard & Fraleigh (1973), ISBN 0-395-14017-X, p. 307)
template<class Type>
void test_eigenvalues(Type)
{
    Random rndGen(1234);
    const label numberOfTests = 20;

    // Non-symmetric
    for (label i = 0; i < numberOfTests; ++i)
    {
        const label mRows = rndGen.position(100, 200);
        const labelPair m(mRows, mRows);

        const SquareMatrix<Type> A
        (
            makeRandomMatrix<SquareMatrix<Type>>(m, rndGen)
        );

        const EigenMatrix<Type> EM(A);
        const DiagonalMatrix<Type>& EValsRe = EM.EValsRe();

        test_eigenvalues_sum(A, EValsRe);

        // LUDecompose does not work with floatScalar at the time of writing,
        // hence det function. Once LUDecompose is refactored, comment out below
        // const DiagonalMatrix<Type>& EValsIm = EM.EValsIm();
        // test_eigenvalues_prod(A, EValsRe, EValsIm);
    }

    // Symmetric
    for (label i = 0; i < numberOfTests; ++i)
    {
        const label mRows = rndGen.position(100, 200);
        const labelPair m(mRows, mRows);

        SquareMatrix<Type> A
        (
            makeRandomMatrix<SquareMatrix<Type>>(m, rndGen)
        );
        // Symmetrise with noise
        for (label n = 0; n < A.n() - 1; ++n)
        {
            for (label m = A.m() - 1; m > n; --m)
            {
                A(n, m) = A(m, n) + SMALL;
            }
        }

        const bool symmetric = true;
        const EigenMatrix<Type> EM(A, symmetric);
        const DiagonalMatrix<Type>& EValsRe = EM.EValsRe();

        test_eigenvalues_sum(A, EValsRe);
    }
}


// Test the relation: "(A & EVec - EVal*EVec) = 0"
template<class Type>
void test_characteristic_eq
(
    const SquareMatrix<Type>& Aorig,
    const DiagonalMatrix<Type>& EValsRe,
    const DiagonalMatrix<Type>& EValsIm,
    const SquareMatrix<complex>& EVecs
)
{
    SquareMatrix<complex> A(Aorig.m());
    auto convertToComplex = [&](const scalar& val) { return complex(val); };
    std::transform
    (
        Aorig.cbegin(),
        Aorig.cend(),
        A.begin(),
        convertToComplex
    );

    for (label i = 0; i < A.m(); ++i)
    {
        const RectangularMatrix<complex>& EVec(EVecs.subColumn(i));
        const complex EVal(EValsRe[i], EValsIm[i]);
        const RectangularMatrix<complex> leftSide(A*EVec);
        const RectangularMatrix<complex> rightSide(EVal*EVec);

        cmp
        (
            "  # (A & EVec - EVal*EVec) = 0:",
            flt(leftSide),
            flt(rightSide),
            getTol(Type(0))
        );
    }
}


// Test eigenvectors in eigendecomposition relations
template<class Type>
void test_eigenvectors(Type)
{
    Random rndGen(1234);
    const label numberOfTests = 20;

    // Non-symmetric
    for (label i = 0; i < numberOfTests; ++i)
    {
        const label mRows = rndGen.position(100, 200);
        const labelPair m(mRows, mRows);

        const SquareMatrix<Type> A
        (
            makeRandomMatrix<SquareMatrix<Type>>(m, rndGen)
        );

        const EigenMatrix<Type> EM(A);
        const DiagonalMatrix<Type>& EValsRe = EM.EValsRe();
        const DiagonalMatrix<Type>& EValsIm = EM.EValsIm();
        const SquareMatrix<complex> EVecs(EM.complexEVecs());

        test_characteristic_eq(A, EValsRe, EValsIm, EVecs);
    }

    // Symmetric
    for (label i = 0; i < numberOfTests; ++i)
    {
        const label mRows = rndGen.position(100, 200);
        const labelPair m(mRows, mRows);

        SquareMatrix<Type> A
        (
            makeRandomMatrix<SquareMatrix<Type>>(m, rndGen)
        );
        // Symmetrise with noise
        for (label n = 0; n < A.n() - 1; ++n)
        {
            for (label m = A.m() - 1; m > n; --m)
            {
                A(n, m) = A(m, n) + SMALL;
            }
        }

        const bool symmetric = true;
        const EigenMatrix<Type> EM(A, symmetric);
        const DiagonalMatrix<Type>& EValsRe = EM.EValsRe();
        const DiagonalMatrix<Type>& EValsIm = EM.EValsIm();
        const SquareMatrix<complex> EVecs(EM.complexEVecs());

        test_characteristic_eq(A, EValsRe, EValsIm, EVecs);
    }
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

    Info<< nl << "    ## Test eigenvalues: "<< typeID[I] <<" ##" << nl;
    test_eigenvalues(std::get<I>(types));

    Info<< nl << "    ## Test eigenvectors: "<< typeID[I] <<" ##" << nl;
    test_eigenvectors(std::get<I>(types));

    run_tests<I + 1, Tp...>(types, typeID);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    Info<< setprecision(15);

    const std::tuple<floatScalar, doubleScalar> types
    (
        std::make_tuple(Zero, Zero)
    );

    const List<word> typeID
    ({
        "SquareMatrix<floatScalar>",
        "SquareMatrix<doubleScalar>"
    });

    run_tests(types, typeID);


    Info<< nl << "    ## Test corner cases ##" << endl;
    {
        Info<< nl  << "  ## Rosser et al. (1951) matrix: ##" << nl;
        // Rosser, J. B., Lanczos, C., Hestenes, M. R., & Karush, W. (1951).
        // Separation of close eigenvalues of a real symmetric matrix.
        // Jour. Research of the National Bureau of Standards, 47(4), 291-297.
        // DOI:10.6028/jres.047.037
        //
        // 8x8 symmetric square matrix consisting of close real eigenvalues
        // ibid, p. 294
        // {
        //      1020.04901843, 1020.000, 1019.90195136, 1000.000,
        //      1000.000, 0.09804864072, 0.000, -1020.0490
        // }
        // Note that prod(eigenvalues) != determinant(A) for this matrix
        // via the LAPACK routine z/dgetrf

        SquareMatrix<doubleScalar> A(8, Zero);
        assignMatrix
        (
            A,
            {
                611, 196, -192, 407, -8, -52, -49, 29,
                196, 899, 113, -192, -71, -43, -8, -44,
                -192, 113, 899, 196, 61, 49, 8, 52,
                407, -192, 196, 611, 8, 44, 59, -23,
                -8, -71, 61, 8, 411, -599, 208, 208,
                -52, -43, 49, 44, -599, 411, 208, 208,
                -49, -8, 8, 59, 208, 208, 99, -911,
                29, -44, 52, -23, 208, 208, -911, 99
            }
        );
        const EigenMatrix<doubleScalar> EM(A);
        const DiagonalMatrix<doubleScalar>& EValsRe = EM.EValsRe();
        const DiagonalMatrix<doubleScalar>& EValsIm = EM.EValsIm();

        test_eigenvalues_sum(A, EValsRe);

        cmp
        (
            "  # Rosser et al. (1951) case, EValsRe = ",
            EValsRe,
            List<doubleScalar> // theoretical EValsRe
            ({
                -1020.0490, 0.000, 0.09804864072, 1000.000,
                1000.000, 1019.90195136, 1020.000, 1020.04901843
            }),
            1e-3
        );
        cmp
        (
            "  # Rosser et al. (1951) case, EValsIm = ",
            EValsIm,
            List<doubleScalar>(8, Zero)
        );
    }
    {
        Info<< nl  << "  ## Test eigenvector unpacking: ##" << nl;

        SquareMatrix<doubleScalar> A(3, Zero);
        assignMatrix
        (
            A,
            {
                1, 2, 3,
                -4, -5, 6,
                7, -8, 9
            }
        );
        const EigenMatrix<doubleScalar> EM(A);
        const SquareMatrix<complex> complexEVecs(EM.complexEVecs());

        SquareMatrix<complex> B(3, Zero);
        assignMatrix
        (
            B,
            {
                complex(-0.373220280),
                complex(0.417439996, 0.642691344),
                complex(0.417439996, -0.642691344),
                complex(-0.263919251),
                complex(-1.165275867, 0.685068715),
                complex(-1.165275867, -0.685068715),
                complex(-0.889411744),
                complex(-0.89990601, -0.3672785281),
                complex(-0.89990601, 0.3672785281),
            }
        );

        cmp
        (
            "  # ",
            flt(complexEVecs),
            flt(B)
        );
    }
    {
        Info<< nl << "  ## Test matrices with small values: ##" << nl;

        const List<doubleScalar> epsilons
        ({
            0, SMALL, Foam::sqrt(SMALL), sqr(SMALL), Foam::cbrt(SMALL),
            -SMALL, -Foam::sqrt(SMALL), -sqr(SMALL), -Foam::cbrt(SMALL)
        });

        Random rndGen(1234);
        const label numberOfTests = 20;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(100, 200);

            for (const auto& eps : epsilons)
            {
                const SquareMatrix<doubleScalar> A(mRows, eps);

                const EigenMatrix<doubleScalar> EM(A);
                const DiagonalMatrix<doubleScalar>& EValsRe = EM.EValsRe();
                const DiagonalMatrix<doubleScalar>& EValsIm = EM.EValsIm();
                const SquareMatrix<complex> EVecs(EM.complexEVecs());

                test_eigenvalues_sum(A, EValsRe);
                test_characteristic_eq(A, EValsRe, EValsIm, EVecs);
            }
        }
    }
    {
        Info<< nl << "  ## Test matrices with repeating eigenvalues: ##" << nl;
        SquareMatrix<doubleScalar> A(3, Zero);
        assignMatrix
        (
            A,
            {
                0, 1, 1,
                1, 0, 1,
                1, 1, 0
            }
        );
        const EigenMatrix<doubleScalar> EM(A);
        const DiagonalMatrix<doubleScalar>& EValsRe = EM.EValsRe();
        const DiagonalMatrix<doubleScalar>& EValsIm = EM.EValsIm();
        const SquareMatrix<complex> EVecs(EM.complexEVecs());

        test_eigenvalues_sum(A, EValsRe);
        test_characteristic_eq(A, EValsRe, EValsIm, EVecs);
    }


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
