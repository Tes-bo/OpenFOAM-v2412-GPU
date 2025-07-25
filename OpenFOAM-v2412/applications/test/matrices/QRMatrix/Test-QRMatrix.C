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
    Test-QRMatrix

Description
    Tests for \c QRMatrix constructors, and member functions
    using \c doubleScalar and \c complex base types.

\*---------------------------------------------------------------------------*/

#include "MatrixTools.H"
#include "QRMatrix.H"
#include "complex.H"
#include "IOmanip.H"
#include "TestTools.H"

#define PRINTALL false

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Verify if Q is a unitary matrix for rectangular matrices
template<class Type>
void verify_Q(const RectangularMatrix<Type>& Q)
{
    // non-square unitary matrices are not possible - do nothing
}


// Verify if Q is a unitary matrix for square matrices
template<class Type>
void verify_Q
(
    const SquareMatrix<Type>& Q
)
{
    // mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:18-06-19)
    const SquareMatrix<Type> QTQ(Q&Q);
    const SquareMatrix<Type> QQT(Q^Q);
    SquareMatrix<Type> IMatrix(Q.sizes(), Zero);
    for (label i = 0; i < IMatrix.n(); ++i)
    {
        IMatrix(i,i) = pTraits<Type>::one;
    }
    cmp("# Q.T()*Q = I: ", flt(QTQ), flt(IMatrix), getTol(Type(0)));
    cmp("# Q*Q.T() = I: ", flt(QQT), flt(IMatrix), getTol(Type(0)));
    cmp("# QTQ = QQT: ", flt(QTQ), flt(QQT), getTol(Type(0)));
}


// Verify if R is an upper-triangular matrix
template<class MatrixType>
void verify_R
(
    const MatrixType& R
)
{
    DynamicList<scalar> ltriR;
    label k = 0;
    // mathworld.wolfram.com/UpperTriangularMatrix.html (Retrieved:18-06-19)
    for (label i = 0; i < R.m(); ++i)
    {
        for (label j = 0; j < R.n(); ++j)
        {
            if (j < i)
            {
                ltriR.append(mag(R(i, j)));
                ++k;
            }
        }
    }
    const List<scalar> ltri(k, Zero);
    cmp("# R(i, j) = 0 for i > j: ", ltri, ltriR);
}


// Verify if R is an upper-triangular matrix with diag elems non-increasing
template<class MatrixType>
void verify_R_pivoting
(
    const MatrixType& R
)
{
    // w.wiki/54m (Retrieved:20-06-19)
    const List<scalar> diag0(mag(R.diag()));
    List<scalar> diag1(diag0);
    Foam::sort(diag1, std::greater<scalar>());
    cmp
    (
        "# Diag elems of R non-increasing: |R11| >= |R22| >= ... >= |Rnn|: ",
        diag0,
        diag1
    );
}


// Verify if the given matrix can be reconstructed by Q and R matrices
template<class MatrixType>
void verify_QR
(
    const MatrixType& A,
    const MatrixType& Q,
    const MatrixType& R
)
{
    // mathworld.wolfram.com/QRDecomposition.html (Retrieved:18-06-19)
    const MatrixType AReconstruct(Q*R);
    cmp("# Q*R = A: ", flt(A), flt(AReconstruct));
}


// Verify if the given matrix can be reconstructed by Q, R and P matrices
template<class MatrixType>
void verify_QR_pivoting
(
    const MatrixType& A,
    const MatrixType& Q,
    const MatrixType& R,
    const SquareMatrix<typename MatrixType::cmptType>& P
)
{
    // w.wiki/54m (Retrieved:20-06-19)
    const MatrixType AReconstruct(Q*R);
    const MatrixType AP(A*P);
    cmp("# Q*R = A*P: ", flt(AReconstruct), flt(AP));
}


// Verify the given QR decomposition
template<class MatrixType>
void verify
(
    const QRMatrix<MatrixType>& QRM,
    const MatrixType& A
)
{
    const MatrixType& Q = QRM.Q();
    const MatrixType& R = QRM.R();

    #if PRINTALL
    InfoNumPyFormat(A, "A");
    InfoNumPyFormat(Q, "Q");
    InfoNumPyFormat(R, "R");
    #endif

    verify_Q(Q);
    verify_R(R);
    verify_QR(A, Q, R);
    Info<< endl;
}


// Verify the given QR decomposition with column pivoting
template<class MatrixType>
void verify_pivoting
(
    const QRMatrix<MatrixType>& QRM,
    const MatrixType& A
)
{
    const MatrixType& Q = QRM.Q();
    const MatrixType& R = QRM.R();
    const SquareMatrix<typename MatrixType::cmptType>& P = QRM.P();

    #if PRINTALL
    InfoNumPyFormat(A, "A");
    InfoNumPyFormat(Q, "Q");
    InfoNumPyFormat(R, "R");
    InfoNumPyFormat(P, "P");
    #endif

    verify_Q(Q);
    verify_R(R);
    verify_R_pivoting(R);
    verify_QR_pivoting(A, Q, R, P);
    Info<< endl;
}


// Verify various QR decompositions
template<class MatrixType>
void verify_decomposition
(
    const MatrixType& A
)
{
    Info<< "## Verify various QR decompositions" << nl << endl;
    {
        Info<< "# modes::FULL" << nl;
        QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::FULL
        );
        verify(QRM, A);
    }
    {
        Info<< "# modes::ECONOMY" << nl;
        QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::ECONOMY
        );
        verify(QRM, A);
    }
    {
        Info<< "# modes::FULL, outputs::BOTH_QR, column-pivoting" << nl;
        QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::FULL,
            QRMatrix<MatrixType>::outputs::BOTH_QR,
            QRMatrix<MatrixType>::pivoting::TRUE
        );
        verify_pivoting(QRM, A);
    }
    {
        Info<< "# modes::ECONOMY, outputs::BOTH_QR, column-pivoting" << nl;
        QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::ECONOMY,
            QRMatrix<MatrixType>::outputs::BOTH_QR,
            QRMatrix<MatrixType>::pivoting::TRUE
        );
        verify_pivoting(QRM, A);
    }
}


// Verify the Moore-Penrose pseudo-inverse function
template<class MatrixType>
void verify_pinv
(
    const MatrixType& A
)
{
    const scalar absTol = 1e-6;

    Info<< "## Verify Moore-Penrose pseudo-inverse: pinv()" << nl << endl;
    const MatrixType APlus(MatrixTools::pinv(A));

    #if PRINTALL
    InfoNumPyFormat(A, "A");
    InfoNumPyFormat(APlus, "A^+");
    #endif

    {
        const MatrixType APlusPlus(MatrixTools::pinv(APlus));
        cmp("# (A^+)^+ = A: ", flt(APlusPlus), flt(A), absTol);
    }

    // The Moore-Penrose conditions
    // w.wiki/5RL (Retrieved: 29-06-19)
    {
        const MatrixType AAPlusA((A*APlus)*A);
        cmp("# A*(A^+)*A = A: ", flt(AAPlusA), flt(A), absTol);
    }

    {
        const MatrixType APlusAAPlus((APlus*A)*APlus);
        cmp("# (A^+)*A*(A^+) = A: ", flt(APlusAAPlus), flt(APlus), absTol);
    }

    {
        const MatrixType AAPlus(A*APlus);
        const MatrixType AAPlusT(AAPlus.T());
        cmp("# (A*(A^+)).T() = A*(A^+): ", flt(AAPlusT), flt(AAPlus), absTol);
    }

    {
        const MatrixType APlusA(APlus*A);
        const MatrixType APlusAT(APlusA.T());
        cmp("# ((A^+)*A = ((A^+)*A).T(): ", flt(APlusA), flt(APlusAT), absTol);
    }
}


// Verify the direct solution of a linear system by QRMatrix
void verify_solve
(
    const SquareMatrix<complex>& A
)
{
    // do nothing
}


// Verify the direct solution of a linear system by QRMatrix
void verify_solve
(
    const SquareMatrix<scalar>& A
)
{
    Info<< "## Verify direct solution of A*x = b with A = Q*R:" << nl << endl;

    typedef SquareMatrix<scalar> SMatrix;
    const scalar absTol = 1e-6;

    // Direct solution of a linear system by QR decomposition
    QRMatrix<SMatrix> QRM
    (
        A,
        QRMatrix<SMatrix>::modes::FULL,
        QRMatrix<SMatrix>::outputs::BOTH_QR
    );

    const scalarField residual(A.m(), 0);
    const scalarField source(A.m(), 1);

    const scalarField x(QRM.solve(source));
    const scalarField AxMinusSource(A*x - source);
    const scalarField xMinusQRinvSource(x - QRM.inv()*source);
    const SMatrix QRinvA(QRM.inv()*A);
    const SMatrix identity(A.m(), I);

    cmp("# A*x - b = residual: ", AxMinusSource, residual, absTol);
    cmp("# x - inv(QR)*b = residual: ", xMinusQRinvSource, residual, absTol);
    cmp("# inv(QR)*A = I: ", flt(QRinvA), flt(identity), absTol);
}


// Verify the parallel tall-skinny QR decomposition
template<class Type>
void verify_tsqr
(
    Random& rndGen
)
{
    Info<< "## Verify parallel direct tall-skinny QR decomposition:"
        << nl << endl;

    typedef RectangularMatrix<Type> RMatrix;

    // Size of the full matrix and its partitions
    const label nColsSum = rndGen.position<label>(1, 100);
    const label qParts = rndGen.position<label>(10, 30);
    List<label> mRowsList(qParts, Zero);

    label mRowsSum = 0;
    for (label i = 0; i < qParts; ++i)
    {
        const label mRows = rndGen.position<label>(nColsSum, 10*nColsSum);
        mRowsList[i] = mRows;
        mRowsSum += mRows;
    }

    #if PRINTALL
    Info<< "mRowsSum = " << mRowsSum << nl;
    Info<< "nColsSum = " << nColsSum << nl;
    Info<< "mRowsList = " << mRowsList << nl;
    #endif

    if (mRowsSum < nColsSum)
    {
        FatalErrorInFunction
            << "Tall skinny QR decomposition cannot be executed for matrices "
            << "with number of columns (nCols) = " << nColsSum << " > "
            << "than number of rows (mRows) = " << mRowsSum << "."
            << "Yet nCols > mRows is allowed for a submatrix."
            << exit(FatalError);
    }

    // Perform the QR decomposition for the full matrix
    const RMatrix A(makeRandomMatrix<RMatrix>({mRowsSum, nColsSum}, rndGen));
    QRMatrix<RMatrix> QRM
    (
        A,
        QRMatrix<RMatrix>::modes::ECONOMY,
        QRMatrix<RMatrix>::outputs::ONLY_R,
        QRMatrix<RMatrix>::pivoting::FALSE
    );
    const RMatrix& R = QRM.R();

    RMatrix RMag(R);
    for (auto& val : RMag)
    {
        val = mag(val);
    }

    // Partition the full matrix,
    // and perform the QR decomposition on each partition
    RMatrix subRs({qParts*nColsSum, nColsSum}, Zero);
    label mRowsIndex = 0;
    for (label i = 0; i < qParts; ++i)
    {
        RMatrix subA({mRowsList[i], nColsSum}, Zero);
        subA = A.subMatrix(mRowsIndex, 0, mRowsList[i]);
        mRowsIndex += mRowsList[i];

        QRMatrix<RMatrix> subQRM
        (
            subA,
            QRMatrix<RMatrix>::modes::ECONOMY,
            QRMatrix<RMatrix>::outputs::ONLY_R,
            QRMatrix<RMatrix>::pivoting::FALSE
        );
        const RMatrix& subR = subQRM.R();

        // Append subR
        subRs.subMatrix(i*nColsSum, 0, nColsSum) = subR;
    }

    // Perform the tall-skinny QR decomposition
    QRMatrix<RMatrix> TSQR
    (
        subRs,
        QRMatrix<RMatrix>::modes::ECONOMY,
        QRMatrix<RMatrix>::outputs::ONLY_R,
        QRMatrix<RMatrix>::pivoting::FALSE
    );
    const RMatrix& TSQRR = TSQR.R();

    RMatrix TSQRRMag(TSQRR);
    for (auto& val : TSQRRMag)
    {
        val = mag(val);
    }

    // Compare magnitude of R, since R is not unique up to the sign
    cmp("# TSQR: ", flt(RMag), flt(TSQRRMag));

    #if PRINTALL
    InfoNumPyFormat(R, "R");
    InfoNumPyFormat(TSQRR, "TSQRR");
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class MatrixType>
void test_constructors(MatrixType)
{
    {
        Info<< "# Construct from a Matrix and modes" << nl;
        MatrixType A({5,5}, Zero);
        const QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::FULL
        );
    }
    {
        Info<< "# Construct from a const Matrix and modes" << nl;
        const MatrixType A({5,5}, Zero);
        const QRMatrix<MatrixType> QRM
        (
            A,
            QRMatrix<MatrixType>::modes::FULL
        );
    }
}


template<class MatrixType>
void test_decomposition(MatrixType)
{
    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> SMatrix;
    typedef RectangularMatrix<cmptType> RMatrix;

    Random rndGen(1234);
    const label szTests = 20;
    const label min = 10;
    const label max = 50;

    {
        for (label i = 0; i < szTests; ++i)
        {
            const label m = rndGen.position(min, max);
            const label n = rndGen.position(min, max);

            const RMatrix A
            (
                makeRandomMatrix<RMatrix>({m,n}, rndGen)
            );

            verify_decomposition(A);
            verify_pinv(A);
        }
    }

    {
        for (label i = 0; i < szTests; ++i)
        {
            const label m = rndGen.position(min, max);

            const SMatrix A
            (
                makeRandomMatrix<SMatrix>({m,m}, rndGen)
            );

            verify_decomposition(A);
            verify_pinv(A);
            verify_solve(A);
        }
    }

    for (label i = 0; i < szTests; ++i)
    {
        verify_tsqr<cmptType>(rndGen);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    Info<< nl << "    ## Test decomposition: "<< typeID[I] <<" ##" << nl;
    test_decomposition(std::get<I>(types));

    run_tests<I + 1, Tp...>(types, typeID);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    Info<< setprecision(15);

    const std::tuple
    <
        RectangularMatrix<doubleScalar>,
        RectangularMatrix<complex>,
        SquareMatrix<doubleScalar>,
        SquareMatrix<complex>
    > types
    (
        std::make_tuple(Zero, Zero, Zero, Zero)
    );

    const List<word> typeID
    ({
        "RectangularMatrix<doubleScalar>",
        "RectangularMatrix<complex>",
        "SquareMatrix<doubleScalar>",
        "SquareMatrix<complex>"
    });

    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();

    run_tests(types, typeID);

    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    Info<< "Time elapsed = "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count()
        << "[µs]" << endl;

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
