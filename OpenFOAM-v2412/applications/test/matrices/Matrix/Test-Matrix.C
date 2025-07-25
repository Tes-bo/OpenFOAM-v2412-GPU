/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "MatrixTools.H"
#include "LUscalarMatrix.H"
#include "LLTMatrix.H"
#include "Random.H"
#include "SortList.H"
#include "Switch.H"

using namespace Foam;
using namespace Foam::MatrixTools;

#define RUNALL true
#define isEqual MatrixTools::equal

void horizontalLine()
{
    Info<< "+---------+---------+---------+---------+---------+" << nl;
}


// Copy values into matrix
template<class Form, class Type>
void assignMatrix
(
    Matrix<Form, Type>& mat,
    std::initializer_list<typename Matrix<Form, Type>::cmptType> list
)
{
    const label nargs = list.size();

    if (nargs != mat.size())
    {
        FatalErrorInFunction
            << "Mismatch in matrix dimension ("
            << mat.m() << ", "
            << mat.n() << ") and number of args (" << nargs << ')' << nl
            << exit(FatalError);
     }

    std::copy(list.begin(), list.end(), mat.begin());
}


template<class MatrixType>
MatrixType makeRandomMatrix
(
    const labelPair& dims,
    Random& rndGen
)
{
    MatrixType mat(dims);

    std::generate
    (
        mat.begin(),
        mat.end(),
        [&]{return rndGen.GaussNormal<typename MatrixType::cmptType>();}
    );

    return mat;
}


template<class MatrixType>
Ostream& operator<<(Ostream& os, const ConstMatrixBlock<MatrixType>& mat)
{
    return MatrixTools::printMatrix(os, mat);
}


template<class MatrixType>
Ostream& operator<<(Ostream& os, const MatrixBlock<MatrixType>& mat)
{
    return MatrixTools::printMatrix(os, mat);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    typedef SquareMatrix<scalar> SMatrix;
    typedef SquareMatrix<complex> SCMatrix;
    typedef RectangularMatrix<scalar> RMatrix;
    typedef RectangularMatrix<complex> RCMatrix;

    Random rndGen(1234);

    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Matrix Constructors & Assignments:" << nl << nl;

        Info<< "# SquareMatrix<scalar> examples:" << nl;
        {
            const label mRows = 2;
            const label nCols = 2;
            Matrix<SMatrix, scalar> MS0(mRows, nCols);
            Matrix<SMatrix, scalar> MS1(mRows, nCols, Zero);
            Matrix<SMatrix, scalar> MS2(mRows, nCols, 17.44);
            Matrix<SMatrix, scalar> MS3(MS2);

            SMatrix S0(mRows);
            SMatrix S1(mRows, Zero);
            SMatrix S2(mRows, I);
            SMatrix S3(mRows, 17.44);
            SMatrix S4(2, Zero);
            SMatrix S5(labelPair{mRows, nCols});
            SMatrix S6(mRows, nCols, Zero);
            SMatrix S7({mRows, nCols}, Zero);
            SMatrix S8({mRows, nCols}, 17.44);
            assignMatrix
            (
                S8,
                {
                    1, 2,    // 0th row
                    3, 4     // 1st row and so on
                }
            );

            S1 = Zero;
            S2 = 13.13;
            S3 = I;
        }

        Info<< "# RectangularMatrix<scalar> examples:" << nl;
        {
            const label mRows = 2;
            const label nCols = 3;
            Matrix<RMatrix, scalar> MR0(mRows, nCols);
            Matrix<RMatrix, scalar> MR1(mRows, nCols, Zero);
            Matrix<RMatrix, scalar> MR2(mRows, nCols, 17.44);

            RMatrix R0(mRows, nCols);
            RMatrix R1(labelPair{mRows, nCols});
            RMatrix R2(mRows, nCols, Zero);
            RMatrix R3({mRows, nCols}, Zero);
            RMatrix R4(mRows, nCols, -17.44);
            RMatrix R5({mRows, nCols}, -17.44);
            RMatrix R6(3, 4, Zero);
            RMatrix R7({3, 4}, Zero);
            assignMatrix
            (
                R7,
                {
                    1, 2, 3, 4,
                    5, 6, 7, 8.8,
                    9, 10, 11, 12
                }
            );
            R0 = Zero;
            R2 = -13.13;
        }

        Info<< "# SquareMatrix<complex> examples:" << nl;
        {
            const label mRows = 2;
            const label nCols = 2;
            Matrix<SCMatrix, complex> MS0(mRows, nCols);
            Matrix<SCMatrix, complex> MS1(mRows, nCols, Zero);
            Matrix<SCMatrix, complex> MS2(mRows, nCols, complex(17.44, 44.17));
            Matrix<SCMatrix, complex> MS3(MS2);

            SCMatrix S0(mRows);
            SCMatrix S1(mRows, Zero);
            SCMatrix S2(mRows, I);
            SCMatrix S3(mRows, complex(17.44,0));
            SCMatrix S4(2, Zero);
            SCMatrix S5(labelPair{mRows, nCols});
            SCMatrix S6(mRows, nCols, Zero);
            SCMatrix S7({mRows, nCols}, Zero);
            SCMatrix S8({mRows, nCols}, complex(17.44,0));
            assignMatrix
            (
                S8,
                {
                    complex(1,0), complex(2,2),
                    complex(4,2), complex(5,0)
                }
            );

            S1 = Zero;
            S2 = complex(13.13,0);
            S3 = I;
        }

        Info<< nl;
        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Access Functions:" << nl;

        SMatrix S(3, Zero);
        assignMatrix
        (
            S,
            {
                -3, 11, -4,
                2, 3, 10,
                2, 6, 1
            }
        );
        S(0,1) = 10;

        Info<< "# SquareMatrix<scalar> example:" << nl;
        printMatrix(Info, S) << nl;

        Info
            << "Number of rows =" << tab << S.m() << nl
            << "Number of columns = " << tab << S.n() << nl
            << "Number of elements = " << tab << S.size() << nl
            << "Number of rows/columns = " << tab << S.sizes() << nl
            << "Matrix is empty = " << tab << Switch(S.empty()) << nl
            << "Constant pointer = " << tab << name(S.cdata()) << nl
            << "Pointer = " << tab << name(S.data()) << nl
            << nl;

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Block Access Functions:" << nl;

        RMatrix R(4, 3, Zero);
        assignMatrix
        (
            R,
            {
                -3, 11, -4,
                2, 3, 10,
                2, 6, 1,
                1, 2, 3
            }
        );

        Info<< "# RectangularMatrix<scalar> example:" << nl;
        printMatrix(Info, R) << nl;

        // Indices begin at 0
        const label columnI = 1;
        const label rowI = 2;
        const label colStartI = 1;
        const label rowStartI = 1;
        const label szRows = 2;
        const label szCols = 2;

        Info
            << "column: " << R.subColumn(columnI) << nl
            << "sub-column: " << R.subColumn(columnI, rowStartI) << nl
            << "sub-sub-column: " << R.subColumn(columnI, rowStartI, szRows) << nl
            << "row: " << R.subRow(rowI) << nl
            << "sub-row: " << R.subRow(rowI, colStartI) << nl
            << "sub-sub-row: " << R.subRow(rowI, colStartI, szCols) << nl
            << "sub-block: " << R.subMatrix(rowStartI, colStartI) << nl
            << "sub-block with row size: " << R.subMatrix(rowStartI, colStartI, szRows) << nl
            << "sub-block with row|col size: " << R.subMatrix(rowStartI, colStartI, szRows, szCols) << nl;

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Edit Functions:" << nl;

        RMatrix R(2, 2, Zero);
        assignMatrix
        (
            R,
            {
                1, 2,
                3, 4
            }
        );

        Info<< "# Before Matrix.clear():" << nl << R << nl;
        R.clear();
        Info<< "# After Matrix.clear():" << nl << R << nl;

        Info<< "# Before Matrix.resize(m, n):" << nl << R;
        R.resize(3, 2);
        Info<< "# After Matrix.resize(m, n):" << nl << R << nl << nl;

        // Steal matrix contents
        Info<< "# Before Matrix.release():" << nl << R << nl;
        List<scalar> newOwner(R.release());
        Info<< "# After Matrix.release():" << nl << R
            << "& corresponding List content:" << nl << newOwner << nl << nl;

        R.resize(3, 2);
        assignMatrix
        (
            R,
            {
                1, 2,
                5, 6,
                9e-19, 1e-17
            }
        );

        Info<< "# Before Matrix.round():" << nl << R << nl;
        R.round();
        Info<< "# After Matrix.round():" << nl << R << nl << nl;

        Info<< "# Before Matrix.T() (Transpose):" << nl << R << nl;
        RMatrix RT = R.T();
        Info<< "# After Matrix.T():" << nl << RT << nl << nl;

        RCMatrix RC(3, 2, Zero);
        assignMatrix
        (
            RC,
            {
                complex(1, 0), complex(2,5),
                complex(4, 3), complex(1,-1),
                complex(1, 2), complex(2,9)
            }
        );

        Info<< "# Before Matrix.T() (Hermitian transpose):" << nl << RC << nl;
        RCMatrix RCH = RC.T();
        Info<< "# After Matrix.T():" << nl << RCH << nl << nl;

        Info<< "# Diagonal and trace:" << nl;
        RMatrix R1(5, 7, 3.4);
        {
            List<scalar> diagR = R1.diag();
            scalar traceR = R1.trace();
            Info<< "RectangularMatrix<scalar>:"
                << nl << R1 << nl
                << "Major diagonal of RectangularMatrix:"
                << nl << diagR << nl
                << "Trace of RectangularMatrix:"
                << nl << traceR << nl << nl;
        }

        Info<< "# List assignment to the main diagonal of Matrix:" << nl;
        {
            SMatrix S1(5, I);
            Info<< "Before List assignment:" << nl;
            printMatrix(Info, S1);
            List<scalar> lst(5, 2.0);
            S1.diag(lst);
            Info<< "After List assignment:" << nl;
            printMatrix(Info, S1);
        }

        Info<< "# Diagonal sort of Matrix:" << nl;
        {
            auto S1
            (
                makeRandomMatrix<SMatrix>
                (
                    {3, 3}, rndGen
                )
            );

            List<scalar> diag = S1.diag();
            SortList<scalar> incrDiag(diag);
            incrDiag.sort();

            Info<< "Diagonal list:" << diag << nl
                << "Ascending diagonal list: " << incrDiag << nl << nl;
        }

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Transpose Verifications:" << nl;

        const label numberOfTests = 10;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(1,10);
            const label nCols = rndGen.position(1,10);

            Info << "# SquareMatrix<scalar> example:" << nl;
            {
                Info<< "A(mxm) where m = " << mRows << nl;
                auto A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));
                auto B(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));

                Info<< nl << "# (A.T()).T() = A:" << nl;
                isEqual((A.T()).T(), A, true);

                Info<< nl << "# (A + B).T() = A.T() + B.T():" << nl;
                isEqual((A + B).T(), A.T() + B.T(), true);

                Info<< nl << "# (A*B).T() = B.T()*A.T():" << nl;
                isEqual((A*B).T(), B.T()*A.T(), true);
                Info<< nl << nl;
            }

            Info << "# RectangularMatrix<scalar> example:" << nl;
            {
                Info<< "A(mxn) where"
                    << " m = " << mRows << ","
                    << " n = " << nCols << nl;
                auto A(makeRandomMatrix<RMatrix>({mRows, nCols}, rndGen));
                auto B(makeRandomMatrix<RMatrix>({mRows, nCols}, rndGen));
                auto C(makeRandomMatrix<RMatrix>({nCols, mRows}, rndGen));

                Info<< nl << "# (A.T()).T() = A:" << nl;
                isEqual((A.T()).T(), A, true);

                Info<< nl << "# (A + B).T() = A.T() + B.T():" << nl;
                isEqual((A + B).T(), A.T() + B.T(), true);

                Info<< nl << "# (A*C).T() = C.T()*C.T():" << nl;
                isEqual((A*C).T(), C.T()*A.T(), true);
                Info<< nl << nl;
            }

            Info << "# SquareMatrix<complex> example:" << nl;
            {
                Info<< "A(mxm) where m = " << mRows << nl;
                SCMatrix A(mRows);
                SCMatrix B(mRows);

                for (auto& val : A)
                {
                    val.real(rndGen.GaussNormal<scalar>());
                    val.imag(rndGen.GaussNormal<scalar>());
                }
                for (auto& val : B)
                {
                    val.real(rndGen.GaussNormal<scalar>());
                    val.imag(rndGen.GaussNormal<scalar>());
                }

                Info<< nl << "# (A.T()).T() = A:" << nl;
                isEqual((A.T()).T(), A, true);

                Info<< nl << "# (A + B).T() = A.T() + B.T():" << nl;
                isEqual((A + B).T(), A.T() + B.T(), true);

                Info<< nl << "# (A*B).T() = B.T()*A.T():" << nl;
                isEqual((A*B).T(), B.T()*A.T(), true);
                Info<< nl << nl;
            }
        }

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Scalar Arithmetic Operations:" << nl;

        if (true)
        {
            SMatrix S1(3, Zero);
            assignMatrix
            (
                S1,
                {
                    4.1, 12.5, -16.3,
                    -192.3, -9.1, -3.0,
                    1.0, 5.02, -4.4
                }
            );
            SMatrix S2 = S1;
            SMatrix S3(3, Zero);

            Info<< "SquareMatrix<scalar> S1 = " << S1 << nl;
            Info<< "SquareMatrix<scalar> S2 = " << S2 << nl;
            Info<< "SquareMatrix<scalar> S3 = " << S3 << nl;

            S3 = S1*S2;
            Info<< "S1*S2 = " << S3 << nl;

            S3 = S1 - S2;
            Info<< "S1 - S2 = " << S3 << nl;

            S3 = S1 + S2;
            Info<< "S1 + S2 = " << S3 << nl;

//          S3 *= S1; // Not Implemented

            S3 -= S1;
            Info<< "S3 -= S1; S3 = " << S3 << nl;

            S3 += S1;
            Info<< "S3 += S1; S3 = " << S3 << nl;


            Info<< nl << "# Scalar broadcasting:" << nl;

            S3 *= 5.0;
            Info<< "S3 *= 5.0; S3 = " << S3 << nl;

            S3 /= 5.0;
            Info<< "S3 /= 5.0; S3 = " << S3 << nl;

            S3 -= 1.0;
            Info<< "S3 -= 1.0; S3 = " << S3 << nl;

            S3 += 1.0;
            Info<< "S3 += 1.0; S3 = " << S3 << nl;


            Info<< nl << "# Operations between different matrix types:" << nl;
            RMatrix R1 = S1;
            Info<< "RectangularMatrix<scalar> R1 = " << R1 << nl;

            // RectangularMatrix*SquareMatrix returns RectangularMatrix
            // S3 = S3*R1; // Not implemented
            R1 = S3*R1;
            Info<< "R1 = S3*R1; R1 = " << R1 << nl;

            S3 = S3 - R1;
            Info<< "S3 = S3 - R1; S3 = " << S3 << nl;

            S3 = S3 + R1;
            Info<< "S3 = S3 + R1; S3 = " << S3 << nl;

            R1 = S3 + R1;
            Info<< "R1 = S3 + R1; R1 = " << R1 << nl;


            Info<< nl << "# Extrema operations:" << nl;

            Info<< "min(S3) = " << min(S3) << nl
                << "max(S3) = " << max(S3) << nl
                << "minMax(S3) = " << minMax(S3) << nl;
        }

        if (true)
        {
            Info<< nl << "# Operations with ListTypes<scalar>:" << nl;

            SMatrix A(3, Zero);
            assignMatrix
            (
                A,
                {
                    4.1, 12.5, -16.3,
                    -192.3, -9.1, -3.0,
                    1.0, 5.02, -4.4
                }
            );
            const List<scalar> lst({1, 2, 3});
            RMatrix colVector(3, 1, Zero);
            assignMatrix
            (
                colVector,
                {1, 2, 3}
            );
            RMatrix rowVector(1, 3, Zero);
            assignMatrix
            (
                rowVector,
                {1, 2, 3}
            );

            Info<< "A = " << A << nl;
            Info<< "colVector = " << colVector << nl;
            Info<< "rowVector = " << rowVector << nl;

            const Field<scalar> field1(A.Amul(lst));
            const Field<scalar> field2(A*lst);
            const Field<scalar> field3(A.Tmul(lst));
            const Field<scalar> field4(lst*A);

            Info
                << "Field1 = A*lst = A.Amul(lst):" << nl << field1 << nl
                << "Field2 = A*lst:" << nl << field2 << nl
                << "A*colVector:" << nl << A*colVector << nl
                << "Field3 = lst*A = A.Tmul(lst):" << nl << field3 << nl
                << "Field4 = lst*A:" << nl << field4 << nl
                << "rowVector*A:" << nl << rowVector*A << nl
                << nl;
        }

        if (true)
        {
            Info<< nl << "# Implicit inner/outer products:" << nl;

            const scalar r = 2.0;
            RMatrix A(3, 2, Zero);
            assignMatrix
            (
                A,
                {
                    1, 2,
                    3, 4,
                    5, 6
                }
            );
            const RMatrix B(A);
            const RMatrix C(B);

            Info<< nl << "# Inner product:" << nl;
            {
                Info<< "- Explicit vs Implicit => A.T()*B == A&B:" << nl;
                isEqual(A.T()*B, A&B, true);

                Info<< "- Commutative => A&B == B&A:" << nl;
                isEqual(A&B, B&A, true);

                Info<< "- Distributive => A&(B+C) == A&B + A&C:" << nl;
                isEqual(A&(B + C), (A&B) + (A&C), true);

                Info<< "- Bilinear => A&(rB+C) == r(A&B) + (A&C):" << nl;
                isEqual(A&(r*B + C), r*(A&B) + (A&C), true);

                Info<< "- Scalar multiplication => (rA)&(rB) == rr(A&B):" << nl;
                isEqual((r*A)&(r*B), r*r*(A&B), true);
            }

            Info<< nl << "# Outer product:" << nl;
            {
                Info<< "- Explicit vs Implicit => A*B.T() == A^B:" << nl;
                isEqual(A*B.T(), A^B, true);

                Info<< "- Commutative => A^B == B^A:" << nl;
                isEqual(A^B, B^A, true);

                Info<< "- Distributive => A^(B+C) == A^B + A^C:" << nl;
                isEqual(A^(B + C), (A^B) + (A^C), true);

                Info<< "- Scalar multiplication => r*(A^B) == (r*A)^B:" << nl;
                isEqual(r*(A^B), (r*A)^B, true);
            }
        }

        Info<< nl;
        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Complex Arithmetic Operations:" << nl;

        if (true)
        {
            SCMatrix S1(3, Zero);
            assignMatrix
            (
                S1,
                {
                    complex(4.1, 1.0), complex(12.5, -1), complex(-16.3,-3),
                    complex(-192.3, 5), complex(-9.1, 3), complex(-3.0, 1),
                    complex(1.0, 3), complex(5.02, 0.3), complex(-4.4, 1)
                }
            );
            SCMatrix S2 = S1;
            SCMatrix S3(3, Zero);

            Info<< "SquareMatrix<complex> S1 = " << S1 << nl;
            Info<< "SquareMatrix<complex> S2 = " << S2 << nl;
            Info<< "SquareMatrix<complex> S3 = " << S3 << nl;

            S3 = S1*S2;
            Info<< "S1*S2 = " << S3 << nl;

            S3 = S1 - S2;
            Info<< "S1 - S2 = " << S3 << nl;

            S3 = S1 + S2;
            Info<< "S1 + S2 = " << S3 << nl;

//          S3 *= S1; // Not Implemented

            S3 -= S1;
            Info<< "S3 -= S1; S3 = " << S3 << nl;

            S3 += S1;
            Info<< "S3 += S1; S3 = " << S3 << nl;


            Info<< nl << "# Complex broadcasting:" << nl;

            S3 *= complex(5, 0);
            Info<< "S3 *= complex(5, 0); S3 = " << S3 << nl;

            S3 /= complex(5, 0);
            Info<< "S3 /= complex(5, 0); S3 = " << S3 << nl;

            S3 -= complex(1, 0);
            Info<< "S3 -= complex(1, 0); S3 = " << S3 << nl;

            S3 += complex(1, 0);
            Info<< "S3 += complex(1, 0); S3 = " << S3 << nl;


            Info<< nl << "# Operations between different matrix types:" << nl;
            RCMatrix R1 = S1;
            Info<< "RectangularMatrix<complex> R1 = " << R1 << nl;

            R1 = S3*R1;
            Info<< "R1 = S3*R1; R1 = " << R1 << nl;

            S3 = S3 - R1;
            Info<< "S3 = S3 - R1; S3 = " << S3 << nl;

            S3 = S3 + R1;
            Info<< "S3 = S3 + R1:" << S3 << nl;

            // Extrama operations // Not Implemented
        }

        if (true)
        {
            Info<< nl << "# Operations with ListTypes<complex>:" << nl;

            SCMatrix A(3, Zero);
            assignMatrix
            (
                A,
                {
                    complex(4.1, 1.0), complex(12.5, -1), complex(-16.3,-3),
                    complex(-192.3, 5), complex(-9.1, 3), complex(-3.0, 1),
                    complex(1.0, 3), complex(5.02, 0.3), complex(-4.4, 1)
                }
            );
            const List<complex> lst({complex(1,1), complex(2,2), complex(3,3)});
            RCMatrix colVector(3, 1, Zero);
            assignMatrix
            (
                colVector,
                {complex(1,1), complex(2,2), complex(3,3)}
            );
            RCMatrix rowVector(1, 3, Zero);
            assignMatrix
            (
                rowVector,
                {complex(1,1), complex(2,2), complex(3,3)}
            );

            Info<< "A = " << A << nl;
            Info<< "colVector = " << colVector << nl;
            Info<< "rowVector = " << rowVector << nl;

            const Field<complex> field1(A.Amul(lst));
            const Field<complex> field2(A*lst);
            const Field<complex> field3(A.Tmul(lst));
            const Field<complex> field4(lst*A);

            Info
                << "Field1 = A*lst = A.Amul(lst):" << nl << field1 << nl
                << "Field2 = A*lst:" << nl << field2 << nl
                << "A*colVector:" << nl << A*colVector << nl
                << "Field3 = lst*A = A.Tmul(lst):" << nl << field3 << nl
                << "Field4 = lst*A:" << nl << field4 << nl
                << "rowVector*A:" << nl << rowVector*A << nl
                << nl;
        }

        #if (0 | RUNALL)
        {
            Info<< nl << "# Implicit inner/outer products:" << nl;

            const complex r(2.0,1.0);
            RCMatrix A(3, 2, Zero);
            assignMatrix
            (
                A,
                {
                    complex(1,0), complex(2,1),
                    complex(3,-0.3), complex(4,-1),
                    complex(0.5,-0.1), complex(6,0.4)
                }
            );
            const RCMatrix B(A);
            const RCMatrix C(B);

            Info<< nl << "# Inner product:" << nl;
            {
                Info<< "- Explicit vs Implicit => A.T()*B == A&B:" << nl;
                isEqual(A.T()*B, A&B, true);

                Info<< "- Commutative => A&B == B&A:" << nl;
                isEqual(A&B, B&A, true);

                Info<< "- Distributive => A&(B+C) == A&B + A&C:" << nl;
                isEqual(A&(B + C), (A&B) + (A&C), true);

                Info<< "- Bilinear => A&(rB+C) == r(A&B) + (A&C):" << nl;
                isEqual(A&(r*B + C), r*(A&B) + (A&C), true);
            }

            Info<< nl << "# Outer product:" << nl;
            {
                Info<< "- Explicit vs Implicit => A*B.T() == A^B:" << nl;
                isEqual(A*B.T(), A^B, true);

                Info<< "- Commutative => A^B == B^A:" << nl;
                isEqual(A^B, B^A, true);

                Info<< "- Distributive => A^(B+C) == A^B + A^C:" << nl;
                isEqual(A^(B + C), (A^B) + (A^C), true);

                Info<< "- Scalar multiplication => r*(A^B) == (r*A)^B:" << nl;
                isEqual(r*(A^B), (r*A)^B, true);
            }
        }
        #endif

        Info<< nl;
        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## Query to determine if a square matrix is symmetric:" << nl;

        {
            const label mRows = 100;
            SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));

            Info<< "# SquareMatrix<scalar>.symmetric():" << nl
                << A.symmetric() << nl;

            // Symmetrise
            for (label n = 0; n < A.n() - 1; ++n)
            {
                for (label m = A.m() - 1; n < m; --m)
                {
                    A(n, m) = A(m, n);
                }
            }

            Info<< "# SquareMatrix<scalar>.symmetric():" << nl
                << A.symmetric() << nl;
        }

        {
            const label mRows = 100;

            SCMatrix A(mRows);

            for (auto& val : A)
            {
                val.real(rndGen.GaussNormal<scalar>());
                val.imag(rndGen.GaussNormal<scalar>());
            }

            Info<< "# SquareMatrix<complex>.symmetric():" << nl
                << A.symmetric() << nl;

            // Symmetrise
            for (label n = 0; n < A.n() - 1; ++n)
            {
                for (label m = A.m() - 1; n < m; --m)
                {
                    A(n, m) = A(m, n);
                }
            }

            Info<< "# SquareMatrix<complex>.symmetric():" << nl
                << A.symmetric() << nl;
        }

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## SymmetricSquareMatrix<scalar> algorithms:" << nl;

        scalarSymmetricSquareMatrix symm(3, Zero);

        symm(0, 0) = 4;
        symm(1, 0) = 12;
        symm(1, 1) = 37;
        symm(2, 0) = -16;
        symm(2, 1) = -43;
        symm(2, 2) = 98;

        Info<< "SymmetricSquareMatrix = " << nl << symm << nl
            << "Inverse = " << nl << inv(symm) << nl
            << "Determinant = " << det(symm) << nl;

        scalarSymmetricSquareMatrix symm2(symm);
        LUDecompose(symm2);
        Info<< "LU decomposition:" << nl
            << "Inverse = " << nl << invDecomposed(symm2) << nl
            << "Determinant = " << detDecomposed(symm2) << nl;

        scalarDiagonalMatrix rhs(3, 0);
        rhs[0] = 1;
        rhs[1] = 2;
        rhs[2] = 3;

        LUsolve(symm, rhs);
        Info<< "Solving linear system through LU decomposition:" << nl
            << "Decomposition = " << nl << symm << nl
            << "Solution = " << rhs << nl;
    }
    #endif


    #if (0 | RUNALL)
    {
        scalarSquareMatrix squareMatrix(3, Zero);

        scalarSquareMatrix S(3, Zero);
        assignMatrix
        (
            S,
            {
                4, 12, -16,
                12, 37, -43,
                -16, -43, 98
            }
        );

        const scalarField source(3, 1);

        Info<< nl << "Square Matrix = " << S << endl;

        if (true)
        {
            {
                scalarSquareMatrix S2(S);
                Info<< "Determinant = " << det(S2) << nl;
            }

            scalarSquareMatrix S2(S);
            labelList rhs(3, 0);
            label sign;
            LUDecompose(S2, rhs, sign);

            Info<< "LU decomposition = " << S2 << nl
                << "Pivots = " << rhs << nl
                << "Sign = " << sign << nl
                << "Determinant = " << detDecomposed(S2, sign) << nl;
        }

        if (true)
        {
            LUscalarMatrix LU(S);
            scalarField x(LU.solve(source));
            scalarSquareMatrix inv(3);
            LU.inv(inv);
            Info<< "LU.solve(source) residual: " << (S*x - source)
                << nl << "LU inv " << inv << nl
                << "LU inv*S " << (inv*S) << nl;
        }

        if (true)
        {
            LLTMatrix<scalar> LLT(S);
            scalarField x(LLT.solve(source));
            Info<< "LLT solve residual " << (S*x - source) << nl;
        }

        horizontalLine();
    }
    #endif


    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
