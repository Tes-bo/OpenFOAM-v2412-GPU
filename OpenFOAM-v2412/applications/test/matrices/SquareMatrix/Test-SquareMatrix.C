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
    Test-SquareMatrix

Description
    Tests for \c SquareMatrix constructors, member functions, member
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
#include "scalar.H"
#include "complex.H"
#include "IOmanip.H"
#include "TestTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create each constructor of SquareMatrix<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct for given size (rows == cols):" << nl;
        const SquareMatrix<Type> A(2);
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing all elements to zero:" << nl;
        const SquareMatrix<Type> A(2, Zero);
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing all elements to a value:" << nl;
        const SquareMatrix<Type> A(2, Type(3));
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) "
            << "initializing to the identity matrix:" << nl;
        const SquareMatrix<Type> A(2, I);
        Info<< A << endl;
    }
    {
        Info<< "# Construct for given size (rows == cols) by using a labelPair"
            << "initializing to the identity matrix:" << nl;
        const SquareMatrix<Type> A(labelPair(2, 2), I);
        Info<< A << endl;
    }
    {
        Info<< "# Construct given number of rows/columns "
            << "by using a label pair (checked to be equal):" << nl;
        const SquareMatrix<Type> A(labelPair(2, 2));
        Info<< A << endl;
    }
    {
        Info<< "# Construct given number of rows/columns "
            << "by using a label pair (checked to be equal) "
            << "and initializing all elements to zero:" << nl;
        const SquareMatrix<Type> A(labelPair(2, 2), Zero);
        Info<< A << endl;
    }
    {
        Info<< "# Construct given number of rows/columns "
            << "by using a label pair (checked to be equal) "
            << "and initializing all elements to the given value:" << nl;
        const SquareMatrix<Type> A(labelPair(2, 2), Type(3));
        Info<< A << endl;
    }
    {
        Info<< "# Construct given number of rows/columns "
            << "(checked to be equal) initializing all elements to zero" << nl;
        const SquareMatrix<Type> A(2, 2, Zero);
        Info<< A << endl;
    }
    {
        Info<< "# Construct from const sub-matrix block:" << nl;
        const SquareMatrix<Type> B(labelPair(3, 3), Type(3));
        const SquareMatrix<Type> A(B.subMatrix(1, 1));
        Info<< A << endl;
    }
    {
        Info<< "# Construct from sub-matrix block:" << nl;
        SquareMatrix<Type> B(labelPair(3, 3), Type(3));
        const SquareMatrix<Type> A(B.subMatrix(1, 1));
        Info<< A << endl;
    }
    {
        Info<< "#: Construct as copy of a RectangularMatrix:" << nl;
        const RectangularMatrix<Type> B(2, 2, Zero);
        const SquareMatrix<Type> A(B);
        Info<< A << endl;
    }

    //  ## Pending ##
    //  Construct from Istream
    //  Clone
    //  #############
}


// Execute each member function of SquareMatrix<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SquareMatrix<Type> A(labelPair(3, 3), Zero);
    assignMatrix
    (
        A,
        {
            Type(1), Type(-2.2), Type(-3.4),
            Type(-0.35), Type(1), Type(5),
            Type(-0.35), Type(1), Type(5)
        }
    );

    Info<< "# Operand: " << nl
        << "  SquareMatrix = " << A << endl;


    {
        Info<< "# Access:" << nl;

        cmp("  The number of rows = ", Type(A.m()), Type(3));
        cmp("  The number of columns = ", Type(A.n()), Type(3));
        cmp("  The number of elements in Matrix = ", Type(A.size()), Type(9));
        cmp("  Return row/column sizes = ", A.sizes(), labelPair(3, 3));
        cmp("  Return true if Matrix is empty = ", A.empty(), false);
        cmp
        (
            "  Return const pointer to the first data elem = ",
            *(A.cdata()),
            Type(1)
        );
        cmp
        (
            "  Return pointer to the first data elem = ",
            *(A.data()),
            Type(1)
        );
        cmp
        (
            "  Return const pointer to data in the specified row = ",
            *(A.rowData(1)),
            Type(-0.35)
        );
        cmp
        (
            "  Return pointer to data in the specified row = ",
            *(A.rowData(1)),
            Type(-0.35)
        );
        const Type& a = A.at(4);
        cmp("  Linear addressing const element access = ", a, Type(1));
        Type& b = A.at(4);
        cmp("  Linear addressing element access = ", b, Type(1));
    }
    {
        Info<< "# Block access (const):" << nl;

        const RectangularMatrix<Type> Acol(A.subColumn(1));
        cmp
        (
            "  Return const column or column's subset of Matrix = ",
            flt(Acol),
            List<Type>({Type(-2.2), Type(1), Type(1)})
        );

        const RectangularMatrix<Type> Arow(A.subRow(1));
        cmp
        (
            "  Return const row or row's subset of Matrix = ",
            flt(Arow),
            List<Type>({Type(-0.35), Type(1), Type(5)})
        );

        const RectangularMatrix<Type> Amat(A.subMatrix(1, 1));
        cmp
        (
            "  Return const sub-block of Matrix = ",
            flt(Amat),
            List<Type>({Type(1), Type(5), Type(1), Type(5)})
        );
    }
    {
        Info<< "# Block access (non-const):" << nl;

        RectangularMatrix<Type> Acol(A.subColumn(1));
        cmp
        (
            "  Return column or column's subset of Matrix = ",
            flt(Acol),
            List<Type>({Type(-2.2), Type(1), Type(1)})
        );

        RectangularMatrix<Type> Arow(A.subRow(1));
        cmp
        (
            "  Return row or row's subset of Matrix = ",
            flt(Arow),
            List<Type>({Type(-0.35), Type(1), Type(5)})
        );

        RectangularMatrix<Type> Amat(A.subMatrix(1, 1));
        cmp
        (
            "  Return sub-block of Matrix = ",
            flt(Amat),
            List<Type>({Type(1), Type(5), Type(1), Type(5)})
        );
    }
    {
        Info<< "# Check:" << nl;

        A.checki(0);
        A.checkj(1);
        A.checkSize();
        cmp("  Check all entries have identical values = ", A.uniform(), false);

        {
            Random rndGen(1234);
            const label mRows = 10;
            SquareMatrix<Type> B
            (
                makeRandomMatrix<SquareMatrix<Type>>(labelPair(mRows, mRows), rndGen)
            );
            // Symmetrise
            for (label n = 0; n < B.n() - 1; ++n)
            {
                for (label m = B.m() - 1; n < m; --m)
                {
                    B(n, m) = B(m, n);
                }
            }
            cmp
            (
                "  Return true if the square matrix is "
                "effectively symmetric/Hermitian",
                B.symmetric(),
                true
            );
        }
        //  ## Pending ##
        //  Return true if the square matrix is reduced tridiagonal
        //  #############
    }
    {
        Info<< "# Edit:" << nl;

        SquareMatrix<Type> cpA0(A);
        cpA0.clear();
        cmp("  Clear Matrix, i.e. set sizes to zero = ", cpA0.empty(), true);

        SquareMatrix<Type> cpA1(A);
        cmp
        (
            "  Release storage management of Matrix contents by transferring "
            "management to a List = ",
            (cpA1.release()),
            List<Type>
            ({
                Type(1), Type(-2.2), Type(-3.4),
                Type(-0.35), Type(1), Type(5),
                Type(-0.35), Type(1), Type(5)
            })
        );

        SquareMatrix<Type> cpA2(A);
        cpA2.swap(cpA1);
        cmp("  Swap contents = ", flt(cpA1), flt(A));

        cpA2.transfer(cpA1);
        cmp
        (
            "  Transfer the contents of the argument Matrix into this Matrix "
            "and annul the argument MatrixSwap contents = ",
            flt(cpA2),
            flt(A)
        );

        cpA2.resize(2, 2);
        cmp
        (
            "  Change Matrix dimensions, preserving the elements = ",
            flt(cpA2),
            List<Type>({Type(1), Type(-2.2), Type(-0.35), Type(1)})
        );

        cpA2.resize(1);
        cmp
        (
            "  Resize the matrix preserving the elements = ",
            flt(cpA2),
            List<Type>({Type(1)})
        );

        SquareMatrix<Type> cpA3(A);
        cpA3.setSize(2);
        cmp
        (
            "  Resize the matrix preserving the elements = ",
            flt(cpA3),
            List<Type>({Type(1), Type(-2.2), Type(-0.35), Type(1)})
        );

        SquareMatrix<Type> cpA4(A);
        cpA4.shallowResize(2);
        cmp
        (
            "  Resize the matrix without reallocating storage (unsafe) = ",
            flt(cpA4),
            List<Type>({Type(1), Type(-2.2), Type(-3.4), Type(-0.35)})
        );

        SquareMatrix<Type> smallA(2, Type(VSMALL));
        smallA.round();
        cmp
        (
            "  Round elements with mag smaller than tol (SMALL) to zero = ",
            flt(smallA),
            List<Type>(4, Zero)
        );
    }
    {
        Info<< "# Sort:" << nl;

        SquareMatrix<Type> cpA(A);
        auto descend = [&](Type a, Type b){ return mag(a) > mag(b); };
        const List<label> sortPermutation(cpA.sortPermutation(descend));
        cmp
        (
            "  Return a sort permutation labelList according to "
            "a given comparison on the diagonal entries",
            sortPermutation,
            List<label>({2, 0, 1})
        );

        SquareMatrix<Type> sortedA0(3, Zero);
        assignMatrix
        (
            sortedA0,
            {
                Type(-3.4), Type(1),     Type(-2.2),
                Type(5),    Type(-0.35), Type(1),
                Type(5),    Type(-0.35), Type(1)
            }
        );
        const SquareMatrix<Type> sortedA1(applyPermutation(A, sortPermutation));
        cmp
        (
            "  Return Matrix column-reordered according to "
            "a given permutation labelList",
            flt(sortedA0),
            flt(sortedA1)
        );

        cpA.applyPermutation(sortPermutation);
        cmp
        (
            "  Column-reorder this Matrix according to "
            "a given permutation labelList",
            flt(sortedA0),
            flt(cpA)
        );
    }
    {
        Info<< "# Operations:" << nl;

        cmp("  Transpose = ", flt((A.T()).T()), flt(A));

        SquareMatrix<complex> cA(labelPair(2, 2), Zero);
        assignMatrix
        (
            cA,
            {
                complex(2, -1.2), complex(4.1, 0.4),
                complex(8.1, -1.25), complex(7.3, -1.4)
            }
        );
        SquareMatrix<complex> cAT(labelPair(2, 2), Zero);
        assignMatrix
        (
            cAT,
            {
                complex(2, 1.2), complex(8.1, 1.25),
                complex(4.1, -0.4), complex(7.3, 1.4)
            }
        );
        cmp("  Hermitian transpose = ", flt(cA.T()), flt(cAT));

        SquareMatrix<Type> B(3, 3, Zero);
        assignMatrix
        (
            B,
            {
                Type(4.1), Type(12.5), Type(-16.3),
                Type(-192.3), Type(-9.1), Type(-3.0),
                Type(1.0), Type(5.02), Type(-4.4)
            }
        );
        const SquareMatrix<Type> cpB(B);
        const List<Type> lst({Type(1), Type(2), Type(3)});
        const Field<Type> out1(B*lst);
        const Field<Type> out2(B.Amul(lst));
        const Field<Type> out3(lst*B);
        const Field<Type> out4(B.Tmul(lst));
        const Field<Type> rsult1({Type(-19.8), Type(-219.5), Type(-2.16)});
        const Field<Type> rsult2({Type(-377.5), Type(9.36), Type(-35.5)});
        cmp
        (
            "  Right-multiply Matrix by a List (A * x) = ",
            out1,
            rsult1,
            1e-4
        );
        cmp
        (
            "  Right-multiply Matrix by a column vector A.Amul(x) = ",
            out1,
            out2,
            1e-4
        );
        cmp
        (
            "  Left-multiply Matrix by a List (x * A) = ",
            out3,
            rsult2,
            1e-4
        );
        cmp
        (
            "  Left-multiply Matrix by a row vector A.Tmul(x) = ",
            out4,
            rsult2,
            1e-4
        );

        List<Type> diagB1({Type(4.1), Type(-9.1), Type(-4.4)});
        cmp
        (
            "  Extract the diagonal elements = ",
            B.diag(),
            diagB1
        );

        List<Type> diagB2({Type(-100), Type(-100), Type(-100)});
        B.diag(diagB2);
        cmp
        (
            "  Assign diagonal of Matrix = ",
            B.diag(),
            diagB2
        );

        B = cpB;
        cmp("  Trace = ", B.trace(), Type(-9.4), 1e-4);
        cmp
        (
            "  Return L2-Norm of chosen column = ",
            B.columnNorm(0),
            192.34630227794867,
            1e-4
        );
        cmp
        (
            "  Return L2-Norm of chosen column (noSqrt=true) = ",
            B.columnNorm(0, true),
            36997.1,
            1e-2
        );
        cmp
        (
            "  Return Frobenius norm of Matrix = ",
            B.norm(),
            193.7921835369012,
            1e-4
        );


    }
}


// Execute each member operators of SquareMatrix<Type>, and print output
template<class Type>
void test_member_opers(Type)
{
    SquareMatrix<Type> A(3, I);
    const SquareMatrix<Type> cpA(A);

    Info<< "# Operand: " << nl
        << "  SquareMatrix = " << A << endl;


    A = Zero;
    cmp("  Assign all elements to zero =", flt(A), List<Type>(9, Zero));
    A = cpA;

    A = Type(5);
    cmp("  Assign all elements to value =", flt(A), List<Type>(9, Type(5)));
    A = cpA;

    A = I;
    cmp("  Set to identity matrix =", flt(A), flt(cpA));
    A = cpA;

    const Type* a = A[1];
    cmp
    (
        "  Return const pointer to data in the specified row = ",
        *a,
        Type(0)
    );

    Type* b = A[1];
    cmp
    (
        "  Return pointer to data in the specified row = ",
        *b,
        Type(0)
    );

    const Type& c = A(1, 1);
    cmp
    (
        "  (i, j) const element access operator = ",
        c,
        Type(1)
    );

    Type& d = A(1, 1);
    cmp
    (
        "  (i, j) element access operator = ",
        d,
        Type(1)
    );

    //  ## Pending ##
    //  Copy assignment
    //  Move assignment
    //  Assignment to a block of another Matrix
    //  Assignment to a block of another Matrix
    //  #############

    A = Zero;
    cmp
    (
        "  Assignment of all elements to zero = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Zero))
    );
    A = cpA;

    A = Type(-1.2);
    cmp
    (
        "  Assignment of all elements to the given value = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Type(-1.2)))
    );
    A = cpA;

    A += cpA;
    cmp
    (
        "  Matrix addition =",
        flt(A),
        flt(Type(2)*cpA)
    );
    A = cpA;

    A -= cpA;
    cmp
    (
        "  Matrix subtraction = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Zero))
    );
    A = cpA;

    A = Zero;
    A += Type(5);
    cmp
    (
        "  Matrix scalar addition = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Type(5)))
    );
    A = cpA;

    A = Zero;
    A -= Type(5);
    cmp
    (
        "  Matrix scalar subtraction = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Type(-5)))
    );
    A = cpA;

    A = Type(1);
    A *= Type(5);
    cmp
    (
        "  Matrix scalar multiplication = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Type(5)))
    );
    A = cpA;

    A = Type(1);
    A /= Type(5);
    cmp
    (
        "  Matrix scalar division = ",
        flt(A),
        flt(SquareMatrix<Type>(3, Type(0.2)))
    );
    A = cpA;

    A(0,0) = Type(-10.5);
    Type* i0 = A.begin();
    cmp
    (
        "  Return an iterator to begin traversing a Matrix = ",
        *i0,
        Type(-10.5)
    );

    A(2,2) = Type(20);
    Type* i1 = A.end();
    cmp
    (
        "  Return an iterator to end traversing a Matrix = ",
        *(--i1),
        Type(20)
    );

    const Type* i2 = A.cbegin();
    cmp
    (
        "  Return const iterator to begin traversing a Matrix = ",
        *i2,
        Type(-10.5)
    );

    const Type* i3 = A.cend();
    cmp
    (
        "  Return const iterator to end traversing a Matrix = ",
        *(--i3),
        Type(20)
    );

    const Type* i4 = A.begin();
    cmp
    (
        "  Return const iterator to begin traversing a Matrix = ",
        *i4,
        Type(-10.5)
    );

    const Type* i5 = A.end();
    cmp
    (
        "  Return const iterator to end traversing a Matrix = ",
        *(--i5),
        Type(20)
    );

    //  ## Pending ##
    //  Read Matrix from Istream, discarding existing contents
    //  Write Matrix, with line-breaks in ASCII when length exceeds shortLen
    //  #############
}


// Execute each friend function of SquareMatrix<Type>, and print output
template<class Type>
void test_friend_funcs(Type)
{
    const Field<Type> F
    ({
        Type(1), Type(-2.2),
        Type(10), Type(-0.1)
    });

    Info<< "# Operand: " << nl
        << "  Field = " << F << endl;

    {
        SquareMatrix<Type> A(labelPair(4, 4), Zero);
        assignMatrix
        (
            A,
            {
                Type(1), Type(-2.2), Type(10), Type(-0.1),
                Type(-2.2), Type(4.84), Type(-22), Type(0.22),
                Type(10), Type(-22), Type(100), Type(-1),
                Type(-0.1), Type(0.22), Type(-1), Type(0.01)
            }
        );
        cmp
        (
            "  Return the outer product of Field-Field as RectangularMatrix = ",
            flt(outer(F, F)),
            flt(A),
            1e-6
        );
        cmp
        (
            "  Return the outer product of Field-Field as SquareMatrix = ",
            flt(symmOuter(F, F)),
            flt(A),
            1e-6
        );
    }
}


// Execute each global function of SquareMatrix<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    SquareMatrix<Type> A(3, Zero);
    assignMatrix
    (
        A,
        {
            Type(1), Type(-2.2), Type(-3.4),
            Type(-0.35), Type(10.32), Type(5),
            Type(-0.35), Type(10.32), Type(-5)
        }
    );
    const SquareMatrix<Type> cpA(A);

    Info<< "# Operand: " << nl
        << "  SquareMatrix = " << A << endl;

    //  ## Pending ##: 'LUDecompose' does not operate with 'complex' base type
    //  cmp
    //  (
    //      "  Return the determinant of SquareMatrix",
    //      det(A),
    //      -95.5
    //  );
    //  #############

    //  ## Pending ##
    //  cmp("  Find max value in Matrix = ", max(A), Type(10.32));
    //  cmp("  Find min value in Matrix = ", min(A), Type(-3.4));
    //  cmp
    //  (
    //      "  Find min-max value in Matrix = ",
    //      minMax(A),
    //      MinMax<Type>(10.32, -3.4)
    //  );
    //  #############
}


// Execute each global operators of SquareMatrix<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    SquareMatrix<Type> A(3, Zero);
    assignMatrix
    (
        A,
        {
            Type(1), Type(-2.2), Type(-3.4),
            Type(-0.35), Type(10.32), Type(5),
            Type(-0.35), Type(10.32), Type(-5)
        }
    );
    const SquareMatrix<Type> cpA(A);

    Info<< "# Operand: " << nl
        << "  SquareMatrix = " << A << endl;


    cmp("  Matrix negation = ", flt(-A), flt(Type(-1)*A));
    cmp("  Matrix addition = ", flt(A + A), flt(Type(2)*A));
    cmp("  Matrix subtraction = ", flt(A - A), flt(Type(0)*A));
    cmp("  Scalar multiplication = ", flt(Type(2)*A), flt(A + A));
    cmp("  Scalar multiplication = ", flt(A*Type(2)), flt(A + A));
    cmp
    (
        "  Scalar addition = ",
        flt(Type(0)*A + Type(1)),
        flt(SquareMatrix<Type>(A.sizes(), Type(1)))
    );
    cmp
    (
        "  Scalar addition = ",
        flt(Type(1) + Type(0)*A),
        flt(SquareMatrix<Type>(A.sizes(), Type(1)))
    );
    cmp
    (
        "  Scalar subtraction = ",
        flt(Type(0)*A - Type(1)),
        flt(SquareMatrix<Type>(A.sizes(), Type(-1)))
    );
    cmp
    (
        "  Scalar subtraction = ",
        flt(Type(1) - Type(0)*A),
        flt(SquareMatrix<Type>(A.sizes(), Type(1)))
    );
    cmp
    (
        "  Scalar division = ",
        flt((Type(1) + Type(0)*A)/Type(2.5)),
        flt(SquareMatrix<Type>(A.sizes(), Type(0.4)))
    );

    SquareMatrix<Type> A1(3, Zero);
    assignMatrix
    (
        A1,
        {
            Type(2.9600), Type(-59.9920), Type(2.6000),
            Type(-5.71199), Type(158.87239), Type(27.789997),
            Type(-2.211999), Type(55.672393), Type(77.789993)
        }
    );
    cmp
    (
        "  Matrix-Matrix multiplication using ikj-algorithm = ",
        flt(A*A),
        flt(A1),
        1e-4
    );
    SquareMatrix<Type> A2(3, Zero);
    assignMatrix
    (
        A2,
        {
            Type(1.245), Type(-9.424), Type(-3.4),
            Type(-9.424), Type(217.8448), Type(7.47999),
            Type(-3.4), Type(7.47999), Type(61.56)
        }
    );
    cmp
    (
        "  Implicit A.T()*B = ",
        flt(A&A),
        flt(A2),
        1e-4
    );
    SquareMatrix<Type> A3(3, Zero);
    assignMatrix
    (
        A3,
        {
            Type(17.4), Type(-40.054), Type(-6.054),
            Type(-40.054), Type(131.6249), Type(81.6249),
            Type(-6.054), Type(81.6249), Type(131.6249)
        }
    );
    cmp
    (
        "  Implicit A*B.T() = ",
        flt(A^A),
        flt(A3),
        1e-4
    );

    const List<Type> lst({Type(1), Type(2), Type(-3)});
    const List<Type> out1(A*lst);
    const List<Type> out2(lst*A);
    const List<Type> rslt1({Type(6.7999), Type(5.2900), Type(35.29)});
    const List<Type> rslt2({Type(1.35), Type(-12.52), Type(21.6)});
    cmp
    (
        "  Matrix-vector multiplication (A * x), where x is a column vector = ",
        out1,
        rslt1,
        1e-2
    );
    cmp
    (
        "  Vector-matrix multiplication (x * A), where x is a row vector = ",
        out2,
        rslt2,
        1e-4
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
        "SquareMatrix<floatScalar>",
        "SquareMatrix<doubleScalar>",
        "SquareMatrix<complex>"
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
