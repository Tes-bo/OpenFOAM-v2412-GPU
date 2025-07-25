/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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

Application
    Test-Tensor

Description
    Tests for \c Tensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c tensor, i.e. Tensor<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "tensor.H"
#include "transform.H"
#include "Random.H"
#include "scalar.H"
#include "complex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Create a random tensor
tensor makeRandomContainer(Random& rnd)
{
    tensor A(Zero);
    std::generate(A.begin(), A.end(), [&]{ return rnd.GaussNormal<scalar>(); });
    return A;
}


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
    const scalar absTol = 0,     //<! useful for cmps near zero
    const scalar relTol = 1e-8   //<! are values the same within 8 decimals
)
{
    Info<< msg << x << "?=" << y << endl;

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
    const scalar absTol = 0,
    const scalar relTol = 1e-8
)
{
    Info<< msg << x << "?=" << y << endl;

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


// Create each constructor of Tensor<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const Tensor<Type> T(Zero);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given MatrixSpace of the same rank:" << nl;
        const MatrixSpace<Tensor<Type>, Type, 3, 3> M(Zero);
        const Tensor<Type> T(M);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given VectorSpace of the same rank:" << nl;
        const VectorSpace<Tensor<Type>, Type, 9> V(Zero);
        const Tensor<Type> T(V);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given SphericalTensor<Type>:" << nl;
        const SphericalTensor<Type> Sp(Type(5));
        const Tensor<Type> T(Sp);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given SymmTensor<Type>:" << nl;
        const SymmTensor<Type> S
        (
            Type(1), Type(2), Type(3),
                     Type(5), Type(6),
                              Type(9)
        );
        const Tensor<Type> T(S);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given triad of row vectors,"
            << " optionally treated as transposed (ie, column vectors)" << nl;
        const Vector<Vector<Type>> vecs
        (
            Vector<Type>(Type(1), Type(2), Type(3)),
            Vector<Type>(Type(4), Type(5), Type(6)),
            Vector<Type>(Type(7), Type(8), Type(9))
        );
        const Tensor<Type> T(vecs);
        Info<< T << nl;

        const Tensor<Type> transposedT(vecs, true);
        Info<< transposedT << endl;
    }
    {
        Info<< "# Construct given the three row vectors,"
            << " optionally treated as transposed (ie, column vectors)" << nl;
        const Vector<Type> a(Type(1), Type(2), Type(3));
        const Vector<Type> b(Type(4), Type(5), Type(6));
        const Vector<Type> c(Type(7), Type(8), Type(9));
        const Tensor<Type> T(a, b, c);
        Info<< T << nl;

        const Tensor<Type> transposedT(a, b, c, true);
        Info<< transposedT << endl;
    }
    {
        Info<< "# Construct given the nine components:" << nl;
        const Tensor<Type> T
        (
            Type(1), Type(2), Type(-3),
            Type(4), Type(5), Type(-6),
            Type(7), Type(8), Type(-9)
        );
        Info<< T << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const Tensor<Type> T(Zero);
        const Tensor<Type> copyT(T);
        Info<< T << tab << copyT << endl;
    }
}


// Execute each member function of Tensor<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    Tensor<Type> T
    (
        Type(1), Type(2), Type(-3),
        Type(4), Type(5), Type(-6),
        Type(7), Type(8), Type(-9)
    );
    Tensor<Type> Tbak = T;
    const Tensor<Type> cT
    (
        Type(-9), Type(8), Type(7),
        Type(-6), Type(5), Type(4),
        Type(-3), Type(2), Type(1)
    );

    Info<< "# Operand: " << nl
        << "  Tensor = " << T << endl;


    {
        Info<< "# Component access:" << nl;

        Tensor<Type> cpT
        (
            T.xx(), T.xy(), T.xz(),
            T.yx(), T.yy(), T.yz(),
            T.zx(), T.zy(), T.zz()
        );
        cmp("  'Tensor' access:", T, cpT);

        const Tensor<Type> cpcT
        (
            cT.xx(), cT.xy(), cT.xz(),
            cT.yx(), cT.yy(), cT.yz(),
            cT.zx(), cT.zy(), cT.zz()
        );
        cmp("  'const Tensor' access:", cT, cpcT);
    }
    {
        Info<< "# Column-vector access:" << nl;
        cmp("  cx():", T.cx(), Vector<Type>(Type(1), Type(4), Type(7)));
        cmp("  cy():", T.cy(), Vector<Type>(Type(2), Type(5), Type(8)));
        cmp("  cz():", T.cz(), Vector<Type>(Type(-3), Type(-6), Type(-9)));
        cmp("  col(0):", T.col(0), Vector<Type>(Type(1), Type(4), Type(7)));
        cmp("  col(1):", T.col(1), Vector<Type>(Type(2), Type(5), Type(8)));
        cmp("  col(2):", T.col(2), Vector<Type>(Type(-3), Type(-6), Type(-9)));
        cmp
        (
            "  col<0>:",
            T.template col<0>(),
            Vector<Type>(Type(1), Type(4), Type(7))
        );
        cmp
        (
            "  col<1>:",
            T.template col<1>(),
            Vector<Type>(Type(2), Type(5), Type(8))
        );
        cmp
        (
            "  col<2>:",
            T.template col<2>(),
            Vector<Type>(Type(-3), Type(-6), Type(-9))
        );
        // Compilation error:  Info << "  col<3> = " << T.col<3>() << nl;


        Info<< "# Column-vector manipulation:" << nl;
        T.col(1, Vector<Type>(Type(0), Type(1), Type(99)));
        cmp
        (
            "  col(1, Vector):",
            T.col(1),
            Vector<Type>(Type(0), Type(1), Type(99))
        );

        T.cols
        (
            Vector<Type>(Type(1), Type(1), Type(1)),
            Vector<Type>(Type(-1), Type(1), Type(2)),
            Vector<Type>(Type(1), Type(1), Type(3))
        );
        cmp
        (
            "  cols(Vectors):",
            T,
            Tensor<Type>
            (
                Type(1), Type(-1), Type(1),
                Type(1), Type(1),  Type(1),
                Type(1), Type(2),  Type(3)
            )
        );
    }
    {
        Info<< "# Row-vector access:" << nl;
        T = Tbak;
        cmp("  x():", T.x(), Vector<Type>(Type(1), Type(2), Type(-3)));
        cmp("  y():", T.y(), Vector<Type>(Type(4), Type(5), Type(-6)));
        cmp("  z():", T.z(), Vector<Type>(Type(7), Type(8), Type(-9)));
        cmp("  row(0):", T.row(0), Vector<Type>(Type(1), Type(2), Type(-3)));
        cmp("  row(1):", T.row(1), Vector<Type>(Type(4), Type(5), Type(-6)));
        cmp("  row(2):", T.row(2), Vector<Type>(Type(7), Type(8), Type(-9)));
        cmp
        (
            "  row<0>:",
            T.template row<0>(),
            Vector<Type>(Type(1), Type(2), Type(-3))
        );
        cmp
        (
            "  row<1>:",
            T.template row<1>(),
            Vector<Type>(Type(4), Type(5), Type(-6))
        );
        cmp
        (
            "  row<2>:",
            T.template row<2>(),
            Vector<Type>(Type(7), Type(8), Type(-9))
        );
        // Compilation error:  Info << "  row<3> = " << T.row<3>() << nl;


        Info<< "# Row-vector manipulation:" << nl;
        T.row(1, Vector<Type>(Type(0), Type(1), Type(99)));
        cmp
        (
            "  row(1, Vector):",
            T.row(1),
            Vector<Type>(Type(0), Type(1), Type(99))
        );

        T.rows
        (
            Vector<Type>(Type(1), Type(1), Type(1)),
            Vector<Type>(Type(-1), Type(1), Type(2)),
            Vector<Type>(Type(1), Type(1), Type(3))
        );
        cmp
        (
            "  rows(Vectors):",
            T,
            Tensor<Type>
            (
                Type(1), Type(1), Type(1),
                Type(-1), Type(1), Type(2),
                Type(1), Type(1), Type(3)
            )
        );
    }
    {
        Info<< "# Diagonal access:" << nl;

        T = Tbak;
        cmp
        (
            "  'Tensor'.diag():",
            T.diag(),
            Vector<Type>(Type(1), Type(5), Type(-9))
        );
        cmp
        (
            "  'const Tensor'.diag():",
            cT.diag(),
            Vector<Type>(Type(-9), Type(5), Type(1))
        );


        Info<< "# Diagonal manipulation:" << nl;

        T.diag(Vector<Type>(Type(-10), Type(-15), Type(-20)));
        cmp
        (
            "  'Tensor'.diag('Vector'):",
            T.diag(),
            Vector<Type>(Type(-10), Type(-15), Type(-20))
        );
    }
    {
        Info<< "# Tensor operations:" << nl;

        T = Tbak;
        cmp("  Transpose:", T, (T.T()).T());
        cmp    // Singular matrix
        (
            "  Inverse:",
            T.inv(),
            Tensor<Type>
            (
            Type(-4.50359963e+15), Type(9.00719925e+15), Type(-4.50359963e+15),
            Type(9.00719925e+15), Type(-1.80143985e+16), Type(9.00719925e+15),
            Type(4.50359963e+15), Type(-9.00719925e+15), Type(4.50359963e+15)
            )
        );
        cmp
        (
            "  Inner-product:",
            T.inner(T),
            Tensor<Type>
            (
                Type(-12), Type(-12), Type(12),
                Type(-18), Type(-15), Type(12),
                Type(-24), Type(-18), Type(12)
            )
        );
        cmp
        (
            "  Schur-product:",
            T.schur(T),
            Tensor<Type>
            (
                Type(1),  Type(4),  Type(9),
                Type(16), Type(25), Type(36),
                Type(49), Type(64), Type(81)
            )
        );
    }
    {
        Info<< "# Member operators:" << nl;

        T = Tbak;
        T &= T;
        cmp
        (
            "  Assign inner-product of this with another Tensor:",
            T,
            (Tbak & Tbak)
        );

        T = VectorSpace<Tensor<Type>, Type, 9>(Zero);
        cmp
        (
            "  Assign to an equivalent vector space:",
            T,
            Tensor<Type>(Zero)
        );

        T = SphericalTensor<Type>(Type(5));
        cmp
        (
            "  Assign to a SphericalTensor:",
            T,
            Tensor<Type>
            (
                Type(5), Zero,    Zero,
                Zero,    Type(5), Zero,
                Zero,    Zero,    Type(5)
            )
        );

        T = SymmTensor<Type>
        (
            Type(1), Type(2), Type(-3),
                     Type(5), Type(-6),
                              Type(-9)
        );
        cmp
        (
            "  Assign to a SymmTensor:",
            T,
            Tensor<Type>
            (
                Type(1),  Type(2),  Type(-3),
                Type(2),  Type(5),  Type(-6),
                Type(-3), Type(-6), Type(-9)
            )
        );

        T = Vector<Vector<Type>>
            (
                Vector<Type>(Type(-1), Type(2), Type(3)),
                Vector<Type>(Type(-4), Type(5), Type(6)),
                Vector<Type>(Type(4), Type(-5), Type(6))
            );
        cmp
        (
            "  Assign to a triad of row vectors:",
            T,
            Tensor<Type>
            (
                Type(-1), Type(2), Type(3),
                Type(-4), Type(5), Type(6),
                Type(4), Type(-5), Type(6)
            )
        );
    }
}


// Execute each global function of Tensor<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const Tensor<Type> T
    (
        Type(-1), Type(2), Type(-3),
        Type(4), Type(5), Type(-6),
        Type(7), Type(8), Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(-1), Type(2), Type(-3),
                  Type(5), Type(-6),
                           Type(-9)
    );

    Info<< "# Operands: " << nl
        << "  Tensor = " << T << nl
        << "  SymmTensor = " << sT << endl;


    cmp("  Trace = ", tr(T), Type(-5));
    cmp("  Spherical part = ", sph(T), SphericalTensor<Type>(tr(T)/Type(3)));
    cmp
    (
        "  Symmetric part = ",
        symm(T),
        SymmTensor<Type>
        (
            Type(-1), Type(3), Type(2),
                      Type(5), Type(1),
                               Type(-9)
        )
    );
    cmp
    (
        "  Twice the symmetric part = ",
        twoSymm(T),
        SymmTensor<Type>
        (
            Type(-2), Type(6),  Type(4),
                      Type(10), Type(2),
                                Type(-18)
        )
    );
    cmp
    (
        "  Skew-symmetric part = ",
        skew(T),
        Tensor<Type>
        (
            Type(0),  Type(-1), Type(-5),
            Type(1),  Type(0),  Type(-7),
            Type(5),  Type(7),  Type(0)
        )
    );
    /*
    // Complex-type is not supported for this function.
    cmp
    (
        "  Skew-symmetric part of a SymmTensor = ",
        skew(sT),
        Tensor<Type>(Zero)
    );
    */
    cmp
    (
        "  Deviatoric part = ",
        dev(T),
        Tensor<Type>
        (
            Type(0.66666667), Type(2),          Type(-3),
            Type(4),          Type(6.66666667), Type(-6),
            Type(7),          Type(8),          Type(-7.33333333)
        ),
        1e-7
    );
    cmp("  Two-third deviatoric part = ", dev2(T), T - 2*sph(T));
    cmp("  Determinant = ", det(T), Type(-6.000000000000005));
    cmp
    (
        "  Cofactor tensor = ",
        cof(T),
        Tensor<Type>
        (
            Type(3),  Type(-6),  Type(-3),
            Type(-6), Type(30),  Type(22),
            Type(3),  Type(-18), Type(-13)
        )
    );
    cmp
    (
        "  Inverse = ",
        inv(T, det(T)),
        Tensor<Type>
        (
            Type(-0.5),  Type(1),            Type(-0.5),
            Type(1),     Type(-5),           Type(3),
            Type(0.5),   Type(-3.66666667),  Type(2.16666667)
        ),
        1e-8
    );
    cmp
    (
        "  Inverse (another) = ",
        inv(T),
        Tensor<Type>
        (
            Type(-0.5),  Type(1),            Type(-0.5),
            Type(1),     Type(-5),           Type(3),
            Type(0.5),   Type(-3.66666667),  Type(2.16666667)
        ),
        1e-8
    );
    cmp
    (
        "  Inverse (another) = ",
        T.inv(),
        Tensor<Type>
        (
            Type(-0.5),  Type(1),            Type(-0.5),
            Type(1),     Type(-5),           Type(3),
            Type(0.5),   Type(-3.66666667),  Type(2.16666667)
        ),
        1e-8
    );
    cmp("  First invariant = ", invariantI(T), Type(-5));
    cmp("  Second invariant = ", invariantII(T), Type(20));
    cmp("  Third invariant = ", invariantIII(T), Type(-6.000000000000005));
}


// Execute each global operator of Tensor<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor<Type> T
    (
        Type(-1), Type(2), Type(-3),
        Type(4),  Type(5),  Type(-6),
        Type(7),  Type(8),  Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(-1), Type(2), Type(-3),
                 Type(5), Type(-6),
                          Type(-9)
    );
    const SphericalTensor<Type> spT(Type(1));
    const Vector<Type> v(Type(3), Type(2), Type(1));
    const Type x(4);

    Info<< "# Operands:" << nl
        << "  Tensor = " << T << nl
        << "  SymmTensor = " << sT << nl
        << "  SphericalTensor = " << spT << nl
        << "  Vector = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Sum of SpTensor-Tensor = ",
        (spT + T),
        Tensor<Type>
        (
            Type(0), Type(2), Type(-3),
            Type(4), Type(6), Type(-6),
            Type(7), Type(8), Type(-8)
        )
    );
    cmp
    (
        "  Sum of Tensor-SpTensor = ",
        (T + spT),
        Tensor<Type>
        (
            Type(0), Type(2), Type(-3),
            Type(4), Type(6), Type(-6),
            Type(7), Type(8), Type(-8)
        )
    );
    cmp
    (
        "  Sum of SymmTensor-Tensor = ",
        (sT + T),
        Tensor<Type>
        (
            Type(-2), Type(4),  Type(-6),
            Type(6),  Type(10), Type(-12),
            Type(4),  Type(2),  Type(-18)
        )
    );
    cmp
    (
        "  Sum of Tensor-SymmTensor = ",
        (T + sT),
        Tensor<Type>
        (
            Type(-2), Type(4),  Type(-6),
            Type(6),  Type(10), Type(-12),
            Type(4),  Type(2),  Type(-18)
        )
    );
    cmp
    (
        "  Subtract Tensor from SpTensor = ",
        (spT - T),
        Tensor<Type>
        (
            Type(2),  Type(-2), Type(3),
            Type(-4), Type(-4), Type(6),
            Type(-7), Type(-8), Type(10)
        )
    );
    cmp
    (
        "  Subtract SpTensor from Tensor = ",
        (T - spT),
        Tensor<Type>
        (
            Type(-2), Type(2), Type(-3),
            Type(4),  Type(4), Type(-6),
            Type(7),  Type(8), Type(-10)
        )
    );
    cmp
    (
        "  Subtract Tensor from SymmTensor = ",
        (sT - T),
        Tensor<Type>
        (
            Type(0),   Type(0),   Type(0),
            Type(-2),  Type(0),   Type(0),
            Type(-10), Type(-14), Type(0)
        )
    );
    cmp
    (
        "  Subtract SymmTensor from Tensor = ",
        (T - sT),
        Tensor<Type>
        (
            Type(0),  Type(0),  Type(0),
            Type(2),  Type(0),  Type(0),
            Type(10), Type(14), Type(0)
        )
    );
    cmp
    (
        "  Hodge dual of a Tensor = ",
        *T,
        Vector<Type>(T.yz(), -T.xz(), T.xy())
    );
    cmp
    (
        "  Hodge dual of a Vector = ",
        *v,
        Tensor<Type>
        (
             Zero,  -v.z(),  v.y(),
             v.z(),  Zero,  -v.x(),
            -v.y(),  v.x(),  Zero
        )
    );
    /*cmp
    (
        "  Division of Vector by Tensor = ",
        (v/T),
        Tensor<Type>
        (
            Type(-3),          Type(1),    Type(-0.33333333),
            Type(0.75),        Type(0.4),  Type(-0.16666667),
            Type(0.42857143),  Type(0.25), Type(-0.11111111)
        )
    );*/
    cmp
    (
        "  Division of Tensor by Type = ",
        (T/x),
        Tensor<Type>
        (
            Type(-0.25), Type(0.5),  Type(-0.75),
            Type(1),     Type(1.25), Type(-1.5),
            Type(1.75),  Type(2),    Type(-2.25)
        )
    );
    cmp
    (
        "  Inner-product of Tensor-Tensor = ",
        (T & T),
        Tensor<Type>
        (
            Type(-12), Type(-16), Type(18),
            Type(-26), Type(-15), Type(12),
            Type(-38), Type(-18), Type(12)
        )
    );
    cmp
    (
        "  Inner-product of SpTensor-Tensor = ",
        (spT & T),
        Tensor<Type>
        (
            Type(-1), Type(2),  Type(-3),
            Type(4),  Type(5),  Type(-6),
            Type(7),  Type(8),  Type(-9)
        )
    );
    cmp
    (
        "  Inner-product of Tensor-SpTensor = ",
        (T & spT),
        Tensor<Type>
        (
            Type(-1), Type(2), Type(-3),
            Type(4),  Type(5), Type(-6),
            Type(7),  Type(8), Type(-9)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor-Tensor = ",
        (sT & T),
        Tensor<Type>
        (
            Type(-12), Type(-16),  Type(18),
            Type(-24), Type(-19),  Type(18),
            Type(-84), Type(-108), Type(126)
        )
    );
    cmp
    (
        "  Inner-product of Tensor-SymmTensor = ",
        (T & sT),
        Tensor<Type>
        (
            Type(14), Type(26),  Type(18),
            Type(24), Type(69),  Type(12),
            Type(36), Type(108), Type(12)
        )
    );
    cmp
    (
        "  Inner-product of Tensor-Vector = ",
        (T & v),
        Vector<Type>(Type(-2), Type(16),  Type(28)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector-Tensor = ",
        (v & T),
        Vector<Type>(Type(12), Type(24),  Type(-30)) // Row-vector
    );
    cmp("  D-inner-product of SpTensor-Tensor = ", (spT && T), Type(-5));
    cmp("  D-inner-product of Tensor-SpTensor = ", (T && spT), Type(-5));
    cmp("  D-inner-product of SymmTensor-Tensor = ", (sT && T), Type(95));
    cmp("  D-inner-product of Tensor-SymmTensor = ", (T && sT), Type(95));
    cmp
    (
        "  Outer-product of Vector-Vector = ",
        (v*v),
        Tensor<Type>
        (
            Type(9), Type(6), Type(3),
            Type(6), Type(4), Type(2),
            Type(3), Type(2), Type(1)
        )
    );
}


// Return false if given eigenvalues fail to satisy eigenvalue relations
// Relations: (Beauregard & Fraleigh (1973), ISBN 0-395-14017-X, p. 307)
void test_eigenvalues
(
    const tensor& T,
    const Vector<complex>& EVals,
    const bool prod = true
)
{
    if (prod)
    {
        const scalar determinant = det(T);
        // In case of complex EVals, the production is effectively scalar
        // due to the (complex*complex conjugate) results in zero imag part
        const scalar EValsProd = ((EVals.x()*EVals.y()*EVals.z()).real());
        cmp("# Product of eigenvalues = det(T):", EValsProd, determinant, 1e-8);
    }

    {
        const scalar trace = tr(T);
        scalar EValsSum = 0.0;
        // In case of complex EVals, the summation is effectively scalar
        // due to the (complex+complex conjugate) results in zero imag part
        for (const auto& val : EVals)
        {
            EValsSum += val.real();
        }
        cmp("# Sum of eigenvalues = trace(T):", EValsSum, trace);
    }
}


// Return false if a given eigenvalue-eigenvector pair
// fails to satisfy the characteristic equation
void test_characteristic_equation
(
    const tensor& T,
    const Vector<complex>& EVals,
    const Tensor<complex>& EVecs
)
{
    Info<< "# Characteristic equation:" << nl;

    Tensor<complex> Tc(Zero);
    forAll(T, i)
    {
        Tc[i] = complex(T[i], 0);
    }

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        const Vector<complex> leftSide(Tc & EVecs.row(dir));
        const Vector<complex> rightSide(EVals[dir]*EVecs.row(dir));
        const Vector<complex> X(leftSide - rightSide);

        for (const auto x : X)
        {
            cmp("  (sT & EVec - EVal*EVec) = 0:", mag(x), 0.0, 1e-5);
        }
    }
}


// Return false if the eigen functions fail to satisfy relations
void test_eigen_funcs(const tensor& T, const bool prod = true)
{
    Info<< "# Operand:" << nl
        << "  tensor = " << T << nl;


    Info<< "# Return eigenvalues of a given tensor:" << nl;
    const Vector<complex> EVals(eigenValues(T));
    Info<< EVals << endl;
    test_eigenvalues(T, EVals, prod);

    Info<< "# Return an eigenvector of a given tensor in a given direction"
        << " corresponding to a given eigenvalue:" << nl;
    const Vector<complex> standardBasis1(Zero, pTraits<complex>::one, Zero);
    const Vector<complex> standardBasis2(Zero, Zero, pTraits<complex>::one);
    const Vector<complex> EVec
    (
        eigenVector(T, EVals.x(), standardBasis1, standardBasis2)
    );
    Info<< EVec << endl;

    Info<< "# Return eigenvectors of a given tensor corresponding to"
        << " given eigenvalues:" << nl;
    const Tensor<complex> EVecs0(eigenVectors(T, EVals));
    Info<< EVecs0 << endl;
    test_characteristic_equation(T, EVals, EVecs0);

    Info<< "# Return eigenvectors of a given tensor by computing"
        << " the eigenvalues of the tensor in the background:" << nl;
    const Tensor<complex> EVecs1(eigenVectors(T));
    Info<< EVecs1 << endl;
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
        "Tensor<floatScalar>",
        "Tensor<doubleScalar>",
        "Tensor<complex>"
    });

    run_tests(types, typeID);


    Info<< nl << "    ## Test tensor eigen functions: ##" << nl;
    const label numberOfTests = 10000;
    Random rndGen(1234);

    for (label i = 0; i < numberOfTests; ++i)
    {
        const tensor T(makeRandomContainer(rndGen));
        test_eigen_funcs(T);
    }

    {
        Info<< nl << "    ## A zero tensor: ##"<< nl;
        const tensor zeroT(Zero);
        test_eigen_funcs(zeroT);
    }
    {
        Info<< nl
            << "    ## A skew-symmetric tensor with no-real eigenvalues: ##"
            << nl;
        const tensor T
        (
            0,  1,  1,
           -1,  0,  1,
           -1, -1,  0
        );
        test_eigen_funcs(T);
    }
    {
        Info<< nl
            << "    ## A stiff tensor: ##"
            << nl;
        const tensor stiff
        (
            pow(10.0, 10), pow(10.0, 8), pow(10.0, 7),
            pow(10.0, -8), pow(10.0, 9), pow(10.0, -8),
            pow(10.0, 10), pow(10.0, 8), pow(10.0, 7)
        );
        // Although eigendecomposition is successful for the stiff tensors,
        // cross-check between prod(eigenvalues) ?= det(stiff) is inherently
        // problematic; therefore, eigenvalues of the stiff tensors are
        // cross-checked by only sum(eigenvalues) ?= trace(stiff)
        const bool testProd = false;
        test_eigen_funcs(stiff, testProd);
    }
    {
        Info<< nl
            << "    ## Random tensor with tiny off-diag elements: ##"
            << nl;

        const List<scalar> epsilons
        ({
            0, SMALL, Foam::sqrt(SMALL), sqr(SMALL), Foam::cbrt(SMALL),
            -SMALL, -Foam::sqrt(SMALL), -sqr(SMALL), -Foam::cbrt(SMALL)
        });

        for (label i = 0; i < numberOfTests; ++i)
        {
            for (const auto& eps : epsilons)
            {
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    T.zx() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    T.zx() = eps*rndGen.GaussNormal<scalar>();
                    T.zy() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = 0;
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = 0;
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    T.zx() = eps*rndGen.GaussNormal<scalar>();
                    T.zy() = 0;
                    test_eigen_funcs(T);
                }
                {
                    tensor T(makeRandomContainer(rndGen));
                    T.xy() = eps;
                    T.xz() = eps;
                    T.yz() = eps;
                    T.yx() = eps;
                    T.zx() = eps;
                    T.zy() = eps;
                    test_eigen_funcs(T);
                }
            }
        }
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
