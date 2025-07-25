/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Test-Tensor2D

Description
    Tests for \c Tensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c tensor2D, i.e. Tensor2D<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "vector2DField.H"
#include "tensor2D.H"
#include "symmTensor2D.H"
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


// Create a random tensor2D
tensor2D makeRandomContainer(Random& rnd)
{
    tensor2D A(Zero);
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


// Create each constructor of Tensor2D<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const Tensor2D<Type> T(Zero);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given VectorSpace:" << nl;
        const VectorSpace<Tensor2D<Type>, Type, 4> V(Zero);
        const Tensor2D<Type> T(V);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given SymmTensor2D:" << nl;
        const SymmTensor2D<Type> S
        (
            Type(1), Type(2),
                     Type(3)
        );
        const Tensor2D<Type> T(S);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given SphericalTensor2D:" << nl;
        const SphericalTensor2D<Type> Sp(Type(5));
        const Tensor2D<Type> T(Sp);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given the two row vectors:" << nl;
        const Vector2D<Type> x(Type(1), Type(2));
        const Vector2D<Type> y(Type(3), Type(4));
        const Tensor2D<Type> T(x, y);
        Info<< T << endl;
    }
    {
        Info<< "# Construct given the four components:" << nl;
        const Tensor2D<Type> T
        (
            Type(1), Type(2),
            Type(3), Type(4)
        );
        Info<< T << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const Tensor2D<Type> T(Zero);
        const Tensor2D<Type> Tcopy(T);
        Info<< T << endl;
    }
}


// Execute each member function of Tensor2D<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    Tensor2D<Type> T
    (
        Type(1), Type(2),
        Type(4), Type(5)
    );
    Tensor2D<Type> Tbak = T;
    const Tensor2D<Type> cT
    (
        Type(-9), Type(8),
        Type(-6), Type(5)
    );

    Info<< "# Operand: " << nl
        << "  Tensor2D = " << T << endl;


    {
        Info<< "# Component access:" << nl;

        Tensor2D<Type> cpT
        (
            T.xx(), T.xy(),
            T.yx(), T.yy()
        );
        cmp("  'Tensor2D' access:", T, cpT);

        const Tensor2D<Type> cpcT
        (
            cT.xx(), cT.xy(),
            cT.yx(), cT.yy()
        );
        cmp("  'const Tensor2D' access:", cT, cpcT);
    }
    {
        Info<< "# Column-vector access:" << nl;
        cmp("  cx():", T.cx(), Vector2D<Type>(Type(1), Type(4)));
        cmp("  cy():", T.cy(), Vector2D<Type>(Type(2), Type(5)));
        cmp("  col(0):", T.col(0), Vector2D<Type>(Type(1), Type(4)));
        cmp("  col(1):", T.col(1), Vector2D<Type>(Type(2), Type(5)));
        cmp
        (
            "  col<0>:",
            T.template col<0>(),
            Vector2D<Type>(Type(1), Type(4))
        );
        cmp
        (
            "  col<1>:",
            T.template col<1>(),
            Vector2D<Type>(Type(2), Type(5))
        );
        // Compilation error:  Info << "  col<2> = " << T.col<2>() << nl;


        Info<< "# Column-vector manipulation:" << nl;
        T.col(1, Vector2D<Type>(Type(0), Type(1)));
        cmp
        (
            "  col(1, Vector):",
            T.col(1),
            Vector2D<Type>(Type(0), Type(1))
        );

        T.cols
        (
            Vector2D<Type>(Type(1), Type(1)),
            Vector2D<Type>(Type(-1), Type(1))
        );
        cmp
        (
            "  cols(Vectors):",
            T,
            Tensor2D<Type>
            (
                Type(1), Type(-1),
                Type(1), Type(1)
            )
        );
    }
    {
        Info<< "# Row-vector access:" << nl;
        T = Tbak;
        cmp("  x():", T.x(), Vector2D<Type>(Type(1), Type(2)));
        cmp("  y():", T.y(), Vector2D<Type>(Type(4), Type(5)));
        cmp("  row(0):", T.row(0), Vector2D<Type>(Type(1), Type(2)));
        cmp("  row(1):", T.row(1), Vector2D<Type>(Type(4), Type(5)));
        cmp
        (
            "  row<0>:",
            T.template row<0>(),
            Vector2D<Type>(Type(1), Type(2))
        );
        cmp
        (
            "  row<1>:",
            T.template row<1>(),
            Vector2D<Type>(Type(4), Type(5))
        );
        // Compilation error:  Info << "  row<2> = " << T.row<2>() << nl;


        Info<< "# Row-vector manipulation:" << nl;
        T.row(1, Vector2D<Type>(Type(0), Type(1)));
        cmp
        (
            "  row(1, Vector):",
            T.row(1),
            Vector2D<Type>(Type(0), Type(1))
        );

        T.rows
        (
            Vector2D<Type>(Type(1), Type(1)),
            Vector2D<Type>(Type(-1), Type(1))
        );
        cmp
        (
            "  rows(Vectors):",
            T,
            Tensor2D<Type>
            (
                Type(1), Type(1),
                Type(-1), Type(1)
            )
        );
    }
    {
        Info<< "# Diagonal access:" << nl;

        T = Tbak;
        cmp
        (
            "  'Tensor2D'.diag():",
            T.diag(),
            Vector2D<Type>(Type(1), Type(5))
        );
        cmp
        (
            "  'const Tensor2D'.diag():",
            cT.diag(),
            Vector2D<Type>(Type(-9), Type(5))
        );


        Info<< "# Diagonal manipulation:" << nl;

        T.diag(Vector2D<Type>(Type(-10), Type(-15)));
        cmp
        (
            "  'Tensor2D'.diag('Vector'):",
            T.diag(),
            Vector2D<Type>(Type(-10), Type(-15))
        );
    }
    {
        Info<< "# Tensor operations:" << nl;

        T = Tbak;
        cmp("  Transpose:", T, (T.T()).T());
        cmp
        (
            "  Inner-product:",
            T.inner(T),
            Tensor2D<Type>
            (
                Type(9),  Type(12),
                Type(24), Type(33)
            )
        );
        cmp
        (
            "  Schur-product:",
            T.schur(T),
            Tensor2D<Type>
            (
                Type(1),  Type(4),
                Type(16), Type(25)
            )
        );
    }
    {
        Info<< "# Member operators:" << nl;

        T = SphericalTensor2D<Type>(Type(5));
        cmp
        (
            "  Assign to a SphericalTensor2D:",
            T,
            Tensor2D<Type>
            (
                Type(5), Zero,
                Zero,    Type(5)
            )
        );

        T = SymmTensor2D<Type>
        (
            Type(1), Type(2),
                     Type(5)
        );
        cmp
        (
            "  Assign to a SymmTensor2D:",
            T,
            Tensor2D<Type>
            (
                Type(1),  Type(2),
                Type(2),  Type(5)
            )
        );
    }
}


// Execute each global function of Tensor2D<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const Tensor2D<Type> T
    (
        Type(-1), Type(2),
        Type(4), Type(5)
    );

    Info<< "# Operand: " << nl
        << "  Tensor2D = " << T << endl;


    cmp("  Trace = ", tr(T), Type(4));
    cmp("  Spherical part = ", sph(T), SphericalTensor2D<Type>(tr(T)/Type(2)));
    cmp
    (
        "  Symmetric part = ",
        symm(T),
        SymmTensor2D<Type>
        (
            Type(-1), Type(3),
                      Type(5)
        )
    );
    cmp
    (
        "  Twice the symmetric part = ",
        twoSymm(T),
        SymmTensor2D<Type>
        (
            Type(-2), Type(6),
                      Type(10)
        )
    );
    cmp
    (
        "  Skew-symmetric part = ",
        skew(T),
        Tensor2D<Type>
        (
            Type(0),  Type(-1),
            Type(1),  Type(0)
        )
    );
    cmp
    (
        "  Deviatoric part = ",
        dev(T),
        Tensor2D<Type>
        (
            Type(-3), Type(2),
            Type(4),  Type(3)
        )
    );
    cmp("  Two-third deviatoric part = ", dev2(T), T - 2*sph(T));
    cmp("  Determinant = ", det(T), Type(-13));
    cmp
    (
        "  Cofactor tensor2D = ",
        cof(T),
        Tensor2D<Type>
        (
            Type(5),  Type(-4),
            Type(-2), Type(-1)
        )
    );
    cmp
    (
        "  Inverse = ",
        inv(T, det(T)),
        Tensor2D<Type>
        (
            Type(-0.38461538),  Type(0.15384615),
            Type(0.30769231),   Type(0.07692308)
        ),
        1e-6,
        1e-6
    );
    cmp
    (
        "  Inverse (another) = ",
        inv(T),
        Tensor2D<Type>
        (
            Type(-0.38461538),  Type(0.15384615),
            Type(0.30769231),   Type(0.07692308)
        ),
        1e-6,
        1e-6
    );
    cmp("  First invariant = ", invariantI(T), Type(4));
    cmp("  Second invariant = ", invariantII(T), Type(-13));
}


// Execute each global operator of Tensor2D<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor2D<Type> T
    (
        Type(-1), Type(2),
        Type(4),  Type(5)
    );
    const SymmTensor2D<Type> sT
    (
        Type(1), Type(2),
                 Type(5)
    );
    const SphericalTensor2D<Type> spT(Type(1));
    const Vector2D<Type> v(Type(3), Type(2));
    const Type x(4);

    Info<< "# Operands:" << nl
        << "  Tensor2D = " << T << nl
        << "  SymmTensor2D = " << sT << nl
        << "  SphericalTensor2D = " << spT << nl
        << "  Vector2D = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Sum of SpTensor2D-Tensor2D = ",
        (spT + T),
        Tensor2D<Type>
        (
            Type(0), Type(2),
            Type(4), Type(6)
        )
    );
    cmp
    (
        "  Sum of Tensor2D-SpTensor2D = ",
        (T + spT),
        Tensor2D<Type>
        (
            Type(0), Type(2),
            Type(4), Type(6)
        )
    );
    cmp
    (
        "  Sum of SymmTensor2D-Tensor2D = ",
        (sT + T),
        Tensor2D<Type>
        (
            Type(0), Type(4),
            Type(6),  Type(10)
        )
    );
    cmp
    (
        "  Sum of Tensor2D-SymmTensor2D = ",
        (T + sT),
        Tensor2D<Type>
        (
            Type(0), Type(4),
            Type(6),  Type(10)
        )
    );
    cmp
    (
        "  Subtract Tensor2D from SpTensor2D = ",
        (spT - T),
        Tensor2D<Type>
        (
            Type(2),   Type(-2),
            Type(-4),  Type(-4)
        )
    );
    cmp
    (
        "  Subtract SpTensor2D from Tensor2D = ",
        (T - spT),
        Tensor2D<Type>
        (
            Type(-2), Type(2),
            Type(4),  Type(4)
        )
    );
    cmp
    (
        "  Subtract Tensor2D from SymmTensor2D = ",
        (sT - T),
        Tensor2D<Type>
        (
            Type(2),  Type(0),
            Type(-2), Type(0)
        )
    );
    cmp
    (
        "  Subtract SymmTensor2D from Tensor2D = ",
        (T - sT),
        Tensor2D<Type>
        (
            Type(-2), Type(0),
            Type(2), Type(0)
        )
    );
    cmp
    (
        "  Division of Tensor2D by Type = ",
        (T/x),
        Tensor2D<Type>
        (
            Type(-0.25), Type(0.5),
            Type(1),     Type(1.25)
        )
    );
    cmp
    (
        "  Inner-product of Tensor2D-Tensor2D = ",
        (T & T),
        Tensor2D<Type>
        (
            Type(9),  Type(8),
            Type(16), Type(33)
        )
    );
    cmp
    (
        "  Inner-product of SpTensor2D-Tensor2D = ",
        (spT & T),
        Tensor2D<Type>
        (
            Type(-1), Type(2),
            Type(4),  Type(5)
        )
    );
    cmp
    (
        "  Inner-product of Tensor2D-SpTensor2D = ",
        (T & spT),
        Tensor2D<Type>
        (
            Type(-1), Type(2),
            Type(4),  Type(5)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor2D-Tensor2D = ",
        (sT & T),
        Tensor2D<Type>
        (
            Type(7),  Type(12),
            Type(18), Type(29)
        )
    );
    cmp
    (
        "  Inner-product of Tensor2D-SymmTensor2D = ",
        (T & sT),
        Tensor2D<Type>
        (
            Type(3),  Type(8),
            Type(14), Type(33)
        )
    );
    cmp
    (
        "  Inner-product of Tensor2D-Vector2D = ",
        (T & v),
        Vector2D<Type>(Type(1), Type(22)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector2D-Tensor2D = ",
        (v & T),
        Vector2D<Type>(Type(5), Type(16)) // Row-vector
    );
    cmp("  D-inner-product of SpTensor2D-Tensor2D = ", (spT && T), Type(4));
    cmp("  D-inner-product of Tensor2D-SpTensor2D = ", (T && spT), Type(4));
    cmp("  D-inner-product of SymmTensor2D-Tensor2D = ", (sT && T), Type(36));
    cmp("  D-inner-product of Tensor2D-SymmTensor2D = ", (T && sT), Type(36));
    cmp
    (
        "  Outer-product of Vector2D-Vector2D = ",
        (v*v),
        Tensor2D<Type>
        (
            Type(9), Type(6),
            Type(6), Type(4)
        )
    );
}


// Return false if given eigenvalues fail to satisy eigenvalue relations
// Relations: (Beauregard & Fraleigh (1973), ISBN 0-395-14017-X, p. 307)
void test_eigenvalues(const tensor2D& T, const Vector2D<complex>& EVals)
{
    {
        const scalar determinant = det(T);
        // In case of complex EVals, the production is effectively scalar
        // due to the (complex*complex conjugate) results in zero imag part
        const scalar EValsProd = ((EVals.x()*EVals.y()).real());
        cmp("# Product of eigenvalues = det(T):", EValsProd, determinant, 1e-7);
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
        cmp("# Sum of eigenvalues = trace(T):", EValsSum, trace, 1e-8);
    }
}


// Return false if a given eigenvalue-eigenvector pair
// fails to satisfy the characteristic equation
void test_characteristic_equation
(
    const tensor2D& T,
    const Vector2D<complex>& EVals,
    const Tensor2D<complex>& EVecs
)
{
    Info<< "# Characteristic equation:" << nl;

    Tensor2D<complex> Tc(Zero);
    forAll(T, i)
    {
        Tc[i] = complex(T[i], 0);
    }

    for (direction dir = 0; dir < pTraits<vector2D>::nComponents; ++dir)
    {
        const Vector2D<complex> leftSide(Tc & EVecs.row(dir));
        const Vector2D<complex> rightSide(EVals[dir]*EVecs.row(dir));
        const Vector2D<complex> X(leftSide - rightSide);

        for (const auto x : X)
        {
            cmp("  (Tc & EVec - EVal*EVec) = 0:", mag(x), 0.0, 1e-5);
        }
    }
}


// Return false if the eigen functions fail to satisfy relations
void test_eigen_funcs(const tensor2D& T)
{
    Info<< "# Operand:" << nl
        << "  tensor2D = " << T << nl;


    Info<< "# Return eigenvalues of a given tensor2D:" << nl;
    const Vector2D<complex> EVals(eigenValues(T));
    Info<< EVals << endl;
    test_eigenvalues(T, EVals);

    Info<< "# Return an eigenvector of a given tensor2D in a given direction"
        << " corresponding to a given eigenvalue:" << nl;
    const Vector2D<complex> standardBasis(pTraits<complex>::one, Zero);
    const Vector2D<complex> EVec(eigenVector(T, EVals.x(), standardBasis));
    Info<< EVec << endl;

    Info<< "# Return eigenvectors of a given tensor2D corresponding to"
        << " given eigenvalues:" << nl;
    const Tensor2D<complex> EVecs0(eigenVectors(T, EVals));
    Info<< EVecs0 << endl;
    test_characteristic_equation(T, EVals, EVecs0);

    Info<< "# Return eigenvectors of a given tensor2D by computing"
        << " the eigenvalues of the tensor2D in the background:" << nl;
    const Tensor2D<complex> EVecs1(eigenVectors(T));
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

int main(int argc, char *argv[])
{

    const std::tuple<floatScalar, doubleScalar, complex> types
    (
        std::make_tuple(Zero, Zero, Zero)
    );

    const List<word> typeID
    ({
        "Tensor2D<floatScalar>",
        "Tensor2D<doubleScalar>",
        "Tensor2D<complex>"
    });

    run_tests(types, typeID);


    Info<< nl << "    ## Test tensor2D eigen functions: ##" << nl;
    const label numberOfTests = 10000;
    Random rndGen(1234);

    for (label i = 0; i < numberOfTests; ++i)
    {
        const tensor2D T(makeRandomContainer(rndGen));
        test_eigen_funcs(T);
    }

    {
        Info<< nl << "    ## A zero tensor2D: ##"<< nl;
        const tensor2D zeroT(Zero);
        test_eigen_funcs(zeroT);
    }
    {
        Info<< nl
            << "    ## A tensor2D with repeated eigenvalues: ##"
            << nl;
        const tensor2D T
        (
            -1.0, 2.0,
            0.0, -1.0
        );
        test_eigen_funcs(T);
    }
    {
        Info<< nl
            << "    ## A skew-symmetric tensor2D with no-real eigenvalues: ##"
            << nl;
        const tensor2D T
        (
            0.0, 1.0,
            -1.0, 0.0
        );
        test_eigen_funcs(T);
    }
    {
        Info<< nl
            << "    ## A stiff tensor2D: ##"
            << nl;
        const tensor2D stiff
        (
            pow(10.0, 10), pow(10.0, 8),
            pow(10.0, -8), pow(10.0, 9)
        );
        test_eigen_funcs(stiff);
    }
    {
        Info<< nl
            << "    ## Random tensor2D with tiny off-diag elements: ##"
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
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = 0;
                    T.yx() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.yx() = 0;
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps;
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.yx() = eps;
                    test_eigen_funcs(T);
                }
                {
                    tensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps;
                    T.yx() = eps;
                    test_eigen_funcs(T);
                }
            }
        }
    }


    {
        Info<< "# Pre-v2006 tests:" << nl;
        vector2D v1(1, 2), v2(3, 4);
        tensor2D t3 = v1*v2;

        Info<< v1 << "*" << v2 << " = " << t3 << endl;

        {
            Info<< "rows:" << nl;
            for (direction i = 0; i < 2; ++i)
            {
                Info<< "  (" << i << ") = " << t3.row(i) << nl;
            }
        }

        {
            Info<< "cols:" << nl;
            for (direction i = 0; i < 2; ++i)
            {
                Info<< "  (" << i << ") = " << t3.col(i) << nl;
            }
            Info<< "col<0> = " << t3.col<0>() << nl;
            Info<< "col<1> = " << t3.col<1>() << nl;
            // Compilation error:  Info << "col<3> = " << t3.col<3>() << nl;

            t3.col<0>({0, 2});
            Info<< "replaced col<0> = " << t3.col<0>() << nl;
            Info<< "tensor " << t3 << nl;

            t3.row<1>(Zero);
            Info<< "replaced row<1> = " << t3.row<1>() << nl;
            Info<< "tensor " << t3 << nl;
        }


        {
            vector2DField vfld1(8, Zero);

            forAll(vfld1, i)
            {
                vfld1[i] = (i + 1) * ((i % 2) ? v1 : v2);
            }

            Info<< "vector: " << flatOutput(vfld1) << nl;

            scalarField xvals(8);
            scalarField yvals(8);
            unzip(vfld1, xvals, yvals);

            Info<< "unzip" << nl
                << " x => " << flatOutput(xvals) << nl
                << " y => " << flatOutput(yvals) << nl;

            reverse(xvals);
            zip(vfld1, xvals, yvals);

            Info<< "rezip (with reversed x)" << nl
                << "   => " << flatOutput(vfld1) << nl;
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


// ************************************************************************* //
