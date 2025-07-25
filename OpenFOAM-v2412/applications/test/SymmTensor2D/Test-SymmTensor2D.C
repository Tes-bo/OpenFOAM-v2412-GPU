/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Test-SymmTensor2D

Description
    Tests for \c SymmTensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c symmTensor2D, i.e. SymmTensor2D<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

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


// Create a random symmTensor2D
symmTensor2D makeRandomContainer(Random& rnd)
{
    symmTensor2D A(Zero);
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


// Create each constructor of SymmTensor2D<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const SymmTensor2D<Type> sT(Zero);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given VectorSpace:" << nl;
        const VectorSpace<SymmTensor2D<Type>, Type, 3> V(Zero);
        const SymmTensor2D<Type> sT(V);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given SphericalTensor2D:" << nl;
        const SphericalTensor2D<Type> Sp(Type(5));
        const SymmTensor2D<Type> sT(Sp);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given the three components:" << nl;
        const SymmTensor2D<Type> sT
        (
            Type(1), Type(2),
                     Type(4)
        );
        Info<< sT << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const SymmTensor2D<Type> S(Type(1), Type(2), Type(3));
        const SymmTensor2D<Type> sT(S);
        Info<< sT << endl;
    }
}


// Execute each member function of SymmTensor2D<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SymmTensor2D<Type> sT(Type(1), Type(2), Type(-3));
    const SymmTensor2D<Type> csT(Type(-3), Type(2), Type(1));


    {
        Info<< "# Component access:" << nl;

        SymmTensor2D<Type> cpsT
        (
            sT.xx(), sT.xy(),
                     sT.yy()
        );
        cmp("  'SymmTensor2D' access:", sT, cpsT);
        cmp("  xy()=yx():", sT.xy(), sT.yx());

        const SymmTensor2D<Type> cpcsT
        (
            csT.xx(), csT.xy(),
                      csT.yy()
        );
        cmp("  'const SymmTensor2D' access:", csT, cpcsT);
        cmp("  xy()=yx():", sT.xy(), sT.yx());
    }
    {
        Info<< "# Diagonal access:" << nl;

        cmp
        (
            "  'SymmTensor2D'.diag():",
            sT.diag(),
            Vector2D<Type>(Type(1), Type(-3))
        );
        cmp
        (
            "  'const SymmTensor2D'.diag():",
            csT.diag(),
            Vector2D<Type>(Type(-3), Type(1))
        );


        Info<< "# Diagonal manipulation:" << nl;

        sT.diag(Vector2D<Type>(Type(-10), Type(-15)));
        cmp
        (
            "  'SymmTensor2D'.diag('Vector2D'):",
            sT.diag(),
            Vector2D<Type>(Type(-10), Type(-15))
        );
    }
    {
        Info<< "# Tensor operations:" << nl;

        Info<< "  Transpose:" << nl;
        cmp("  'SymmTensor2D'.T():", sT.T(), sT);
    }
    {
        Info<< "# Member operators:" << nl;

        sT = SphericalTensor2D<Type>(Type(5));
        cmp
        (
            "  Assign to a SphericalTensor2D:",
            sT,
            SymmTensor2D<Type>
            (
                Type(5),  Zero,
                          Type(5)
            )
        );
    }
}


// Execute each global function of SymmTensor2D<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const SymmTensor2D<Type> sT(Type(1), Type(2), Type(-3));
    const Vector2D<Type> v(Type(-3), Type(2));

    Info<< "# Operands: " << nl
        << "  SymmTensor2D<Type> = " << sT << nl
        << "  Vector2D<Type> = " << v << endl;


    cmp("  Trace = ", tr(sT), Type(-2));
    cmp("  Spherical part = ", sph(sT), SphericalTensor2D<Type>(tr(sT)/Type(2)));
    cmp("  Symmetric part = ", symm(sT), sT);
    cmp("  Twice the symmetric part = ", twoSymm(sT), 2*sT);
    cmp
    (
        "  Deviatoric part = ",
        dev(sT),
        SymmTensor2D<Type>
        (
            Type(2), Type(2),
                     Type(-2)
        )
    );
    cmp("  Two-third deviatoric part = ", dev2(sT), sT - 2*sph(sT));
    cmp("  Determinant = ", det(sT), Type(-7.000000000000001));
    cmp
    (
        "  Cofactor tensor = ",
        cof(sT),
        SymmTensor2D<Type>
        (
            Type(-3), Type(-2),
                      Type(1)
        )
    );
    cmp
    (
        "  Inverse = ",
        inv(sT, det(sT)),
        SymmTensor2D<Type>
        (
            Type(0.42857143), Type(0.28571429),
                              Type(-0.14285714)
        ),
        1e-7
    );
    cmp
    (
        "  Inverse (another) = ",
        inv(sT),
        SymmTensor2D<Type>
        (
            Type(0.42857143), Type(0.28571429),
                              Type(-0.14285714)
        ),
        1e-7
    );
    cmp("  First invariant = ", invariantI(sT), Type(-2));
    cmp("  Second invariant = ", invariantII(sT), Type(-7));
    cmp
    (
        "  Inner-product with self = ",
        innerSqr(sT),
        SymmTensor2D<Type>
        (
            Type(5), Type(-4),
                     Type(13)
        )
    );
    cmp("  Square of Frobenius norm = ", magSqr(sT), scalar(17.999999999999996));
    cmp
    (
        "  Outer-product of a Vector2D with itself = ",
        sqr(v),
        SymmTensor2D<Type>
        (
            Type(9), Type(-6),
                     Type(4)
        )
    );
}


// Execute each global operator of SymmTensor2D<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor2D<Type> T
    (
        Type(1), Type(-2),
        Type(3), Type(-4)
    );
    const SymmTensor2D<Type> sT
    (
        Type(1), Type(2),
                 Type(-4)
    );
    const SphericalTensor2D<Type> spT(Type(-1));
    const Vector2D<Type> v(Type(3), Type(-2));
    const Type x(-4);

    Info<< "# Operands:" << nl
        << "  Tensor2D = " << T << nl
        << "  SymmTensor2D = " << sT << nl
        << "  SphericalTensor2D = " << spT << nl
        << "  Vector2D = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Sum of SpTensor2D-SymmTensor2D = ",
        (spT + sT),
        SymmTensor2D<Type>
        (
            Type(0), Type(2),
                     Type(-5)
        )
    );
    cmp
    (
        "  Sum of SymmTensor2D-SpTensor2D = ",
        (sT + spT),
        SymmTensor2D<Type>
        (
            Type(0), Type(2),
                     Type(-5)
        )
    );
    cmp
    (
        "  Subtract SymmTensor2D from SpTensor2D = ",
        (spT - sT),
        SymmTensor2D<Type>
        (
            Type(-2), Type(-2),
                      Type(3)
        )
    );
    cmp
    (
        "  Subtract SpTensor2D from SymmTensor2D = ",
        (sT - spT),
        SymmTensor2D<Type>
        (
            Type(2), Type(2),
                     Type(-3)
        )
    );
    cmp
    (
        "  Division of a SymmTensor2D by a Type",
        sT/x,
        SymmTensor2D<Type>
        (
            Type(-0.25), Type(-0.5),
                         Type(1)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor2D-SymmTensor2D = ",
        (sT & sT),
        Tensor2D<Type>
        (
            Type(5),  Type(-6),
            Type(-6), Type(20)
        )
    );
    cmp
    (
        "  Inner-product of SpTensor2D-SymmTensor2D = ",
        (spT & sT),
        SymmTensor2D<Type>
        (
            Type(-1), Type(-2),
                      Type(4)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor2D-SpTensor2D = ",
        (sT & spT),
        SymmTensor2D<Type>
        (
            Type(-1), Type(-2),
                      Type(4)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor2D-Vector2D = ",
        (sT & v),
        Vector2D<Type>(Type(-1), Type(14)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector2D-SymmTensor2D = ",
        (v & sT),
        Vector2D<Type>(Type(-1), Type(14)) // Row-vector
    );
    cmp("  D-inner-prod of SymmTensor2D-SymmTensor2D = ", (sT && sT), Type(25));
    cmp("  D-inner-prod of SymmTensor2D-SpTensor2D = ", (sT && spT), Type(3));
    cmp("  D-inner-prod of SpTensor2D-SymmTensor2D = ", (spT && sT), Type(3));
}


// Return false if given eigenvalues fail to satisy eigenvalue relations
// Relations: (Beauregard & Fraleigh (1973), ISBN 0-395-14017-X, p. 307)
void test_eigenvalues(const symmTensor2D& T, const vector2D& EVals)
{
    {
        const scalar determinant = det(T);
        const scalar EValsProd = EVals.x()*EVals.y();
        cmp("# Product of eigenvalues = det(T):", EValsProd, determinant, 1e-8);
    }

    {
        const scalar trace = tr(T);
        scalar EValsSum = 0.0;
        for (const auto& val : EVals)
        {
            EValsSum += val;
        }
        cmp("# Sum of eigenvalues = trace(T):", EValsSum, trace);
    }
}


// Return false if a given eigenvalue-eigenvector pair
// fails to satisfy the characteristic equation
void test_characteristic_equation
(
    const symmTensor2D& T,
    const vector2D& EVals,
    const tensor2D& EVecs
)
{
    Info<< "# Characteristic equation:" << nl;
    for (direction dir = 0; dir < pTraits<vector2D>::nComponents; ++dir)
    {
        Info<< "EVal = " << EVals[dir] << nl
            << "EVec = " << EVecs.row(dir) << nl;
        const vector2D leftSide(T & EVecs.row(dir));
        const vector2D rightSide(EVals[dir]*EVecs.row(dir));
        const vector2D X(leftSide - rightSide);

        for (const auto x : X)
        {
            cmp("  (sT & EVec - EVal*EVec) = 0:", x, 0.0, 1e-5);
        }
    }
}


// Return false if the eigen functions fail to satisfy relations
void test_eigen_funcs(const symmTensor2D& T)
{
    Info<< "# Operand:" << nl
        << "  symmTensor2D = " << T << nl;


    Info<< "# Return eigenvalues of a given symmTensor2D:" << nl;
    const vector2D EVals(eigenValues(T));
    Info<< EVals << endl;
    test_eigenvalues(T, EVals);

    Info<< "# Return eigenvectors of a given symmTensor2D corresponding to"
        << " given eigenvalues:" << nl;
    const tensor2D EVecs0(eigenVectors(T, EVals));
    Info<< EVecs0 << endl;
    test_characteristic_equation(T, EVals, EVecs0);

    Info<< "# Return eigenvectors of a given symmTensor2D by computing"
        << " the eigenvalues of the symmTensor2D in the background:" << nl;
    const tensor2D EVecs1(eigenVectors(T));
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
        "SymmTensor2D<floatScalar>",
        "SymmTensor2D<doubleScalar>",
        "SymmTensor2D<complex>"
    });

    run_tests(types, typeID);


    Info<< nl << "    ## Test symmTensor2D eigen functions: ##" << nl;
    const label numberOfTests = 10000;
    Random rndGen(1234);

    for (label i = 0; i < numberOfTests; ++i)
    {
        const symmTensor2D T(makeRandomContainer(rndGen));
        test_eigen_funcs(T);
    }

    {
        Info<< nl << "    ## A symmTensor2D with T.xy = VSMALL" << nl;
        for (label i = 0; i < numberOfTests; ++i)
        {
            symmTensor2D T(makeRandomContainer(rndGen));
            T.xy() = VSMALL;
            test_eigen_funcs(T);
        }
    }
    {
        Info<< nl << "    ## A symmTensor2D with T.xx = T.yy: ##" << nl;
        for (label i = 0; i < numberOfTests; ++i)
        {
            symmTensor2D T(makeRandomContainer(rndGen));
            T.xx() = T.yy();
            test_eigen_funcs(T);
        }
    }
    {
        Info<< nl << "    ## A zero symmTensor2D: ##"<< nl;
        const symmTensor2D zeroT(Zero);
        test_eigen_funcs(zeroT);
    }
    {
        Info<< nl << "    ## A stiff symmTensor2D: ##"<< nl;
        const symmTensor2D stiff
        (
            pow(10.0, 10), pow(10.0, -8),
                           pow(10.0, 9)
        );
        test_eigen_funcs(stiff);
    }
    {
        Info<< nl
            << "    ## Random symmTensor2D with tiny off-diag elements: ##"
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
                    symmTensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor2D T(makeRandomContainer(rndGen));
                    T.xy() = eps;
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


// ************************************************************************* //
