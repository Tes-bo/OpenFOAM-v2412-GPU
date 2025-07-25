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
    Test-SymmTensor

Description
    Tests for \c SymmTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c symmTensor, i.e. SymmTensor<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "symmTensor.H"
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


// Create a random symmTensor
symmTensor makeRandomContainer(Random& rnd)
{
    symmTensor T(Zero);
    std::generate(T.begin(), T.end(), [&]{ return rnd.GaussNormal<scalar>(); });
    return T;
}


// Create a symmTensor based on a given value
template<class Type>
typename std::enable_if
<
    std::is_same<floatScalar, Type>::value ||
    std::is_same<doubleScalar, Type>::value,
    symmTensor
>::type makeContainer(const Type val)
{
    symmTensor T(Zero);
    std::fill(T.begin(), T.end(), val);
    return T;
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
    const scalar absTol = 0,    //<! useful for cmps near zero
    const scalar relTol = 1e-8  //<! are values the same within 8 decimals
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


// Create each constructor of SymmTensor<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const SymmTensor<Type> sT(Zero);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given VectorSpace of the same rank:" << nl;
        const VectorSpace<SymmTensor<Type>, Type, 6> M(Zero);
        const SymmTensor<Type> sT(M);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given SphericalTensor:" << nl;
        const SphericalTensor<Type> Sp(Type(5));
        const SymmTensor<Type> sT(Sp);
        Info<< sT << endl;
    }
    {
        Info<< "# Construct given the six components:" << nl;
        const SymmTensor<Type> sT
        (
            Type(1), Type(2), Type(-3),
                     Type(5), Type(-6),
                              Type(-9)
        );
        Info<< sT << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const SymmTensor<Type> sT(Zero);
        const SymmTensor<Type> copysT(sT);
        Info<< sT << tab << copysT << endl;
    }
}


// Execute each member function of SymmTensor<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SymmTensor<Type> sT
    (
        Type(1), Type(2), Type(-3),
                 Type(5), Type(-6),
                          Type(-9)
    );
    const SymmTensor<Type> csT
    (
        Type(1), Type(2), Type(-3),
                 Type(5), Type(-6),
                          Type(-9)
    );

    Info<< "# Operand: " << nl
        << "  SymmTensor = " << sT << endl;


    {
        Info<< "# Component access:" << nl;

        SymmTensor<Type> cpsT
        (
            sT.xx(), sT.xy(), sT.xz(),
                     sT.yy(), sT.yz(),
                              sT.zz()
        );
        cmp("  'SymmTensor' access:", sT, cpsT);
        cmp("  xy()=yx():", sT.xy(), sT.yx());
        cmp("  xz()=zx():", sT.xz(), sT.zx());
        cmp("  yz()=zy():", sT.yz(), sT.zy());

        const SymmTensor<Type> cpcsT
        (
            csT.xx(), csT.xy(), csT.xz(),
                      csT.yy(), csT.yz(),
                                csT.zz()
        );
        cmp("  'const SymmTensor' access:", csT, cpcsT);
        cmp("  xy()=yx():", sT.xy(), sT.yx());
        cmp("  xz()=zx():", sT.xz(), sT.zx());
        cmp("  yz()=zy():", sT.yz(), sT.zy());
    }
    {
        Info<< "# Diagonal access:" << nl;

        cmp
        (
            "  'SymmTensor'.diag():",
            sT.diag(),
            Vector<Type>(Type(1), Type(5), Type(-9))
        );
        cmp
        (
            "  'const SymmTensor'.diag():",
            csT.diag(),
            Vector<Type>(Type(1), Type(5), Type(-9))
        );


        Info<< "# Diagonal manipulation:" << nl;

        sT.diag(Vector<Type>(Type(-10), Type(-15), Type(-20)));
        cmp
        (
            "  'SymmTensor'.diag('Vector'):",
            sT.diag(),
            Vector<Type>(Type(-10), Type(-15), Type(-20))
        );
    }
    {
        Info<< "# Tensor operations:" << nl;

        Info<< "  Transpose:" << nl;
        cmp("  'SymmTensor'.T():", sT.T(), sT);
    }
    {
        Info<< "# Member operators:" << nl;

        sT = SphericalTensor<Type>(Type(5));
        cmp
        (
            "  Assign to a SphericalTensor:",
            sT,
            SymmTensor<Type>
            (
                Type(5),  Zero,    Zero,
                          Type(5), Zero,
                                   Type(5)
            )
        );
    }
}


// Execute each global function of SymmTensor<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const SymmTensor<Type> sT
    (
        Type(1), Type(2), Type(-3),
                 Type(5), Type(-6),
                          Type(-9)
    );

    Info<< "# Operand: " << nl
        << "  SymmTensor = " << sT << nl << endl;


    cmp("  Trace = ", tr(sT), Type(-3));
    cmp("  Spherical part = ", sph(sT), SphericalTensor<Type>(tr(sT)/Type(3)));
    cmp("  Symmetric part = ", symm(sT), sT);
    cmp("  Twice the symmetric part = ", twoSymm(sT), 2*sT);
    cmp
    (
        "  Deviatoric part = ",
        dev(sT),
        SymmTensor<Type>
        (
            Type(2), Type(2), Type(-3),
                     Type(6), Type(-6),
                              Type(-8)
        )
    );
    cmp("  Two-third deviatoric part = ", dev2(sT), sT - 2*sph(sT));
    cmp("  Determinant = ", det(sT), Type(-17.999999999999996));
    cmp
    (
        "  Cofactor tensor = ",
        cof(sT),
        SymmTensor<Type>
        (
            Type(-81), Type(36),  Type(3),
                       Type(-18), Type(0),
                                  Type(1)
        )
    );
    cmp
    (
        "  Inverse = ",
        inv(sT, det(sT)),
        SymmTensor<Type>
        (
            Type(4.5), Type(-2), Type(-0.16666667),
                       Type(1),  Type(0),
                                 Type(-0.05555556)
        ),
        1e-8
    );
    cmp
    (
        "  Inverse (another) = ",
        inv(sT),
        SymmTensor<Type>
        (
            Type(4.5), Type(-2), Type(-0.16666667),
                       Type(1),  Type(0),
                                 Type(-0.05555556)
        ),
        1e-8
    );
    cmp("  First invariant = ", invariantI(sT), Type(-3));
    cmp("  Second invariant = ", invariantII(sT), Type(-98));
    cmp("  Third invariant = ", invariantIII(sT), Type(-17.999999999999996));
    cmp
    (
        "  Inner-product with self = ",
        innerSqr(sT),
        SymmTensor<Type>
        (
            Type(14), Type(30), Type(12),
                      Type(65), Type(18),
                                Type(126)
        )
    );
    cmp("  Square of Frobenius norm = ", magSqr(sT), scalar(205));
}


// Execute each global operator of SymmTensor<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor<Type> T
    (
        Type(1), Type(2), Type(-3),
        Type(4), Type(5), Type(-6),
        Type(7), Type(8), Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(1), Type(2), Type(-3),
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
        "  Sum of SpTensor-SymmTensor = ",
        (spT + sT),
        SymmTensor<Type>
        (
            Type(2), Type(2), Type(-3),
                     Type(6), Type(-6),
                              Type(-8)
        )
    );
    cmp
    (
        "  Sum of SymmTensor-SpTensor = ",
        (sT + spT),
        SymmTensor<Type>
        (
            Type(2), Type(2), Type(-3),
                     Type(6), Type(-6),
                              Type(-8)
        )
    );
    cmp
    (
        "  Subtract SymmTensor from SpTensor = ",
        (spT - sT),
        SymmTensor<Type>
        (
            Type(0), Type(-2), Type(3),
                     Type(-4), Type(6),
                               Type(10)
        )
    );
    cmp
    (
        "  Subtract SpTensor from SymmTensor = ",
        (sT - spT),
        SymmTensor<Type>
        (
            Type(0), Type(2), Type(-3),
                     Type(4), Type(-6),
                              Type(-10)
        )
    );
    cmp
    (
        "  Hodge dual of a SymmTensor",
        *sT,
        Vector<Type>(Type(-6), Type(3), Type(2))
    );
    cmp
    (
        "  Division of a SymmTensor by a Type",
        sT/x,
        SymmTensor<Type>
        (
            Type(0.25), Type(0.5),  Type(-0.75),
                        Type(1.25), Type(-1.5),
                                    Type(-2.25)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor-SymmTensor = ",
        (sT & sT),
        Tensor<Type>
        (
            Type(14), Type(30), Type(12),
            Type(30), Type(65), Type(18),
            Type(12), Type(18), Type(126)
        )
    );
    cmp
    (
        "  Inner-product of SpTensor-SymmTensor = ",
        (spT & sT),
        SymmTensor<Type>
        (
            Type(1), Type(2), Type(-3),
                     Type(5), Type(-6),
                              Type(-9)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor-SpTensor = ",
        (sT & spT),
        SymmTensor<Type>
        (
            Type(1), Type(2), Type(-3),
                     Type(5), Type(-6),
                              Type(-9)
        )
    );
    cmp
    (
        "  Inner-product of SymmTensor-Vector = ",
        (sT & v),
        Vector<Type>(Type(4), Type(10), Type(-30)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector-SymmTensor = ",
        (v & sT),
        Vector<Type>(Type(4), Type(10), Type(-30)) // Row-vector
    );
    cmp("  D-inner-product of SymmTensor-SymmTensor = ", (sT && sT), Type(205));
    cmp("  D-inner-product of SymmTensor-SpTensor = ", (sT && spT), Type(-3));
    cmp("  D-inner-product of SpTensor-SymmTensor = ", (spT && sT), Type(-3));
}


// Return false if given eigenvalues fail to satisy eigenvalue relations
// Relations: (Beauregard & Fraleigh (1973), ISBN 0-395-14017-X, p. 307)
void test_eigenvalues(const symmTensor& T, const vector& EVals)
{
    {
        const scalar determinant = det(T);
        const scalar EValsProd = EVals.x()*EVals.y()*EVals.z();
        cmp("# Product of eigenvalues = det(T):", EValsProd, determinant, 1e-6);
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
    const symmTensor& T,
    const vector& EVals,
    const tensor& EVecs
)
{
    Info<< "# Characteristic equation:" << nl;
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        Info<< "EVal = " << EVals[dir] << nl
            << "EVec = " << EVecs.row(dir) << endl;
        const vector leftSide(T & EVecs.row(dir));
        const vector rightSide(EVals[dir]*EVecs.row(dir));
        const vector X(leftSide - rightSide);

        for (const auto x : X)
        {
            cmp("  (T & EVec - EVal*EVec) = 0:", x, 0.0, 1e-5);
        }
    }
}


// Return false if the eigen functions fail to satisfy relations
void test_eigen_funcs(const symmTensor& T)
{
    Info<< "# Operand:" << nl
        << "  symmTensor = " << T << nl;


    Info<< "# Return eigenvalues of a given symmTensor:" << nl;
    const vector EVals(eigenValues(T));
    Info<< EVals << endl;
    test_eigenvalues(T, EVals);

    Info<< "# Return eigenvectors of a given symmTensor corresponding to"
        << " given eigenvalues:" << nl;
    const tensor EVecs0(eigenVectors(T, EVals));
    Info<< EVecs0 << endl;
    test_characteristic_equation(T, EVals, EVecs0);

    Info<< "# Return eigenvectors of a given symmTensor by computing"
        << " the eigenvalues of the symmTensor in the background:" << nl;
    const tensor EVecs1(eigenVectors(T));
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
        "SymmTensor<floatScalar>",
        "SymmTensor<doubleScalar>",
        "SymmTensor<complex>"
    });

    run_tests(types, typeID);


    Info<< nl << "    ## Test symmTensor eigen functions: ##" << nl;
    const label numberOfTests = 10000;
    Random rndGen(1234);

    for (label i = 0; i < numberOfTests; ++i)
    {
        const symmTensor T(makeRandomContainer(rndGen));
        test_eigen_funcs(T);
    }

    {
        Info<< nl << "    ## A zero symmTensor: ##"<< nl;
        const symmTensor zeroT(Zero);
        test_eigen_funcs(zeroT);
    }
    {
        Info<< nl
            << "    ## A symmTensor with 2 repeated eigenvalues: ##"
            << nl;
        const symmTensor T
        (
            1.0, 0.0, Foam::sqrt(2.0),
                 2.0, 0.0,
                      0.0
        );
        test_eigen_funcs(T);
    }
    {
        Info<< nl
            << "    ## A symmTensor with 3 repeated eigenvalues: ##"
            << nl;
        const symmTensor T
        (
            0.023215, -5.0739e-09, -7.0012e-09,
                         0.023215, -8.148e-10,
                                    0.023215
        );
        test_eigen_funcs(T);
    }
    {
        Info<< nl << "    ## A stiff symmTensor: ##" << nl;
        const symmTensor stiff
        (
            pow(10.0, 10), pow(10.0, 8),  pow(10.0, -8),
                           pow(10.0, -8), pow(10.0, 8),
                                          pow(10.0, 7)
        );
        test_eigen_funcs(stiff);
    }
    {
        Info<< nl
            << "    ## Random symmTensors with tiny off-diag elements: ##"
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
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps*rndGen.GaussNormal<scalar>();
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = eps*rndGen.GaussNormal<scalar>();
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps;
                    T.xz() = eps;
                    T.yz() = eps;
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = eps;
                    T.xz() = eps;
                    T.yz() = eps;
                    T.zz() = eps;
                    test_eigen_funcs(T);
                }
                {
                    symmTensor T(makeRandomContainer(rndGen));
                    T.xy() = 0;
                    T.xz() = eps*rndGen.GaussNormal<scalar>();
                    T.yz() = 0;
                    test_eigen_funcs(T);
                }
            }
        }
    }
    #if 0
    // Numerical diagonalisation of 2x2 or 3x3 matrices with analytic methods
    // are, like the methods currently being used in OpenFOAM, inherently error
    // prone. Despite its speed, the analytic methods may becomes inaccurate or
    // may even fail completely if the matrix entries differ greatly in
    // magnitude, particularly with large off-diagonal elements.
    // The remedy is to use iterative or hybrid analytic/iterative methods
    // such as published here (for 3x3/2x2 matrices):
    // (Kopp, 2008) arXiv.org: physics/0610206
    // mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
    {
        Info<< nl << "    ## symmTensors consisting machine epsilons: ##" << nl;
        Info<< "    # floatScalar" << nl;
        const List<floatScalar> floatEpsilons
        ({
            floatScalarGREAT, floatScalarVGREAT, floatScalarROOTVGREAT,
            floatScalarSMALL, floatScalarVSMALL, floatScalarROOTVSMALL,
            Foam::sqrt(floatScalarSMALL), 0
        });

        for (const auto& eps : floatEpsilons)
        {
            const symmTensor T(makeContainer(eps));
            test_eigen_funcs(T);
        }

        Info<< "    # doubleScalar" << nl;
        const List<doubleScalar> doubleEpsilons
        ({
            doubleScalarGREAT, doubleScalarROOTVGREAT, // doubleVGREAT fails
            doubleScalarSMALL, doubleScalarVSMALL, doubleScalarROOTVSMALL,
            Foam::sqrt(doubleScalarSMALL), 0
        });

        for (const auto& eps : doubleEpsilons)
        {
            const symmTensor T(makeContainer(eps));
            test_eigen_funcs(T);
        }
    }
    {
        Info<< nl
            << "    ## Random symmTensors with machine eps off-diag elmes: ##"
            << nl;

        const List<floatScalar> floatEpsilons
        ({
            floatScalarGREAT, floatScalarVGREAT, floatScalarROOTVGREAT,
            floatScalarSMALL, floatScalarVSMALL, floatScalarROOTVSMALL
        });

        const List<doubleScalar> doubleEpsilons
        ({
            doubleScalarGREAT, doubleScalarVGREAT, doubleScalarROOTVGREAT,
            doubleScalarSMALL, doubleScalarVSMALL, doubleScalarROOTVSMALL
        });

        for (label i = 0; i < numberOfTests; ++i)
        {
            symmTensor T(makeRandomContainer(rndGen));

            for (const auto& eps : floatEpsilons)
            {
                T.xy() = eps;
                T.xz() = eps;
                T.yz() = eps;
                test_eigen_funcs(T);
            }

            for (const auto& eps : doubleEpsilons)
            {
                T.xy() = eps;
                T.xz() = eps;
                T.yz() = eps;
                test_eigen_funcs(T);
            }
        }
    }
    #endif


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
