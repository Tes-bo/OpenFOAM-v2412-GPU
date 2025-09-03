/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tensor.H"
#include "cubicEqn.H"
#include "mathematicalConstants.H"
#include "SVD.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Vector<Foam::complex> Foam::eigenValues(const tensor& T)
{
    // Return diagonal if T is effectively diagonal tensor
    if
    (
        (
            sqr(T.xy()) + sqr(T.xz()) + sqr(T.yz())
          + sqr(T.yx()) + sqr(T.zx()) + sqr(T.zy())
        ) < ROOTSMALL
    )
    {
        return Vector<complex>
        (
            complex(T.xx()), complex(T.yy()), complex(T.zz())
        );
    }

    // Coefficients of the characteristic cubic polynomial (a = 1)
    const scalar b = - T.xx() - T.yy() - T.zz();
    const scalar c =
        T.xx()*T.yy() + T.xx()*T.zz() + T.yy()*T.zz()
      - T.xy()*T.yx() - T.yz()*T.zy() - T.zx()*T.xz();
    const scalar d =
      - T.xx()*T.yy()*T.zz()
      - T.xy()*T.yz()*T.zx() - T.xz()*T.zy()*T.yx()
      + T.xx()*T.yz()*T.zy() + T.yy()*T.zx()*T.xz() + T.zz()*T.xy()*T.yx();

    // Determine the roots of the characteristic cubic polynomial
    const Roots<3> roots(cubicEqn(1, b, c, d).roots());
    // Check the root types
    bool isComplex = false;
    forAll(roots, i)
    {
        switch (roots.type(i))
        {
            case roots::complex:
                isComplex = true;
                break;
            case roots::posInf:
            case roots::negInf:
            case roots::nan:
                WarningInFunction
                    << "Eigenvalue computation fails for tensor: " << T
                    << "due to the not-a-number root = " << roots[i]
                    << endl;
            case roots::real:
                break;
        }
    }

    if (isComplex)
    {
        return
            Vector<complex>
            (
                complex(roots[0], 0),
                complex(roots[1], roots[2]),
                complex(roots[1], -roots[2])
            );
    }

    return
        Vector<complex>
        (
            complex(roots[0], 0),
            complex(roots[1], 0),
            complex(roots[2], 0)
        );
}


Foam::Vector<Foam::complex> Foam::eigenVector
(
    const tensor& T,
    const complex& eVal,
    const Vector<complex>& standardBasis1,
    const Vector<complex>& standardBasis2
)
{
    // Construct the characteristic equation system for this eigenvalue
    Tensor<complex> A(Zero);
    forAll(A, i)
    {
        A[i] = complex(T[i], 0);
    }
    A.xx() -= eVal;
    A.yy() -= eVal;
    A.zz() -= eVal;

    // Determinants of the 2x2 sub-matrices used to find the eigenvectors
    // Sub-determinants for a unique eigenvenvalue
    complex sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    complex sd1 = A.zz()*A.xx() - A.zx()*A.xz();
    complex sd2 = A.xx()*A.yy() - A.xy()*A.yx();
    scalar magSd0 = mag(sd0);
    scalar magSd1 = mag(sd1);
    scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const Vector<complex> eVec
        (
            complex(1, 0),
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        const Vector<complex> eVec
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            complex(1, 0),
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd2 > SMALL)
    {
        const Vector<complex> eVec
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            complex(1, 0)
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }

    // Sub-determinants for a repeated eigenvalue
    sd0 = A.yy()*standardBasis1.z() - A.yz()*standardBasis1.y();
    sd1 = A.zz()*standardBasis1.x() - A.zx()*standardBasis1.z();
    sd2 = A.xx()*standardBasis1.y() - A.xy()*standardBasis1.x();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const Vector<complex> eVec
        (
            complex(1, 0),
            (A.yz()*standardBasis1.x() - standardBasis1.z()*A.yx())/sd0,
            (standardBasis1.y()*A.yx() - A.yy()*standardBasis1.x())/sd0
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        const Vector<complex> eVec
        (
            (standardBasis1.z()*A.zy() - A.zz()*standardBasis1.y())/sd1,
            complex(1, 0),
            (A.zx()*standardBasis1.y() - standardBasis1.x()*A.zy())/sd1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd2 > SMALL)
    {
        const Vector<complex> eVec
        (
            (A.xy()*standardBasis1.z() - standardBasis1.y()*A.xz())/sd2,
            (standardBasis1.x()*A.xz() - A.xx()*standardBasis1.z())/sd2,
            complex(1, 0)
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }

    // Triple eigenvalue
    return standardBasis1^standardBasis2;
}


Foam::Tensor<Foam::complex> Foam::eigenVectors
(
    const tensor& T,
    const Vector<complex>& eVals
)
{
    Vector<complex> Ux(complex(1, 0), Zero, Zero);
    Vector<complex> Uy(Zero, complex(1, 0), Zero);
    Vector<complex> Uz(Zero, Zero, complex(1, 0));

    Ux = eigenVector(T, eVals.x(), Uy, Uz);
    Uy = eigenVector(T, eVals.y(), Uz, Ux);
    Uz = eigenVector(T, eVals.z(), Ux, Uy);

    return Tensor<complex>(Ux, Uy, Uz);
}


Foam::Tensor<Foam::complex> Foam::eigenVectors(const tensor& T)
{
    const Vector<complex> eVals(eigenValues(T));

    return eigenVectors(T, eVals);
}


Foam::tensor Foam::pinv(const tensor& t)
{
    const scalar detval = t.det();

    if (detval < ROOTVSMALL)
    {
        // Fall back to pseudo inverse
        scalarRectangularMatrix mat(3, 3);

        mat(0,0) = t.xx(); mat(0,1) = t.xy(); mat(0,2) = t.xz();
        mat(1,0) = t.yx(); mat(1,1) = t.yy(); mat(1,2) = t.yz();
        mat(2,0) = t.zx(); mat(2,1) = t.zy(); mat(2,2) = t.zz();

        mat = SVDinv(mat);

        return tensor
        (
            mat(0,0), mat(0,1), mat(0,2),
            mat(1,0), mat(1,1), mat(1,2),
            mat(2,0), mat(2,1), mat(2,2)
        );
    }

    return Foam::inv(t, detval);
}


// ************************************************************************* //
