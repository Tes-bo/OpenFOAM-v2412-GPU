/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "symmTensor.H"
#include "cubicEqn.H"
#include "mathematicalConstants.H"
#include "SVD.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::symmTensor::vsType::typeName = "symmTensor";

template<>
const char* const Foam::symmTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
          "yy", "yz",
                "zz"
};

template<>
const Foam::symmTensor Foam::symmTensor::vsType::vsType::zero
(
    symmTensor::uniform(0)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::one
(
    symmTensor::uniform(1)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::max
(
    symmTensor::uniform(VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::min
(
    symmTensor::uniform(-VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMax
(
    symmTensor::uniform(ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMin
(
    symmTensor::uniform(-ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::I
(
    1, 0, 0,
       1, 0,
          1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const symmTensor& T)
{
    // Return ascending diagonal if T is effectively diagonal tensor
    if ((sqr(T.xy()) + sqr(T.xz()) + sqr(T.yz())) < ROOTSMALL)
    {
        vector eVals(T.diag());

        std::sort(eVals.begin(), eVals.end());

        return eVals;
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
    Roots<3> roots(cubicEqn(1, b, c, d).roots());

    vector eVals(Zero);

    // Check the root types
    forAll(roots, i)
    {
        switch (roots.type(i))
        {
            case roots::real:
                eVals[i] = roots[i];
                break;
            case roots::complex:
            case roots::posInf:
            case roots::negInf:
            case roots::nan:
                WarningInFunction
                    << "Eigenvalue computation fails for symmTensor: " << T
                    << "due to the non-real root = " << roots[i]
                    << endl;
                eVals[i] = roots[i];
                break;
        }
    }

    // Sort the eigenvalues into ascending order
    std::sort(eVals.begin(), eVals.end());

    return eVals;
}


Foam::vector Foam::eigenVector
(
    const symmTensor& T,
    const scalar eVal,
    const vector& standardBasis1,
    const vector& standardBasis2
)
{
    // Construct the characteristic equation system for this eigenvalue
    const tensor A(T - eVal*I);

    {
        // Determinants of the 2x2 sub-matrices used to find the eigenvectors
        // Sub-determinants for a unique eigenvenvalue
        const scalar sd0 = A.yy()*A.zz() - A.yz()*A.zy();
        const scalar sd1 = A.zz()*A.xx() - A.zx()*A.xz();
        const scalar sd2 = A.xx()*A.yy() - A.xy()*A.yx();
        const scalar magSd0 = mag(sd0);
        const scalar magSd1 = mag(sd1);
        const scalar magSd2 = mag(sd2);

        // Evaluate the eigenvector using the largest sub-determinant
        if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
        {
            const vector eVec
            (
                1,
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
            const vector eVec
            (
                (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
                1,
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
            const vector eVec
            (
                (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
                (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
                1
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
    }

    // Sub-determinants for a repeated eigenvalue
    const scalar sd0 = A.yy()*standardBasis1.z() - A.yz()*standardBasis1.y();
    const scalar sd1 = A.zz()*standardBasis1.x() - A.zx()*standardBasis1.z();
    const scalar sd2 = A.xx()*standardBasis1.y() - A.xy()*standardBasis1.x();
    const scalar magSd0 = mag(sd0);
    const scalar magSd1 = mag(sd1);
    const scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const vector eVec
        (
            1,
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
        const vector eVec
        (
            (standardBasis1.z()*A.zy() - A.zz()*standardBasis1.y())/sd1,
            1,
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
        const vector eVec
        (
            (A.xy()*standardBasis1.z() - standardBasis1.y()*A.xz())/sd2,
            (standardBasis1.x()*A.xz() - A.xx()*standardBasis1.z())/sd2,
            1
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


Foam::tensor Foam::eigenVectors
(
    const symmTensor& T,
    const vector& eVals
)
{
    vector Ux(1, 0, 0), Uy(0, 1, 0), Uz(0, 0, 1);

    Ux = eigenVector(T, eVals.x(), Uy, Uz);
    Uy = eigenVector(T, eVals.y(), Uz, Ux);
    Uz = eigenVector(T, eVals.z(), Ux, Uy);

    return tensor(Ux, Uy, Uz);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T)
{
    const vector eVals(eigenValues(T));

    return eigenVectors(T, eVals);
}


Foam::symmTensor Foam::pinv(const symmTensor& st)
{
    const scalar detval = st.det();

    if (detval < ROOTVSMALL)
    {
        // Fall back to pseudo inverse
        scalarRectangularMatrix mat(3, 3);

        mat(0,0) = st.xx(); mat(0,1) = st.xy(); mat(0,2) = st.xz();
        mat(1,0) = st.yx(); mat(1,1) = st.yy(); mat(1,2) = st.yz();
        mat(2,0) = st.zx(); mat(2,1) = st.zy(); mat(2,2) = st.zz();

        mat = SVDinv(mat);

        return symmTensor
        (
            mat(0,0), mat(0,1), mat(0,2),
                      mat(1,1), mat(1,2),
                                mat(2,2)
        );
    }

    return Foam::inv(st, detval);
}


// ************************************************************************* //
