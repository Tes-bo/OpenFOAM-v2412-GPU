/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tensor2D.H"
#include "quadraticEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::tensor2D::vsType::typeName = "tensor2D";

template<>
const char* const Foam::tensor2D::vsType::componentNames[] =
{
    "xx", "xy",
    "yx", "yy"
};

template<>
const Foam::tensor2D Foam::tensor2D::vsType::vsType::zero
(
    tensor2D::uniform(0)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::one
(
    tensor2D::uniform(1)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::max
(
    tensor2D::uniform(VGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::min
(
    tensor2D::uniform(-VGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::rootMax
(
    tensor2D::uniform(ROOTVGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::vsType::rootMin
(
    tensor2D::uniform(-ROOTVGREAT)
);

template<>
const Foam::tensor2D Foam::tensor2D::I
(
    1, 0,
    0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Vector2D<Foam::complex> Foam::eigenValues(const tensor2D& T)
{
    const scalar a = T.xx();
    const scalar b = T.xy();
    const scalar c = T.yx();
    const scalar d = T.yy();

    // Return diagonal if T is effectively diagonal tensor
    if ((sqr(b) + sqr(c)) < ROOTSMALL)
    {
        return Vector2D<complex>(complex(a), complex(d));
    }

    const scalar trace = a + d;

    // (JLM:p. 2246)
    scalar w = b*c;
    scalar e = std::fma(-b, c, w);
    scalar f = std::fma(a, d, -w);
    const scalar determinant = f + e;

    // Square-distance between two eigenvalues
    const scalar gapSqr = std::fma(-4.0, determinant, sqr(trace));

    // (F:Sec. 8.4.2.)
    // Eigenvalues are effectively real
    if (0 <= gapSqr)
    {
        scalar firstRoot = 0.5*(trace - sign(-trace)*Foam::sqrt(gapSqr));

        if (mag(firstRoot) < SMALL)
        {
            WarningInFunction
                << "Zero-valued root is found. Adding SMALL to the root "
                << "to avoid floating-point exception." << nl;
            firstRoot = SMALL;
        }

        Vector2D<complex> eVals
        (
            complex(firstRoot, 0),
            complex(determinant/firstRoot, 0)
        );

        // Sort the eigenvalues into ascending order
        if (eVals.x().real() > eVals.y().real())
        {
            std::swap(eVals.x(), eVals.y());
        }

        return eVals;
    }
    // Eigenvalues are complex
    else
    {
        const complex eVal(0.5*trace, 0.5*Foam::sqrt(mag(gapSqr)));

        return Vector2D<complex>
        (
            eVal,
            eVal.conjugate()
        );
    }
}


Foam::Vector2D<Foam::complex> Foam::eigenVector
(
    const tensor2D& T,
    const complex& eVal,
    const Vector2D<complex>& standardBasis
)
{
    // Construct the linear system for this eigenvalue
    const Tensor2D<complex> A
    (
        complex(T.xx()) - eVal,  complex(T.xy()),
        complex(T.yx()),         complex(T.yy()) - eVal
    );

    // Evaluate the eigenvector using the largest divisor
    if (mag(A.yy()) > mag(A.xx()) && mag(A.yy()) > SMALL)
    {
        Vector2D<complex> eVec(complex(1), -A.yx()/A.yy());

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
    else if (mag(A.xx()) > SMALL)
    {
        Vector2D<complex> eVec(-A.xy()/A.xx(), complex(1));

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
    // (K:p. 47-48)
    else if (mag(T.yx()) > mag(T.xy()) && mag(T.yx()) > SMALL)
    {
        const Vector2D<complex> eVec(eVal - T.yy(), complex(T.yx()));

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
    else if (mag(T.xy()) > SMALL)
    {
        const Vector2D<complex> eVec(complex(T.xy()), eVal - T.xx());

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

    // Repeated eigenvalue
    return Vector2D<complex>(-standardBasis.y(), standardBasis.x());
}


Foam::Tensor2D<Foam::complex> Foam::eigenVectors
(
    const tensor2D& T,
    const Vector2D<complex>& eVals
)
{
    Vector2D<complex> Ux(pTraits<complex>::one, Zero);
    Vector2D<complex> Uy(Zero, pTraits<complex>::one);

    Ux = eigenVector(T, eVals.x(), Uy);
    Uy = eigenVector(T, eVals.y(), Ux);

    return Tensor2D<complex>(Ux, Uy);
}


Foam::Tensor2D<Foam::complex> Foam::eigenVectors(const tensor2D& T)
{
    const Vector2D<complex> eVals(eigenValues(T));

    return eigenVectors(T, eVals);
}


// ************************************************************************* //
