/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "symmTensor2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::symmTensor2D::vsType::typeName = "symmTensor2D";

template<>
const char* const Foam::symmTensor2D::vsType::componentNames[] =
{
    "xx", "xy",
          "yy"
};

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::vsType::zero
(
    symmTensor2D::uniform(0)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::one
(
    symmTensor2D::uniform(1)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::max
(
    symmTensor2D::uniform(VGREAT)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::min
(
    symmTensor2D::uniform(-VGREAT)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::rootMax
(
    symmTensor2D::uniform(ROOTVGREAT)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::vsType::rootMin
(
    symmTensor2D::uniform(-ROOTVGREAT)
);

template<>
const Foam::symmTensor2D Foam::symmTensor2D::I
(
    1, 0,
       1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector2D Foam::eigenValues(const symmTensor2D& T)
{
    // Return diagonal if T is effectively diagonal tensor
    if (sqr(T.xy()) < ROOTSMALL)
    {
        return vector2D(T.diag());
    }

    //(K:Eqs. 3.2-3.3)
    const scalar skewTrace = T.xx() - T.yy();
    const scalar trace = tr(T);
    const scalar gap = sign(skewTrace)*hypot(skewTrace, 2*T.xy());

    return vector2D(0.5*(trace + gap), 0.5*(trace - gap));
}


Foam::vector2D Foam::eigenVector
(
    const symmTensor2D& T,
    const scalar eVal,
    const vector2D& standardBasis
)
{
    // Construct the characteristic equation system for this eigenvalue
    const tensor2D A(T - eVal*tensor2D::I);

    // Evaluate the eigenvector using the largest divisor
    if (mag(A.yy()) > mag(A.xx()) && mag(A.yy()) > SMALL)
    {
        const vector2D eVec(1, -A.yx()/A.yy());

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
        const vector2D eVec(-A.xy()/A.xx(), 1);

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
    return vector2D(-standardBasis.y(), standardBasis.x());
}


Foam::tensor2D Foam::eigenVectors
(
    const symmTensor2D& T,
    const vector2D& eVals
)
{
    // (K:Eq. 3.5)
    const scalar skewTrace = T.xx() - T.yy();

    if (mag(skewTrace) > SMALL)
    {
        const scalar phi = 0.5*atan(2*T.xy()/skewTrace);
        const scalar cphi = cos(phi);
        const scalar sphi = sin(phi);
        return tensor2D(cphi, sphi, -sphi, cphi);
    }
    else if (mag(T.xy()) > SMALL)
    {
        const scalar a = 0.70710678;    // phi ~ 45deg
        return tensor2D(a, sign(T.xy())*a, -1*sign(T.xy())*a, a);
    }

    // (K:p. 3)
    return tensor2D(1, 0, 0, 1);
}


Foam::tensor2D Foam::eigenVectors(const symmTensor2D& T)
{
    const vector2D eVals(eigenValues(T));

    return eigenVectors(T, eVals);
}


// ************************************************************************* //
