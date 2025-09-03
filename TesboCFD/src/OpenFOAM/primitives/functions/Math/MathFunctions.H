/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Namespace
    Foam::Math

Description
    A namespace for various mathematical functions.

    Reference:
    \verbatim
        Inverse error function (tag:W):
            Winitzki, S. (2008).
            A handy approximation for the error function and its inverse.
            A lecture note obtained through private communication.
            URL:https://sites.google.com/site/winitzki/sergei-winitzkis-files
            (Retrieved on: 16 Feb 2021).

        Incomplete gamma functions (tag:DM):
            DiDonato, A. R., & Morris Jr, A. H. (1986).
            Computation of the incomplete gamma
            function ratios and their inverse.
            ACM Transactions on Mathematical Software (TOMS), 12(4), 377-393.
            DOI:10.1145/22721.23109
    \endverbatim

Note
  - The algorithm in \c invIncGamma is described in (DM:Sec. 4).
  - The algorithm in \c incGammaRatio_Q is described in (DM:Sec. 3).
  - The accuracy parameter \c IND is set to a value of 1.

SourceFiles
    erfInv.C
    incGamma.C
    invIncGamma.C

\*---------------------------------------------------------------------------*/

#ifndef MathFunctions_H
#define MathFunctions_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace Math Declaration
\*---------------------------------------------------------------------------*/

namespace Math
{

//- Inverse error function of a real-number argument
//  \param y Real-number argument at which to evaluate. Domain: (-1, 1)
//  \return The inverse of error function of y
scalar erfInv(const scalar y);


// Incomplete gamma functions

//- Inverse of regularised lower incomplete gamma function
//  \param a Real-number argument. Domain: (0, infty]
//  \param P Real-number argument. Domain: [0,1]
scalar invIncGamma(const scalar a, const scalar P);

//- Regularised upper incomplete gamma function
//  \param a Real-number argument. Domain: (0, infty]
//  \param x Real-number argument. Domain: [0, infty]
scalar incGammaRatio_Q(const scalar a, const scalar x);

//- Regularised lower incomplete gamma function
//  \param a Real-number argument. Domain: (0, infty]
//  \param x Real-number argument. Domain: [0, infty]
scalar incGammaRatio_P(const scalar a, const scalar x);

//- Upper incomplete gamma function
//  \param a Real-number argument. Domain: (0, infty]
//  \param x Real-number argument. Domain: [0, infty]
scalar incGamma_Q(const scalar a, const scalar x);

//- Lower incomplete gamma function
//  \param a Real-number argument. Domain: (0, infty]
//  \param x Real-number argument. Domain: [0, infty]
scalar incGamma_P(const scalar a, const scalar x);


} // End namespace Math

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
