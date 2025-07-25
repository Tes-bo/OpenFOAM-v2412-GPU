/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Class
    Foam::GammaWeight

Description
    Class with operator() which returns the weighting factors for the
    Gamma differencing scheme.  Used in conjunction with the template class
    NVDscheme.

SourceFiles
    GammaMake.C

\*---------------------------------------------------------------------------*/

#ifndef Gamma_H
#define Gamma_H

#include "scalar.H"
#include "vector.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GammaWeight Declaration
\*---------------------------------------------------------------------------*/

class GammaWeight
{
    // Private Data

        //- Model coefficient [0,1]
        scalar k_;


public:

    GammaWeight(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorInFunction(is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Rescale k_ to be >= 0 and <= 0.5 (TVD conformant)
        // and avoid the /0 when k_ = 0
        k_ = max(0.5*k_, SMALL);
    }


    scalar weight
    (
        scalar cdWeight,
        scalar faceFlux,
        scalar phiP,
        scalar phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d
    ) const
    {
        const vector dHat(normalised(d));

        // Choose gradc based on faceFlux
        const vector& gradcPN = (faceFlux > 0) ? gradcP : gradcN;
        const scalar udWeight = (faceFlux > 0) ? 1 : 0;

        // Stabilise for division
        const scalar gradcf = stabilise(dHat & gradcPN, SMALL);

        const scalar gradf = (phiN - phiP)/mag(d);

        const scalar phict = 1 - 0.5*gradf/gradcf;
        const scalar limiter = clamp(phict/k_, zero_one{});

        return lerp(udWeight, cdWeight, limiter);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
