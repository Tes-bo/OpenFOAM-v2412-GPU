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
    Foam::linearUpwindWeight

Description
    Class with operator() which returns the weighting factors for the
    linear-upwind differencing scheme.  Note that the weighting factors are
    not bounded between upwind and central-differencing, some downwind
    contribution is possible although the interpolate is limited to be between
    the upwind and downwind cell values.

    Used in conjunction with the template class NVDscheme although this scheme
    is not NVD.

SourceFiles
    linearUpwindMake.C

\*---------------------------------------------------------------------------*/

#ifndef linearUpwind_H
#define linearUpwind_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearUpwindWeight Declaration
\*---------------------------------------------------------------------------*/

class linearUpwindWeight
{
public:

    // Generated Methods

        //- No copy construct
        linearUpwindWeight(const linearUpwindWeight&) = delete;

        //- No copy assignment
        void operator=(const linearUpwindWeight&) = delete;


    linearUpwindWeight(Istream&)
    {}

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
        scalar phif;
        if (faceFlux > 0)
        {
            phif = phiP + (1 - cdWeight)*(d & gradcP);
        }
        else
        {
            phif = phiN - cdWeight*(d & gradcN);
        }

        // Limit the estimated face value between the upwind and downwind cell
        // values
        phif = min(phif, max(phiN, phiP));
        phif = max(phif, min(phiN, phiP));

        return (phif - phiN)/stabilise(phiP - phiN, SMALL);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
