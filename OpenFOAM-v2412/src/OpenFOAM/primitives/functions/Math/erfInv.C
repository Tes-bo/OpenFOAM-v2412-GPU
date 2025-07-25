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

\*---------------------------------------------------------------------------*/

#include "MathFunctions.H"
#include "mathematicalConstants.H"
#include "error.H"
#include <cmath>
#include <limits>

// * * * * * * * * * * * * * * * Global Functions * * * * * * * * * * * * * * //

Foam::scalar Foam::Math::erfInv(const scalar y)
{
    #ifdef FULLDEBUG
    if (mag(y) >= scalar(1))
    {
        WarningInFunction
            << "The domain of inverse error function argument "
            << "(i.e. y) should be limited to (-1, 1):" << nl
            << "    y = " << y
            << endl;

        return std::numeric_limits<scalar>::infinity();
    }
    #endif

    // (W:p. 2) to reduce the max relative error to O(1e-4)
    constexpr scalar a = 0.147;

    const scalar k =
        scalar(2)/(a*constant::mathematical::pi) + 0.5*log(scalar(1) - sqr(y));

    const scalar h = log(scalar(1) - sqr(y))/a;

    // (W:Eq. 7)
    const scalar x = sqrt(-k + sqrt(sqr(k) - h));

    if (y < 0)
    {
        return -x;
    }

    return x;
}


// ************************************************************************* //
