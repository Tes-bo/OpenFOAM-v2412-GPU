/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "MULES.H"
#include "profiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::MULES::limitSum(UPtrList<scalarField>& phiPsiCorrs)
{
    const label nPhases = phiPsiCorrs.size();

    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumPos = 0;
        scalar sumNeg = 0;

        for (label phasei = 0; phasei < nPhases; ++phasei)
        {
            if (phiPsiCorrs[phasei][facei] > 0)
            {
                sumPos += phiPsiCorrs[phasei][facei];
            }
            else
            {
                sumNeg += phiPsiCorrs[phasei][facei];
            }
        }

        scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > VSMALL)
        {
            scalar lambda = -sumNeg/sumPos;

            for (label phasei = 0; phasei < nPhases; ++phasei)
            {
                if (phiPsiCorrs[phasei][facei] > 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -VSMALL)
        {
            scalar lambda = -sumPos/sumNeg;

            for (label phasei = 0; phasei < nPhases; ++phasei)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}


void Foam::MULES::limitSum
(
    const UPtrList<const scalarField>& alphas,
    UPtrList<scalarField>& phiPsiCorrs,
    const labelHashSet& fixed
)
{
    labelHashSet notFixed(identity(phiPsiCorrs.size()));
    notFixed -= fixed;

    forAll(phiPsiCorrs[0], facei)
    {
        scalar alphaNotFixed = 0, corrNotFixed = 0;
        for (const label phasei : notFixed)
        {
            alphaNotFixed += alphas[phasei][facei];
            corrNotFixed += phiPsiCorrs[phasei][facei];
        }

        scalar corrFixed = 0;
        for (const label phasei : fixed)
        {
            corrFixed += phiPsiCorrs[phasei][facei];
        }

        const scalar sumCorr = corrNotFixed + corrFixed;

        const scalar lambda = - sumCorr/alphaNotFixed;

        for (const label phasei : notFixed)
        {
            phiPsiCorrs[phasei][facei] += lambda*alphas[phasei][facei];
        }
    }
}

// ************************************************************************* //
