/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "normal.H"
#include "mathematicalConstants.H"
#include "MathFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(normal, 0);
    addToRunTimeSelectionTable(distributionModel, normal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::normal::normal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    mu_
    (
        distributionModelDict_.getCompat<scalar>
        (
            "mu",
            {{"expectation", 2106}}
        )
    ),
    sigma_
    (
        distributionModelDict_.getCompat<scalar>
        (
            "sigma",
            {{"variance", 2106}}
        )
    )
{
    if (mag(sigma_) == 0)
    {
        FatalErrorInFunction
            << "Standard deviation cannot be zero." << nl
            << "    sigma = " << sigma_ << nl
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::normal::normal(const normal& p)
:
    distributionModel(p),
    mu_(p.mu_),
    sigma_(p.sigma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::normal::sample() const
{
    const scalar a = (minValue_ - mu_)/sigma_;
    const scalar b = (maxValue_ - mu_)/sigma_;

    const scalar aPhi = 0.5*(scalar(1) + erf(a/Foam::sqrt(scalar(2))));
    const scalar bPhi = 0.5*(scalar(1) + erf(b/Foam::sqrt(scalar(2))));

    const scalar u = rndGen_.sample01<scalar>();
    const scalar p = u*(bPhi - aPhi) + aPhi;

    // (B:p. 20-24)
    const scalar x =
        mu_ + sigma_*Foam::sqrt(scalar(2))*Math::erfInv(scalar(2)*p - scalar(1));

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    return clamp(x, minValue_, maxValue_);
}


Foam::scalar Foam::distributionModels::normal::meanValue() const
{
    const scalar a = (minValue_ - mu_)/sigma_;
    const scalar b = (maxValue_ - mu_)/sigma_;

    // (B:p. 2)
    const scalar aphi =
        scalar(1)/Foam::sqrt(scalar(2)*constant::mathematical::pi)
       *exp(-0.5*sqr(a));
    const scalar bphi =
        scalar(1)/Foam::sqrt(scalar(2)*constant::mathematical::pi)
       *exp(-0.5*sqr(b));

    // (B:p. 4)
    const scalar aPhi = 0.5*(scalar(1) + erf(a/Foam::sqrt(scalar(2))));
    const scalar bPhi = 0.5*(scalar(1) + erf(b/Foam::sqrt(scalar(2))));

    // (B:p. 25)
    return mu_ - sigma_*(bphi - aphi)/(bPhi - aPhi + VSMALL);
}


// ************************************************************************* //
