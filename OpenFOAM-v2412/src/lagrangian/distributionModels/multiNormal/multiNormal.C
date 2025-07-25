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

#include "multiNormal.H"
#include "mathematicalConstants.H"
#include "MathFunctions.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(multiNormal, 0);
    addToRunTimeSelectionTable(distributionModel, multiNormal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::multiNormal::multiNormal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    mu_
    (
        distributionModelDict_.lookupCompat
        (
            "mu",
            {{"expectation", 2112}}
        )
    ),
    sigma_
    (
        distributionModelDict_.lookupCompat
        (
            "sigma",
            {{"variance", 2112}}
        )
    ),
    weight_
    (
        distributionModelDict_.lookupCompat
        (
            "weight",
            {{"strength", 2112}}
        )
    )
{
    check();

    scalar sum = 0;
    for (label i = 0; i < weight_.size(); ++i)
    {
        if (i > 0 && weight_[i] < weight_[i-1])
        {
            FatalErrorInFunction
                << type() << "distribution: "
                << "Weights must be specified in a monotonic order." << nl
                << "Please see the row i = " << i << nl
                << "weight[i-1] = " << weight_[i-1] << nl
                << "weight[i] = " << weight_[i]
                << exit(FatalError);
        }

        sum += weight_[i];
    }

    if (sum < VSMALL)
    {
        FatalErrorInFunction
            << type() << "distribution: "
            << "The sum of weights cannot be zero." << nl
            << "weight = " << weight_
            << exit(FatalError);
    }

    for (label i = 1; i < weight_.size(); ++i)
    {
        weight_[i] += weight_[i-1];
    }

    for (auto& w : weight_)
    {
        w /= sum;
    }
}


Foam::distributionModels::multiNormal::multiNormal(const multiNormal& p)
:
    distributionModel(p),
    mu_(p.mu_),
    sigma_(p.sigma_),
    weight_(p.weight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::multiNormal::sample() const
{
    const scalar u = rndGen_.sample01<scalar>();

    for (label i = 0; i < weight_.size(); ++i)
    {
        if (weight_[i] > u)
        {
            return sample(mu_[i], sigma_[i]);
        }
    }

    const label last = weight_.size() - 1;

    return sample(mu_[last], sigma_[last]);
}


Foam::scalar Foam::distributionModels::multiNormal::sample
(
    const scalar mu,
    const scalar sigma
) const
{
    const scalar a = (minValue_ - mu)/sigma;
    const scalar b = (maxValue_ - mu)/sigma;

    const scalar aPhi = 0.5*(1.0 + erf(a/Foam::sqrt(2.0)));
    const scalar bPhi = 0.5*(1.0 + erf(b/Foam::sqrt(2.0)));

    const scalar u = rndGen_.sample01<scalar>();
    const scalar p = u*(bPhi - aPhi) + aPhi;

    // (B:p. 20-24)
    const scalar x =
        mu + sigma*Foam::sqrt(2.0)*Math::erfInv(2.0*p - 1.0);

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    return clamp(x, minValue_, maxValue_);
}


Foam::scalar Foam::distributionModels::multiNormal::meanValue() const
{
    scalar mean = 0;
    forAll(weight_, i)
    {
        mean += weight_[i]*mu_[i];
    }

    return mean;
}


// ************************************************************************* //
