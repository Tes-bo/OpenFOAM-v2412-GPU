/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "RosinRammler.H"
#include "MathFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(RosinRammler, 0);
    addToRunTimeSelectionTable(distributionModel, RosinRammler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::RosinRammler::RosinRammler
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    lambda_(distributionModelDict_.getCompat<scalar>("lambda", {{"d", 2106}})),
    n_(distributionModelDict_.get<scalar>("n"))
{
    const word parcelBasisType =
        dict.getOrDefault<word>("parcelBasisType", "none");

    if (parcelBasisType == "mass")
    {
        WarningInFunction
            << "Selected parcel basis type: " << parcelBasisType << nl
            << "    Please consider to use massRosinRammler distribution model"
            << endl;
    }

    if (lambda_ < VSMALL || n_ < VSMALL)
    {
        FatalErrorInFunction
            << "Scale/Shape parameter cannot be equal to or less than zero:"
            << "    lambda = " << lambda_
            << "    n = " << n_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::RosinRammler::RosinRammler(const RosinRammler& p)
:
    distributionModel(p),
    lambda_(p.lambda_),
    n_(p.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::RosinRammler::sample() const
{
    const scalar u = rndGen_.sample01<scalar>();
    const scalar qMin = pow(minValue_/lambda_, n_);
    const scalar qMax = pow(maxValue_/lambda_, n_);
    const scalar r = scalar(1) - exp(-qMax + qMin);
    return lambda_*pow(qMin - log(scalar(1) - u*r), scalar(1)/n_);
}


Foam::scalar Foam::distributionModels::RosinRammler::meanValue() const
{
    // (C:Eq. 5)
    const scalar a = scalar(1)/lambda_ + scalar(1);
    const scalar qMax = pow(maxValue_/n_, lambda_);
    const scalar qMin = pow(minValue_/n_, lambda_);
    const scalar gMax = Math::incGamma_P(a, qMax);
    const scalar gMin = Math::incGamma_P(a, qMin);

    return n_/(exp(-qMin) - exp(-qMax))*(gMax - gMin);
}


// ************************************************************************* //
