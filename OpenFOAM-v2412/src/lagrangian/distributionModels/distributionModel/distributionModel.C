/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "distributionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributionModel, 0);
    defineRunTimeSelectionTable(distributionModel, dictionary);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::distributionModel::check() const
{
    if (minValue() < 0)
    {
        FatalErrorInFunction
            << type() << "Distribution: "
            << "Minimum value must be greater than zero." << nl
            << "Supplied minValue = " << minValue()
            << abort(FatalError);
    }

    if (maxValue() < minValue())
    {
        FatalErrorInFunction
            << type() << "Distribution: "
            << "Maximum value cannot be smaller than minimum value" << nl
            << "    maxValue = " << maxValue()
            << ", minValue = " << minValue()
            << abort(FatalError);
    }

    if (maxValue() == minValue())
    {
        WarningInFunction
            << type() << "Distribution: "
            << "Maximum and minimum values are equal to each other" << nl
            << "    maxValue = " << maxValue()
            << ", minValue = " << minValue()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModel::distributionModel
(
    const word& name,
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModelDict_(dict),
    rndGen_(rndGen),
    minValue_(distributionModelDict_.getOrDefault<scalar>("minValue", GREAT)),
    maxValue_(distributionModelDict_.getOrDefault<scalar>("maxValue", -GREAT))
{}


Foam::distributionModel::distributionModel
(
    const distributionModel& p
)
:
    distributionModelDict_(p.distributionModelDict_),
    rndGen_(p.rndGen_),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModel::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModel::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
