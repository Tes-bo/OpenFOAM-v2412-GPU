/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "piecewiseLinearRamp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(piecewiseLinearRamp, 0);
    addToRunTimeSelectionTable
    (
        faceAreaWeightModel,
        piecewiseLinearRamp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::piecewiseLinearRamp::piecewiseLinearRamp
(
    const dictionary& faceAreaWeightDict
)
:
    faceAreaWeightModel(typeName, faceAreaWeightDict),
    lAF_(coeffDict().get<scalar>("lowerAreaFraction")),
    uAF_(coeffDict().get<scalar>("upperAreaFraction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::piecewiseLinearRamp::faceAreaWeight(scalar faceAreaFraction) const
{
    if (faceAreaFraction < lAF_)
    {
        return 0;
    }
    else if (faceAreaFraction < uAF_)
    {
        return faceAreaFraction/((uAF_ - lAF_)) - lAF_/(uAF_ - lAF_);
    }

    return 1;
}


// ************************************************************************* //
