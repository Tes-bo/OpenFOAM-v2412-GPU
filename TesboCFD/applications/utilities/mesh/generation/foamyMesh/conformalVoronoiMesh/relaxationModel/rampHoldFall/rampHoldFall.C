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

#include "rampHoldFall.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rampHoldFall, 0);
    addToRunTimeSelectionTable(relaxationModel, rampHoldFall, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rampHoldFall::rampHoldFall
(
    const dictionary& relaxationDict,
    const Time& runTime
)
:
    relaxationModel(typeName, relaxationDict, runTime),
    rampStartRelaxation_(coeffDict().get<scalar>("rampStartRelaxation")),
    holdRelaxation_(coeffDict().get<scalar>("holdRelaxation")),
    fallEndRelaxation_(coeffDict().get<scalar>("fallEndRelaxation")),
    rampEndFraction_(coeffDict().get<scalar>("rampEndFraction")),
    fallStartFraction_(coeffDict().get<scalar>("fallStartFraction")),
    rampGradient_((holdRelaxation_ - rampStartRelaxation_)/(rampEndFraction_)),
    fallGradient_
    (
        (fallEndRelaxation_ - holdRelaxation_)/(1 - fallStartFraction_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::rampHoldFall::relaxation()
{
    scalar t = runTime_.time().timeOutputValue();

    scalar tStart = runTime_.time().startTime().value();
    scalar tEnd = runTime_.time().endTime().value();
    scalar tSpan = tEnd - tStart;

    if (tSpan < VSMALL)
    {
        return rampStartRelaxation_;
    }

    if (t - tStart < rampEndFraction_*tSpan)
    {
        // Ramp

        return rampGradient_*((t - tStart)/tSpan) + rampStartRelaxation_;
    }
    else if (t - tStart > fallStartFraction_*tSpan)
    {
        // Fall

        return
            fallGradient_*((t - tStart)/tSpan)
          + fallEndRelaxation_ - fallGradient_;
    }
    else
    {
        //Hold

        return holdRelaxation_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
