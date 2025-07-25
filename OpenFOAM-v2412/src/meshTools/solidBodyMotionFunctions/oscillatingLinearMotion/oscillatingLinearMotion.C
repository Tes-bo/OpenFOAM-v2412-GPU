/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "oscillatingLinearMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(oscillatingLinearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingLinearMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingLinearMotion::oscillatingLinearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    omegaPtr_(Function1<scalar>::New("omega", SBMFCoeffs_, &runTime)),
    phaseShiftPtr_
    (
        Function1<scalar>::NewIfPresent
        (
            "phaseShift",
            SBMFCoeffs_,
            word::null,
            &runTime
        )
    ),
    amplitudePtr_(Function1<vector>::New("amplitude", SBMFCoeffs_, &runTime)),
    verticalShiftPtr_
    (
        Function1<vector>::NewIfPresent
        (
            "verticalShift",
            SBMFCoeffs_,
            word::null,
            &runTime
        )
    )
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::oscillatingLinearMotion::transformation() const
{
    const scalar t = time_.value();

    const vector amplitude(amplitudePtr_->value(t));
    const scalar omega = omegaPtr_->value(t);

    scalar phaseShift = 0;
    if (phaseShiftPtr_)
    {
        phaseShift = phaseShiftPtr_->value(t);
    }

    vector verticalShift(Zero);
    if (verticalShiftPtr_)
    {
        verticalShift = verticalShiftPtr_->value(t);
    }

    const vector displacement
    (
        amplitude*sin(omega*(t + phaseShift)) + verticalShift
    );

    quaternion R(1);
    septernion TR(septernion(-displacement)*R);

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::oscillatingLinearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    if (!solidBodyMotionFunction::read(SBMFCoeffs))
    {
        return false;
    }

    omegaPtr_.reset
    (
        Function1<scalar>::New("omega", SBMFCoeffs_, &time_)
    );

    phaseShiftPtr_.reset
    (
        Function1<scalar>::NewIfPresent
        (
            "phaseShift",
            SBMFCoeffs_,
            word::null,
            &time_
        )
    );

    amplitudePtr_.reset
    (
        Function1<vector>::New("amplitude", SBMFCoeffs_, &time_)
    );

    verticalShiftPtr_.reset
    (
        Function1<vector>::NewIfPresent
        (
            "verticalShift",
            SBMFCoeffs_,
            word::null,
            &time_
        )
    );

    return true;
}


// ************************************************************************* //
