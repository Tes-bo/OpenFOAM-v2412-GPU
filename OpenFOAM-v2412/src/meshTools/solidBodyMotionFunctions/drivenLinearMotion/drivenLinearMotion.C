/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "drivenLinearMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(drivenLinearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        drivenLinearMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::drivenLinearMotion::drivenLinearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    cOfGdisplacement_(SBMFCoeffs.get<word>("cOfGdisplacement")),
    CofGdisp_
    (
        IOobject
        (
            cOfGdisplacement_,
            time_.timeName(),
            "uniform",
            time_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector(dimless, Zero)
    )
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::drivenLinearMotion::transformation() const
{

    DebugInFunction << "displacement  :" << CofGdisp_.value() << endl;
    quaternion R(1);
    septernion TR(septernion(-CofGdisp_.value())*R);

    DebugInFunction << "Time = " << time_.value()
                    << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::drivenLinearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
