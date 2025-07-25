/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "multiMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(multiMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        multiMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiMotion::multiMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::multiMotion::transformation() const
{
    scalar t = time_.value();

    septernion TR = SBMFs_[0].transformation();

    for (label i = 1; i < SBMFs_.size(); i++)
    {
        TR *= SBMFs_[i].transformation();
    }

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::multiMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    label i = 0;
    SBMFs_.setSize(SBMFCoeffs_.size());

    for (const entry& dEntry : SBMFCoeffs_)
    {
        if (dEntry.isDict())
        {
            SBMFs_.set
            (
                i,
                solidBodyMotionFunction::New(dEntry.dict(), time_)
            );

            Info<< "Constructed SBMF " << i << " : "
                << dEntry.keyword() << " of type "
                << SBMFs_[i].type() << endl;

            i++;
        }
    }
    SBMFs_.setSize(i);

    return true;
}


// ************************************************************************* //
