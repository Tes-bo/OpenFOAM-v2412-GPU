/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionFunction, 0);
    defineRunTimeSelectionTable(solidBodyMotionFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunction::solidBodyMotionFunction
(
    const dictionary& dict,
    const Time& runTime
)
:
    SBMFCoeffs_
    (
        dict.found("solidBodyMotionFunction")
      ? dict.optionalSubDict
        (
            dict.get<word>("solidBodyMotionFunction") + "Coeffs"
        )
      : dict
    ),
    time_(runTime)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidBodyMotionFunction::read(const dictionary& dict)
{
    SBMFCoeffs_ = dict.optionalSubDict(type() + "Coeffs");

    return true;
}


void Foam::solidBodyMotionFunction::writeData(Ostream& os) const
{
    os << SBMFCoeffs_;
}


// ************************************************************************* //
