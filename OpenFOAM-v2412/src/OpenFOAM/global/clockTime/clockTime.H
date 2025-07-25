/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

Class
    Foam::clockTime

Description
    Starts timing and returns elapsed time from start.
    Uses std::chrono::high_resolution_clock for better resolution
    (2uSec instead of ~20mSec) than cpuTime.

Note
    It has twice the storage requirement of a simple clockValue since
    it tracks both total and incremental elapsed times.
    Additionally, it always invokes a clock query on construction
    which may make it less desirable for arrays of values (for example).

See Also
    Foam::clockValue

SourceFiles
    clockTimeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_clockTime_H
#define Foam_clockTime_H

#include "clockValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class clockTime Declaration
\*---------------------------------------------------------------------------*/

class clockTime
{
    // Private Data

        //- Time point at start, or after resetTime
        clockValue start_;

        //- Time point when elapsedTime or timeIncrement was called.
        //- Also updated by resetTime and resetTimeIncrement.
        mutable clockValue last_;


public:

    // Constructors

        //- Construct with the current clock value for the start point
        inline clockTime();

        //- Implicit construct from the clock value as the start point
        inline clockTime(const clockValue& clockval);


    // Member Functions

        //- Reset to use the current clock value for the start
        //- and increment points
        inline void resetTime();

        //- Reset to use the current clock value for the increment point
        inline void resetTimeIncrement() const;

        //- The time [seconds] since the start point
        inline double elapsedTime() const;

        //- The time [seconds] since the last call to elapsedTime(),
        //- timeIncrement() or resetTime(), resetTimeIncrement()
        inline double timeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "clockTimeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
