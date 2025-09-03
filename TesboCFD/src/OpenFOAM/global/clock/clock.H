/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::clock

Description
    Read access to the system clock with formatting.

SourceFiles
    clock.C

\*---------------------------------------------------------------------------*/

#ifndef clock_H
#define clock_H

#include <ctime>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clock Declaration
\*---------------------------------------------------------------------------*/

class clock
{
    // Private Data

        //- Time point at start (in seconds)
        time_t start_;

        //- Last time when elapsedClockTime or clockTimeIncrement was called
        mutable time_t last_;


public:

    // Constructors

        //- Construct with the current clock time for the start point
        clock();


    // Low-level Member Functions

        //- Get the current clock time in seconds
        static time_t getTime();

        //- The current wall-clock date as a raw struct
        //  \deprecated(2020-05) may be removed in future versions
        static const struct tm rawDate();


    // Static Member Functions

        //- The current wall-clock date/time (in local time) as a string
        //- in ISO-8601 format (yyyy-mm-ddThh:mm:ss).
        //  Without time-zone information.
        static std::string dateTime();

        //- The current wall-clock date as a string formatted as
        //- (MON dd yyyy), where MON is Jan, Feb, etc.
        static std::string date();

        //- The current wall-clock (in local time) as a string formatted as
        //- as (hh:mm:ss).
        //  Without time-zone information.
        static std::string clockTime();


    // Member Functions

        //- Returns wall-clock time since clock instantiation
        double elapsedClockTime() const;

        //- Returns wall-clock time since last clockTimeIncrement() call
        double clockTimeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
