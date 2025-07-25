/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::cpuTimePosix

Description
    Starts timing CPU usage and return elapsed time from start.

See also
    clockTime

SourceFiles
    cpuTimePosix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cpuTimePosix_H
#define Foam_cpuTimePosix_H

#include <string>
#include <ctime>
#include <sys/times.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cpuTimePosix Declaration
\*---------------------------------------------------------------------------*/

class cpuTimePosix
{
    //- Time structure used, with additional methods
    struct value_type : tms
    {
        //- Construct with the current clock time
        value_type();

        //- Update with the current clock time
        void update();
    };


    // Private Data

        //- Start time, at the time of construction
        value_type start_;

        //- Last time when elapsedCpuTime or cpuTimeIncrement was called.
        //- Also affected by resetCpuTime and resetCpuTimeIncrement.
        mutable value_type last_;


    // Private Member Functions

        //- Difference between two times (a - b)
        inline static double diff(const value_type& a, const value_type& b);


public:

    // Constructors

        //- Construct with the current clock time
        cpuTimePosix();


    // Member Functions

        //- Reset to use the current time for the start time
        void resetCpuTime();

        //- Reset to use the current time for the increment point
        void resetCpuTimeIncrement() const;

        //- Return CPU time [seconds] from the start
        double elapsedCpuTime() const;

        //- Return CPU time [seconds] since last call to cpuTimeIncrement(),
        //- resetCpuTimeIncrement().
        double cpuTimeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
