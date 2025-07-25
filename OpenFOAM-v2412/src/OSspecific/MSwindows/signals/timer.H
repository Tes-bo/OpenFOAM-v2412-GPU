/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::timer

Description
    Implements a timeout mechanism via sigalarm.

    Example usage:
    \code
        timer myTimer(5);     // 5 sec
        ..
        if (timedOut(myTimer))
        {
            // timed out
        }
        else
        {
            // do something possible blocking
        }
    \endcode

    Constructor set signal handler on sigalarm and alarm(). Destructor
    clears these.

Warning
    The setjmp restores complete register state so including local vars
    held in regs. So if in blocking part something gets calced in a stack
    based variable make sure it is declared 'volatile'.

Note
    timedOut is macro because setjmp can't be in member function of timer.
    ?something to do with stack frames.

SourceFiles
    timer.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_timer_H
#define Foam_timer_H

#include <csetjmp>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Check if timeout has occurred
//  keep setjmp in same stack frame so no function calls
#undef  timedOut
#define timedOut(x) \
    ((x).timeOut_ ? setjmp(Foam::timer::envAlarm) : false)

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timer Declaration
\*---------------------------------------------------------------------------*/

class timer
{
    // Private Data

        //- Old alarm() value
        static unsigned int oldTimeOut_;


    // Private Member Functions

        //- Alarm handler
        static void sigHandler(int);


public:

    // Public Data

        //- Named/registered debug switch: 'timer'
        static int debug;

        //- The time-out value (seconds). Needed by macro timedOut
        unsigned int timeOut_;

        //- State for setjmp. Needed by macro timedOut
        static jmp_buf envAlarm;


    // Constructors

        //- Construct with specified time-out, a value of 0 makes it a no-op.
        explicit timer(unsigned int seconds);


    //- Destructor. Restores the alarm and signal handler as required.
    ~timer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
