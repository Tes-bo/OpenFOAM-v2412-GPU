/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "timer.H"
#include "error.H"
#include "defineDebugSwitch.H"

#include <unistd.h>

// File-local functions
#include "signalMacros.C"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(timer, "timer", 0);
    registerDebugSwitchWithName(timer, timer, "timer");
}

jmp_buf Foam::timer::envAlarm;

unsigned int Foam::timer::oldTimeOut_ = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::timer::sigHandler(int)
{
    DebugInFunction<< "Timed out. Jumping." << endl;

    longjmp(envAlarm, 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timer::timer(unsigned int seconds)
:
    timeOut_(seconds)
{
    if (!timeOut_)
    {
        return;
    }

    // Singleton since handler is static function
    if (oldTimeOut_)
    {
        FatalErrorInFunction
            << "timer already used."
            << abort(FatalError);
    }

    // Set alarm signal handler
    // - do not block any signals while in it
    // - clear list of signals to mask

    setHandler("SIGALRM", SIGALRM, sigHandler);

    // Set alarm timer
    oldTimeOut_ = ::alarm(timeOut_);

    DebugInFunction
        << "Installing timeout " << int(timeOut_) << " seconds"
        << " (overriding old timeout " << int(oldTimeOut_) << ")." << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timer::~timer()
{
    if (!timeOut_)
    {
        return;
    }

    DebugInFunction
        << "timeOut=" << int(timeOut_)
        << " : resetting timeOut to " << int(oldTimeOut_) << endl;

    // Reset alarm timer
    ::alarm(oldTimeOut_);
    oldTimeOut_ = 0;

    resetHandler("SIGALRM", SIGALRM);
}


// ************************************************************************* //
