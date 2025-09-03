/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2011 Symscape
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
#include "MSwindows.H"
#undef DebugInfo        // Windows name clash with OpenFOAM messageStream

#define WIN32_LEAN_AND_MEAN
#undef  WINVER
#define WINVER 0x0500   // To access CreateTimerQueueTimer
#include <windows.h>

// File-local functions
#include "signalMacros.C"

#undef  SIGALRM
#define SIGALRM 14


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(timer, "timer", 0);
    registerDebugSwitchWithName(timer, timer, "timer");
}

jmp_buf Foam::timer::envAlarm;

unsigned int Foam::timer::oldTimeOut_ = 0;

static HANDLE hTimer_ = nullptr;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static VOID CALLBACK timerExpired(PVOID lpParam, BOOLEAN TimerOrWaitFired)
{
    ::raise(SIGALRM);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::timer::sigHandler(int)
{
    DebugInFunction << "Timed out. Jumping." << endl;

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
    if (hTimer_)
    {
        FatalErrorInFunction
            << "timer already used."
            << abort(FatalError);
    }

    // Set alarm signal handler
    setHandler("SIGALRM", SIGALRM, sigHandler);

    // Set alarm timer
    const bool ok = ::CreateTimerQueueTimer
    (
        &hTimer_,
        nullptr,
        static_cast<WAITORTIMERCALLBACK>(timerExpired),
        nullptr,
        timeOut_ * 1000,
        0,
        0
    );

    if (!ok)
    {
        hTimer_ = nullptr;
        FatalErrorInFunction
            << "CreateTimerQueueTimer, "
            << MSwindows::lastError() << nl
            << abort(FatalError);
    }

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
    const bool ok = ::DeleteTimerQueueTimer(nullptr, hTimer_, nullptr);

    hTimer_ = nullptr;

    if (!ok)
    {
        FatalErrorInFunction
            << "DeleteTimerQueueTimer, "
            << MSwindows::lastError() << nl
            << abort(FatalError);
    }

    resetHandler("SIGALRM", SIGALRM);
}


// ************************************************************************* //
