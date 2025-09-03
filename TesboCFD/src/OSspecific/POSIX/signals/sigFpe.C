/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "sigFpe.H"
#include "error.H"
#include "JobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"
#include "UList.H"
#include "Switch.H"
#include <algorithm>
#include <limits>

// File-local functions
#include "signalMacros.C"

#if defined(__linux__) && defined(__GNUC__)
    #ifndef __USE_GNU
        #define __USE_GNU      // To use feenableexcept()
    #endif
    #include <fenv.h>
    #include <malloc.h>
#endif

// Special handling for APPLE
#ifdef __APPLE__
    #include "feexceptErsatz.H"
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::sigFpe::switchFpe_(Foam::debug::optimisationSwitch("trapFpe", 0));
bool Foam::sigFpe::switchNan_(Foam::debug::optimisationSwitch("setNaN", 0));

bool Foam::sigFpe::sigActive_ = false;
bool Foam::sigFpe::nanActive_ = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Can turn on/off via env variable containing a bool (true|false|on|off ...)
// or by the specified flag
static bool isTrue(const char* envName, bool deflt)
{
    Foam::Switch sw(Foam::Switch::find(Foam::getEnv(envName)));

    if (sw.good())
    {
        return static_cast<bool>(sw);
    }

    // Env was not set or did not contain a valid bool value
    return deflt;
}


#ifdef __linux__
extern "C"
{
    extern void* __libc_malloc(size_t size);

    // Override the GLIBC malloc to support filling with NaN
    void* malloc(size_t size)
    {
        // Call the low-level GLIBC malloc function
        void* ptr = __libc_malloc(size);

        if (Foam::sigFpe::nanActive())
        {
            // Fill with signaling_NaN
            const auto val = std::numeric_limits<Foam::scalar>::signaling_NaN();

            // Can dispatch with
            // - std::execution::parallel_unsequenced_policy
            // - std::execution::unsequenced_policy
            std::fill_n
            (
                reinterpret_cast<Foam::scalar*>(ptr),
                (size/sizeof(Foam::scalar)),
                val
            );
        }

        return ptr;
    }
} // End extern C

#endif  // __linux__


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigFpe::sigHandler(int)
{
    #if (defined(__linux__) && defined(__GNUC__)) || defined(__APPLE__)

    resetHandler("SIGFPE", SIGFPE);

    JobInfo::shutdown();        // From running -> finished
    error::printStack(Perr);
    ::raise(SIGFPE);            // Throw signal (to old handler)

    #endif  // (__linux__ && __GNUC__) || __APPLE__
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    set(false);  // false = non-verbose
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    unset(false);  // false = non-verbose
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sigFpe::requested()
{
    return isTrue("FOAM_SIGFPE", switchFpe_);
}


void Foam::sigFpe::set(bool verbose)
{
    if (!sigActive_ && requested())
    {
        #if (defined(__linux__) && defined(__GNUC__)) || defined(__APPLE__)

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        setHandler("SIGFPE", SIGFPE, sigHandler);

        sigActive_ = true;
        #endif

        if (verbose)
        {
            Info<< "trapFpe: Floating point exception trapping ";

            if (sigActive_)
            {
                Info<< "enabled (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "- not supported on this platform" << endl;
            }
        }
    }


    nanActive_ = false;
    if (isTrue("FOAM_SETNAN", switchNan_))
    {
        #ifdef __linux__
        nanActive_ = true;
        #endif

        if (verbose)
        {
            Info<< "setNaN : Fill allocated memory with NaN ";

            if (nanActive_)
            {
                Info<< "enabled (FOAM_SETNAN)." << endl;
            }
            else
            {
                Info<< " - not supported on this platform" << endl;
            }
        }
    }
}


void Foam::sigFpe::unset(bool verbose)
{
    #if (defined(__linux__) && defined(__GNUC__)) || defined(__APPLE__)
    if (sigActive_)
    {
        if (verbose)
        {
            Info<< "sigFpe : Disabling floating point exception trapping"
                << endl;
        }

        resetHandler("SIGFPE", SIGFPE);

        // Reset exception raising
        const int oldExcept = fedisableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        if (oldExcept == -1)
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }

        sigActive_ = false;
    }
    #endif

    nanActive_ = false;
}


void Foam::sigFpe::fillNan(char* buf, size_t count)
{
    if (!buf || !count) return;

    // Fill with signaling_NaN
    const auto val = std::numeric_limits<scalar>::signaling_NaN();

    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    std::fill_n
    (
        reinterpret_cast<scalar*>(buf), (count/sizeof(scalar)), val
    );
}


void Foam::sigFpe::fillNan(UList<scalar>& list)
{
    if (list.empty()) return;

    // Fill with signaling_NaN
    const auto val = std::numeric_limits<scalar>::signaling_NaN();

    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    std::fill_n
    (
        list.data(), list.size(), val
    );
}


// ************************************************************************* //
