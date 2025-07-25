/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "sigWriteNow.H"
#include "sigStopAtWriteNow.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"
#include "Time.H"

// File-local functions
#include "signalMacros.C"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Signal number to catch
int Foam::sigStopAtWriteNow::signal_
(
    Foam::debug::optimisationSwitch("stopAtWriteNowSignal", -1)
);

// Pointer to Time (file-local variable)
static Foam::Time const* runTimePtr_ = nullptr;


// * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * * //

namespace Foam
{
// Register re-reader
struct addstopAtWriteNowSignalToOpt
:
    public ::Foam::simpleRegIOobject
{
    addstopAtWriteNowSignalToOpt(const addstopAtWriteNowSignalToOpt&) = delete;

    void operator=(const addstopAtWriteNowSignalToOpt&) = delete;

    explicit addstopAtWriteNowSignalToOpt(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
    {}

    virtual ~addstopAtWriteNowSignalToOpt() = default;

    virtual void readData(Foam::Istream& is)
    {
        is >> sigStopAtWriteNow::signal_;
        sigStopAtWriteNow::set(true);
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << sigStopAtWriteNow::signal_;
    }
};

addstopAtWriteNowSignalToOpt addstopAtWriteNowSignalToOpt_
(
    "stopAtWriteNowSignal"
);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigStopAtWriteNow::sigHandler(int)
{
    resetHandler("stopAtWriteNow", signal_);

    JobInfo::shutdown();        // From running -> finished

    if (runTimePtr_)
    {
        Info<< "sigStopAtWriteNow :"
            << " setting up write and stop at end of the next iteration"
            << nl << endl;
        runTimePtr_->stopAt(Time::saWriteNow);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigStopAtWriteNow::sigStopAtWriteNow(const Time& runTime, bool verbose)
{
    runTimePtr_ = &runTime; // Store runTime
    set(verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigStopAtWriteNow::~sigStopAtWriteNow()
{
    if (!active())
    {
        return;
    }

    resetHandler("stopAtWriteNow", signal_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigStopAtWriteNow::set(bool verbose)
{
    if (!active())
    {
        return;
    }


    // Check that the signal is different from the writeNowSignal
    if (sigWriteNow::signalNumber() == signal_)
    {
        FatalErrorInFunction
            << "stopAtWriteNowSignal : " << signal_
            << " cannot be the same as the writeNowSignal."
            << " Please change this in the etc/controlDict."
            << exit(FatalError);
    }

    if (verbose)
    {
        Info<< "sigStopAtWriteNow :"
            << " Enabling writing and stopping upon signal " << signal_
            << endl;
    }

    setHandler("stopAtWriteNow", signal_, sigHandler);
}


// ************************************************************************* //
