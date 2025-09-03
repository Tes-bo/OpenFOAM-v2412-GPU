/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"
#include "Time.H"

// File-local functions
#include "signalMacros.C"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Signal number to catch
int Foam::sigWriteNow::signal_
(
    Foam::debug::optimisationSwitch("writeNowSignal", -1)
);

// Pointer to Time (file-local variable)
static Foam::Time* runTimePtr_ = nullptr;


// * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * * //

namespace Foam
{

// Register re-reader
struct addwriteNowSignalToOpt
:
    public ::Foam::simpleRegIOobject
{
    addwriteNowSignalToOpt(const addwriteNowSignalToOpt&) = delete;

    void operator=(const addwriteNowSignalToOpt&) = delete;

    explicit addwriteNowSignalToOpt(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
    {}

    virtual ~addwriteNowSignalToOpt() = default;

    virtual void readData(Foam::Istream& is)
    {
        is >> sigWriteNow::signal_;
        sigWriteNow::set(true);
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << sigWriteNow::signal_;
    }
};

addwriteNowSignalToOpt addwriteNowSignalToOpt_("writeNowSignal");

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigWriteNow::sigHandler(int)
{
    if (runTimePtr_)
    {
        Info<< "sigWriteNow :"
            << " setting up write at end of the next iteration" << nl << endl;
        runTimePtr_->writeOnce();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigWriteNow::sigWriteNow(Time& runTime, bool verbose)
{
    runTimePtr_ = &runTime;  // Store Time reference
    set(verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigWriteNow::~sigWriteNow()
{
    if (!active())
    {
        return;
    }

    resetHandler("writeNow", signal_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigWriteNow::set(bool verbose)
{
    if (!active())
    {
        return;
    }

    if (verbose)
    {
        Info<< "sigWriteNow :"
            << " Enabling writing upon signal " << signal_ << nl;
    }

    setHandler("writeNow", signal_, sigHandler);
}


// ************************************************************************* //
