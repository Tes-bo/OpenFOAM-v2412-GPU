/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "externalFileCoupler.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "Switch.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalFileCoupler, 0);
}

Foam::word Foam::externalFileCoupler::lockName = "OpenFOAM";


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read file contents and return a stop control as follows:
// - contains "done" (should actually be status=done, but we are generous) :
//   The master (OpenFOAM) has signalled that it is done. Report as <endTime>
//
// - action=writeNow, action=nextWrite action=noWriteNow :
//   The slave has signalled that it is done and wants the master to exit with
//   the specified type of action. Report as corresponding <action>.
//
// Anything else (empty file, no action=, etc) is reported as <unknown>.
//
static enum Time::stopAtControls getStopAction(const std::string& filename)
{
    // Slurp entire input file (must exist) as a single string
    string fileContent;

    std::ifstream is(filename);
    std::getline(is, fileContent, '\0');

    if (fileContent.contains("done"))
    {
        return Time::stopAtControls::saEndTime;
    }

    const auto equals = fileContent.find('=');

    if (equals != std::string::npos)
    {
        const word actionName(word::validate(fileContent.substr(equals+1)));

        return
            Time::stopAtControlNames.lookup
            (
                actionName,
                Time::stopAtControls::saUnknown
            );
    }

    return Time::stopAtControls::saUnknown;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalFileCoupler::externalFileCoupler()
:
    runState_(NONE),
    commsDir_("<case>/comms"),
    statusDone_("done"),
    waitInterval_(1u),
    timeOut_(100u),
    slaveFirst_(false),
    log(false)
{
    commsDir_.expand();
    commsDir_.clean();  // Remove unneeded ".."
}


Foam::externalFileCoupler::externalFileCoupler(const fileName& commsDir)
:
    runState_(NONE),
    commsDir_(commsDir),
    statusDone_("done"),
    waitInterval_(1u),
    timeOut_(100u),
    slaveFirst_(false),
    log(false)
{
    commsDir_.expand();
    commsDir_.clean();  // Remove unneeded ".."

    if (Pstream::master())
    {
        mkDir(commsDir_);
    }
}


Foam::externalFileCoupler::externalFileCoupler(const dictionary& dict)
:
    externalFileCoupler()
{
    readDict(dict);

    if (Pstream::master())
    {
        mkDir(commsDir_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalFileCoupler::~externalFileCoupler()
{
    shutdown();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::externalFileCoupler::readDict(const dictionary& dict)
{
    // Normally cannot change directory or initialization
    // if things have already been initialized
    if (!initialized())
    {
        dict.readEntry("commsDir", commsDir_);
        commsDir_.expand();
        commsDir_.clean();  // Remove unneeded ".."
        statusDone_ = dict.getOrDefault<word>("statusDone", "done");
        slaveFirst_ = dict.getOrDefault("initByExternal", false);

        Info<< type() << ": initialize" << nl
            << "    directory: " << commsDir_ << nl
            << "    slave-first: " << Switch(slaveFirst_) << endl;
    }

    waitInterval_ = dict.getOrDefault("waitInterval", 1u);
    if (!waitInterval_)
    {
        // Enforce non-zero sleep
        waitInterval_ = 1u;
    }

    timeOut_ = dict.getOrDefault("timeOut", 100*waitInterval_);

    log = dict.getOrDefault("log", false);

    return true;
}


enum Foam::Time::stopAtControls
Foam::externalFileCoupler::useMaster(const bool wait) const
{
    const bool wasInit = initialized();
    runState_ = MASTER;

    if (Pstream::master())
    {
        if (!wasInit)
        {
            // First time
            mkDir(commsDir_);
        }

        const fileName lck(lockFile());

        // Create lock file - only if it doesn't already exist
        if (!Foam::isFile(lck))
        {
            Log << type()
                << ": creating lock file with status=openfoam" << endl;

            std::ofstream os(lck);
            os << "status=openfoam\n";
        }
    }

    if (wait)
    {
        return waitForMaster();
    }

    return Time::stopAtControls::saUnknown;
}


enum Foam::Time::stopAtControls
Foam::externalFileCoupler::useSlave(const bool wait) const
{
    const bool wasInit = initialized();
    runState_ = SLAVE;

    if (Pstream::master())
    {
        if (!wasInit)
        {
            // First time
            mkDir(commsDir_);
        }

        Log << type() << ": removing lock file" << endl;

        Foam::rm(lockFile());
    }

    if (wait)
    {
        return waitForSlave();
    }

    return Time::stopAtControls::saUnknown;
}


enum Foam::Time::stopAtControls
Foam::externalFileCoupler::waitForMaster() const
{
    if (!initialized())
    {
        useMaster(); // was not initialized
    }

    auto action = Time::stopAtControls::saUnknown;

    if (Pstream::master())
    {
        const fileName lck(lockFile());

        double prevTime = 0;
        double modTime = 0;

        // Wait until file disappears (modTime == 0)
        // But also check for status=done content in the file
        while ((modTime = highResLastModified(lck)) > 0)
        {
            if (prevTime < modTime)
            {
                prevTime = modTime;

                if (Time::stopAtControls::saEndTime == getStopAction(lck))
                {
                    // Found 'done' - slave should not wait for master
                    action = Time::stopAtControls::saEndTime;
                    break;
                }
            }
            sleep(waitInterval_);
        }
    }

    // Send to sub-ranks. Also acts as an MPI barrier
    int intAction(action);
    Pstream::broadcast(intAction);

    return Time::stopAtControls(intAction);
}


enum Foam::Time::stopAtControls
Foam::externalFileCoupler::waitForSlave() const
{
    if (!initialized())
    {
        useSlave(); // was not initialized
    }

    auto action = Time::stopAtControls::saUnknown;

    if (Pstream::master())
    {
        const fileName lck(lockFile());
        unsigned totalTime = 0;

        Log << type() << ": waiting for lock file to appear " << lck << endl;

        while (!Foam::isFile(lck))
        {
            sleep(waitInterval_);

            if (timeOut_ && (totalTime += waitInterval_) > timeOut_)
            {
                FatalErrorInFunction
                    << "Wait time exceeded timeout of " << timeOut_
                    << " s" << abort(FatalError);
            }

            Log << type() << ": wait time = " << totalTime << endl;
        }

        action = getStopAction(lck);

        Log << type() << ": found lock file " << lck << endl;
    }

    // Send to sub-ranks. Also acts as an MPI barrier
    int intAction(action);
    Pstream::broadcast(intAction);

    return Time::stopAtControls(intAction);
}


void Foam::externalFileCoupler::readDataMaster()
{}


void Foam::externalFileCoupler::readDataSlave()
{}


void Foam::externalFileCoupler::writeDataMaster() const
{}


void Foam::externalFileCoupler::writeDataSlave() const
{}


void Foam::externalFileCoupler::removeDataMaster() const
{}


void Foam::externalFileCoupler::removeDataSlave() const
{}


void Foam::externalFileCoupler::shutdown() const
{
    if (Pstream::master() && runState_ == MASTER && Foam::isDir(commsDir_))
    {
        Log << type() << ": lock file status=" << statusDone_ << endl;

        std::ofstream os(lockFile());
        os << "status=" << statusDone_ << nl;
    }

    runState_ = DONE;   // Avoid re-triggering in destructor
}


// ************************************************************************* //
