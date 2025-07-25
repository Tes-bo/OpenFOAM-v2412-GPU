/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "Time.H"
#include "argList.H"
#include "Pstream.H"
#include "simpleObjectRegistry.H"
#include "dimensionedConstants.H"
#include "profiling.H"
#include "IOdictionary.H"
#include "fileOperation.H"
#include "fstreamPointer.H"

#include <iomanip>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Output seconds as day-hh:mm:ss
static std::ostream& printTimeHMS(std::ostream& os, double seconds)
{
    const unsigned long ss = seconds;

    // days
    const auto dd = (ss / 86400);

    if (dd) os << dd << '-';

    // hours
    const int hh = ((ss / 3600) % 24);

    if (dd || hh)
    {
        os  << std::setw(2) << std::setfill('0')
            << hh << ':';
    }

    // minutes
    os  << std::setw(2) << std::setfill('0')
        << ((ss / 60) % 60) << ':';

    // seconds
    os  << std::setw(2) << std::setfill('0')
        << (ss % 60);


    // 1/100th seconds. As none or 2 decimal places
    const int hundredths = int(100 * (seconds - ss)) % 100;

    if (hundredths)
    {
        os  << '.' << std::setw(2) << std::setfill('0') << hundredths;
    }

    return os;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Time::readDict()
{
    word application;
    if (controlDict_.readIfPresent("application", application))
    {
        // Do not override if already set so external application can override
        setEnv("FOAM_APPLICATION", application, false);
    }

    // Check for local switches and settings

    const dictionary* localDict = nullptr;

    // DebugSwitches
    if
    (
        (localDict = controlDict_.findDict("DebugSwitches")) != nullptr
     && localDict->size()
    )
    {
        DetailInfo
            << "Overriding DebugSwitches according to "
            << controlDict_.name() << nl;

        debug::debugObjects().setValues(*localDict, true);
    }


    // InfoSwitches
    if
    (
        (localDict = controlDict_.findDict("InfoSwitches")) != nullptr
     && localDict->size()
    )
    {
        DetailInfo
            << "Overriding InfoSwitches according to "
            << controlDict_.name() << nl;

        debug::infoObjects().setValues(*localDict, true);
    }

    // OptimisationSwitches
    if
    (
        (localDict = controlDict_.findDict("OptimisationSwitches")) != nullptr
     && localDict->size()
    )
    {
        DetailInfo
            << "Overriding OptimisationSwitches according to "
            << controlDict_.name() << nl;

        debug::optimisationObjects().setValues(*localDict, true);
    }


    // Handle fileHandler explicitly since it affects local dictionary
    // monitoring.
    word fileHandlerName;
    if
    (
        localDict
     && localDict->readIfPresent("fileHandler", fileHandlerName)
     && fileHandler().type() != fileHandlerName
    )
    {
        DetailInfo << "Overriding fileHandler to " << fileHandlerName << nl;

        // Remove old watches since destroying the file
        fileNameList oldWatched(controlDict_.watchIndices().size());
        forAllReverse(controlDict_.watchIndices(), i)
        {
            const label watchi = controlDict_.watchIndices()[i];
            oldWatched[i] = fileHandler().getFile(watchi);
            fileHandler().removeWatch(watchi);
        }
        controlDict_.watchIndices().clear();

        // Reporting verbosity corresponding to detail level
        const bool verbose = (::Foam::infoDetailLevel > 0);

        // The new handler
        refPtr<fileOperation> newHandler
        (
            fileOperation::New(fileHandlerName, verbose)
        );

        // Install the new handler
        (void) fileOperation::fileHandler(newHandler);

        if (TimePaths::distributed())
        {
            newHandler->distributed(true);
        }

        // Reinstall old watches
        fileHandler().addWatches(controlDict_, oldWatched);
    }


    // DimensionedConstants.
    // - special case since it may change both the 'unitSet' and the
    //   individual values
    if
    (

        (localDict = controlDict_.findDict("DimensionedConstants")) != nullptr
     && localDict->size()
    )
    {
        DetailInfo
            << "Overriding DimensionedConstants according to "
            << controlDict_.name() << nl;

        simpleObjectRegistry& objs = debug::dimensionedConstantObjects();

        // Change in-memory
        dimensionedConstants().merge(*localDict);

        ISpanStream dummyIs;

        // Reporting verbosity corresponding to detail level
        const bool verbose = (::Foam::infoDetailLevel > 0);

        forAllConstIters(objs, iter)
        {
            const List<simpleRegIOobject*>& objects = *iter;

            for (simpleRegIOobject* obj : objects)
            {
                obj->readData(dummyIs);

                if (verbose)
                {
                    Info<< "    ";
                    obj->writeData(Info);
                    Info<< nl;
                }
            }
        }
    }


    // DimensionSets
    if
    (
        (localDict = controlDict_.findDict("DimensionSets")) != nullptr
        && localDict->size()
    )
    {
        DetailInfo
            << "Overriding DimensionSets according to "
            << controlDict_.name() << nl;

        simpleObjectRegistry& objs = debug::dimensionSetObjects();

        dictionary dict(Foam::dimensionSystems());
        dict.merge(*localDict);

        simpleObjectRegistryEntry* objPtr = objs.find("DimensionSets");

        if (objPtr)
        {
            DetailInfo << *localDict << nl;

            const List<simpleRegIOobject*>& objects = *objPtr;

            OCharStream os;

            for (simpleRegIOobject* obj : objects)
            {
                os.rewind();
                os  << dict;

                ISpanStream is(os.view());
                obj->readData(is);
            }
        }
    }


    if (!deltaTchanged_)
    {
        controlDict_.readEntry("deltaT", deltaT_);
    }

    writeControlNames.readIfPresent
    (
        "writeControl",
        controlDict_,
        writeControl_
    );

    scalar oldWriteInterval = writeInterval_;

    if (controlDict_.readIfPresent("writeInterval", writeInterval_))
    {
        if (writeControl_ == wcTimeStep && label(writeInterval_) < 1)
        {
            FatalIOErrorInFunction(controlDict_)
                << "writeInterval < 1 for writeControl timeStep"
                << exit(FatalIOError);
        }
    }
    else
    {
        controlDict_.readEntry("writeFrequency", writeInterval_);
    }


    if (oldWriteInterval != writeInterval_)
    {
        switch (writeControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate writeTimeIndex_ to be in units of current
                // writeInterval.
                writeTimeIndex_ = label
                (
                    writeTimeIndex_
                  * oldWriteInterval
                  / writeInterval_
                );
            break;

            default:
            break;
        }
    }

    if (controlDict_.readIfPresent("purgeWrite", purgeWrite_))
    {
        if (purgeWrite_ < 0)
        {
            WarningInFunction
                << "invalid value for purgeWrite " << purgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            purgeWrite_ = 0;
        }
    }

    format_ =
        IOstreamOption::floatFormatEnum("timeFormat", controlDict_, format_);

    controlDict_.readIfPresent("timePrecision", precision_);

    // stopAt at 'endTime' or a specified value
    // if nothing is specified, the endTime is zero
    if (stopAtControlNames.readIfPresent("stopAt", controlDict_, stopAt_))
    {
        if (stopAt_ == saEndTime)
        {
            controlDict_.readEntry("endTime", endTime_);
        }
        else
        {
            endTime_ = GREAT;
        }
    }
    else if (!controlDict_.readIfPresent("endTime", endTime_))
    {
        endTime_ = 0;
    }

    // Adjust the TimeState name
    dimensionedScalar::name() = timeName(value());

    // Specifying the write version doesn't make sense (2023-10)
    if (controlDict_.found("writeVersion"))
    {
        writeStreamOption_.version(controlDict_.get<token>("writeVersion"));
    }

    // FUTURE? optional control to support command-line option to
    // set the write format and ignore this dictionary entry.

    if (controlDict_.found("writeFormat"))
    {
        writeStreamOption_.format(controlDict_.get<word>("writeFormat"));
    }

    if (controlDict_.found("writePrecision"))
    {
        IOstream::defaultPrecision
        (
            controlDict_.get<unsigned int>("writePrecision")
        );

        Sout.precision(IOstream::defaultPrecision());
        Serr.precision(IOstream::defaultPrecision());

        Pout.precision(IOstream::defaultPrecision());
        Perr.precision(IOstream::defaultPrecision());

        FatalError.stream().precision(IOstream::defaultPrecision());
        FatalIOError.stream().precision(IOstream::defaultPrecision());
    }

    if (controlDict_.found("writeCompression"))
    {
        writeStreamOption_.compression
        (
            controlDict_.get<word>("writeCompression")
        );

        if (writeStreamOption_.compression() == IOstreamOption::COMPRESSED)
        {
            if (writeStreamOption_.format() != IOstreamOption::ASCII)
            {
                IOWarningInFunction(controlDict_)
                    << "Disabled output compression for non-ascii format"
                    << " (inefficient/ineffective)"
                    << endl;

                writeStreamOption_.compression(IOstreamOption::UNCOMPRESSED);
            }
            else if (!ofstreamPointer::supports_gz())
            {
                IOWarningInFunction(controlDict_)
                    << "Disabled output compression"
                    << " (missing libz support)"
                    << endl;

                writeStreamOption_.compression(IOstreamOption::UNCOMPRESSED);
            }
        }
    }

    controlDict_.readIfPresent("graphFormat", graphFormat_);
    controlDict_.readIfPresent("runTimeModifiable", runTimeModifiable_);


    if (!runTimeModifiable_ && controlDict_.watchIndices().size())
    {
        forAllReverse(controlDict_.watchIndices(), i)
        {
            fileHandler().removeWatch(controlDict_.watchIndices()[i]);
        }
        controlDict_.watchIndices().clear();
    }
}


bool Foam::Time::read()
{
    if (controlDict_.regIOobject::read())
    {
        // Read contents
        readDict();
        functionObjects_.read();

        if (runTimeModifiable_)
        {
            // For IOdictionary the call to regIOobject::read() would have
            // already updated all the watchIndices via the addWatch but
            // controlDict_ is an unwatchedIOdictionary so will only have
            // stored the dependencies as files.
            fileHandler().addWatches(controlDict_, controlDict_.files());
        }
        controlDict_.files().clear();

        return true;
    }

    return false;
}


void Foam::Time::readModifiedObjects()
{
    if (runTimeModifiable_)
    {
        // Get state of all monitored objects (=registered objects with a
        // valid filePath).
        // Note: requires same ordering in objectRegistries on different
        // processors!
        fileHandler().updateStates
        (
            (
                IOobject::fileModificationChecking == IOobject::inotifyMaster
             || IOobject::fileModificationChecking == IOobject::timeStampMaster
            ),
            Pstream::parRun()
        );
        // Time handling is special since controlDict_ is the one dictionary
        // that is not registered to any database.

        if (controlDict_.readIfModified())
        {
            readDict();
            functionObjects_.read();

            if (runTimeModifiable_)
            {
                // For IOdictionary the call to regIOobject::read() would have
                // already updated all the watchIndices via the addWatch but
                // controlDict_ is an unwatchedIOdictionary so will only have
                // stored the dependencies as files.

                fileHandler().addWatches(controlDict_, controlDict_.files());
            }
            controlDict_.files().clear();
        }

        bool registryModified = objectRegistry::modified();

        if (registryModified)
        {
            objectRegistry::readModifiedObjects();
        }
    }
}


bool Foam::Time::writeTimeDict() const
{
    addProfiling(writing, "objectRegistry::writeObject");

    const word tmName(timeName());

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            tmName,
            "uniform",
            *this,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );

    timeDict.add("value", timeName(timeToUserTime(value()), maxPrecision_));
    timeDict.add("name", string(tmName));
    timeDict.add("index", timeIndex_);
    timeDict.add("deltaT", timeToUserTime(deltaT_));
    timeDict.add("deltaT0", timeToUserTime(deltaT0_));

    return timeDict.regIOobject::writeObject
    (
        IOstreamOption(IOstreamOption::ASCII),
        true
    );
}


bool Foam::Time::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (writeTime())
    {
        bool writeOK = writeTimeDict();

        if (writeOK)
        {
            writeOK = objectRegistry::writeObject(streamOpt, writeOnProc);
        }

        if (writeOK)
        {
            // Does the writeTime trigger purging?
            if (writeTime_ && purgeWrite_)
            {
                if
                (
                    previousWriteTimes_.empty()
                 || previousWriteTimes_.top() != timeName()
                )
                {
                    previousWriteTimes_.push(timeName());
                }

                while (previousWriteTimes_.size() > purgeWrite_)
                {
                    fileHandler().rmDir
                    (
                        fileHandler().filePath
                        (
                            objectRegistry::path(previousWriteTimes_.pop()),
                            false  // No .gz check (is directory)
                        )
                    );
                }
            }
        }

        return writeOK;
    }

    return false;
}


bool Foam::Time::writeNow()
{
    writeTime_ = true;
    return write();
}


bool Foam::Time::writeAndEnd()
{
    stopAt_  = saWriteNow;
    endTime_ = value();

    return writeNow();
}


void Foam::Time::writeOnce()
{
    writeOnce_ = true;
}


Foam::Ostream& Foam::Time::printExecutionTime(OSstream& os) const
{
    switch (printExecutionFormat_)
    {
        case 1:
        {
            os  << "ExecutionTime = ";
            printTimeHMS(os.stdStream(), elapsedCpuTime());

            os  << "  ClockTime = ";
            printTimeHMS(os.stdStream(), elapsedClockTime());
        }
        break;

        default:
        {
            os  << "ExecutionTime = " << elapsedCpuTime() << " s"
                << "  ClockTime = " << elapsedClockTime() << " s";
        }
        break;
    }

    os  << nl << endl;

    return os;
}


// ************************************************************************* //
