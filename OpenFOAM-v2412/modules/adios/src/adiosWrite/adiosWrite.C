/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "adiosWrite.H"
#include "dictionary.H"
#include "FlatOutput.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(adiosWrite, 0);
    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        adiosWrite,
        dictionary
    );
} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::adiosWrite::restart()
{
    if (restartType_ == restartDisabled)
    {
        return true;
    }

    if (restartType_ == restartNone)
    {
        restartType_ = restartDisabled;

        // Not a restart,
        // so write initial conditions (including mesh) on first time
        Info<< type()
            << " : Restart none, write initial state for future restart"
            << endl;

        write();
        return true;
    }


    // Restart from a specified time, or from latest

    instantList adiosTimes;

    if (UPstream::master())
    {
        adiosTimes = adiosFoam::findTimes(dataDir_);
    }
    #if (OPENFOAM < 2206)
    Pstream::scatter(adiosTimes);  // Older version for broadcast
    #else
    Pstream::broadcast(adiosTimes);
    #endif


    restartIndex_ = -1;
    if (restartType_ == restartTime)
    {
        restartType_  = restartDisabled;
        restartIndex_ = Time::findClosestTimeIndex(adiosTimes, restartTime_);

        Info<< type()
            << " :  Restart from " << restartTime_ << endl;

        if (restartIndex_ < 0)
        {
            FatalErrorInFunction
                << "No appropriate ADIOS restart found for time "
                << restartTime_
                << exit(FatalIOError);
        }
    }
    else if (restartType_ == restartLatest)
    {
        restartType_  = restartDisabled;
        restartIndex_ = adiosTimes.size() - 1;

        Info<< type()
            << " : Restart from latestTime" << nl;

        if (restartIndex_  < 0)
        {
            Info<< "    No adios files found" << endl;
            return true;
        }
    }
    else
    {
        // Can only hit this as programmer error
        FatalErrorInFunction
            << "bad value for restart status"
            << exit(FatalIOError);
    }


    Info<< type()
        << " : Restart from adios " << adiosTimes[restartIndex_]
        << endl;


    adiosFoam::adiosTime bpTime = readData(adiosTimes[restartIndex_]);

    if (bpTime.good())
    {
        Time& t = const_cast<Time&>(time());
        t.setTime(bpTime.timeValue(), bpTime.timeIndex());

        // TODO: handle deltaT0 too
        t.setDeltaT(bpTime.deltaT(), false);

        restartIndex_ = bpTime.timeIndex();

        Info<< "# Time = " << t.timeOutputValue()
            << " index = " << t.timeIndex()
            << " deltaT = " << bpTime.deltaT() << nl;

        // pointTime_ = xxx, faceTime_ = xxx
    }
    else
    {
        FatalErrorInFunction
            << "Restart reading failed for time " << restartTime_
            << exit(FatalIOError);
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::adiosWrite::adiosWrite
(
    const word& groupName,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::timeFunctionObject("adiosWrite", runTime),
    adiosCoreWrite(dict),
    dataDir_(),
    restartType_(restartNone),
    stopAt_(stopNone),
    restartTime_(VGREAT),
    stopTime_(VGREAT),
    restartIndex_(-1)
{
    // Make write controls an absolute requirement.
    //
    // The presence of 'writeControl' triggers loading via
    // functionObjects::timeControl

    // Exit with FatalIOError, to avoid it being trapped/ignored by the
    // functionObjectList.

    if (!dict.found("writeControl"))
    {
        FatalIOErrorInFunction(dict)
            << "Invoked " << type()
            << " without specifying a 'writeControl'"
            << exit(FatalIOError);
    }
    if (!dict.found("writeInterval"))
    {
        FatalIOErrorInFunction(dict)
            << "Invoked " << type()
            << " without specifying a 'writeInterval'"
            << exit(FatalIOError);
    }

    // Read dictionary
    read(dict);

    if (postProcess)
    {
        restartType_ = restartDisabled;
        stopAt_      = stopDisabled;
    }
    else if (restartType_ == restartNone)
    {
        restart();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::adiosWrite::~adiosWrite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::adiosWrite::read(const dictionary& dict)
{
    dataDir_.clear();
    if (dict.readIfPresent("dataDirectory", dataDir_))
    {
        dataDir_.expand();
    }
    if (dataDir_.empty())
    {
        dataDir_ = adiosFoam::defaultDirectory;
    }

    regions_.clear();

    // All known regions
    // - sorted order for consistency
    wordList regionNames(time().sortedNames<fvMesh>());

    // Info<< "Known regions: " << regionNames << endl;

    // Top-level dictionaries
    dictionary regionsDict(dict.subOrEmptyDict("regions", keyType::LITERAL));


    // Verify that "regions" only contains sub-dictionaries
    // and detect any unknown region names
    {
        wordList bad(regionsDict.size());
        wordList warn(regionsDict.size());

        label nbad = 0;
        label nwarn = 0;

        for (const entry& dEntry : regionsDict)
        {
            const word& regName = dEntry.keyword();

            if (!regionNames.found(regName))
            {
                warn[nwarn++] = regName;
            }
            if (!dEntry.isDict())
            {
                bad[nbad++] = regName;
            }
        }

        if (nwarn)
        {
            warn.resize(nwarn);

            WarningInFunction
                << warn.size() << " unknown regions specified in adiosDict:"
                << nl
                << "   " << flatOutput(warn) << endl;
        }

        if (nbad)
        {
            bad.resize(nbad);

            FatalErrorInFunction
                << bad.size() << " bad region specifications in adiosDict:"
                << nl
                << "   " << flatOutput(bad) << nl
                << "must be dictionary format with region-name as keyword" << nl
                << exit(FatalIOError);
        }
    }


    // For each region
    // - uses "explicit" entry and "write" and "ignore" sub-dictionaries
    for (const word& regName : regionNames)
    {
        // Cannot fail since the region name came from previous lookup
        const fvMesh& mesh = time().lookupObject<fvMesh>(regName);

        const entry* dictPtr =
            regionsDict.findEntry(regName, keyType::LITERAL);

        if (dictPtr)
        {
            // Make copy and merge in region-specific content
            // (overwrite existing entries)

            dictionary mergedDict(regionsDict);
            mergedDict <<= dictPtr->dict();

            regions_.append(adiosFoam::regionControl(mesh, mergedDict));
        }
        else
        {
            regions_.append(adiosFoam::regionControl(mesh, dict));
        }
    }


    // Print info to terminal

    // Only allowed to change restart status if not already triggered
    restartTime_ = VGREAT;
    if (restartType_ != restartDisabled)
    {
        // Default is none (no restart)
        const word when = dict.getOrDefault<word>("restartFrom", "none");

        if (when == "none")
        {
            restartType_ = restartNone;

            Info<< type()
                << " : No restart requested" << endl;
        }
        else if (when == "restartTime")
        {
            restartType_ = restartTime;
            dict.readEntry("restartTime", restartTime_);

            Info<< type()
                << " : Restart requested, Time = " << restartTime_ << endl;
        }
        else if (when == "latestTime")
        {
            restartType_ = restartLatest;

            Info<< type()
                << " restart from latestTime" << nl;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Expected 'restartFrom' as none, restartTime, or latestTime"
                << " found '" << when << "'"
                << exit(FatalIOError);
        }
    }

    stopTime_ = VGREAT;
    if (stopAt_ != stopDisabled)
    {
        // Default is none (no stopping control)
        const word when = dict.getOrDefault<word>("stopAt", "none");

        if (when == "none")
        {
            stopAt_ = stopNone;

            Info<< type()
                << " : No stopAt requested" << endl;
        }
        else if (when == "now")
        {
            stopAt_ = stopNow;

            Info<< type()
                << ": Stop as soon as possible" << nl;
        }
        else if (when == "stopTime")
        {
            stopAt_ = stopTime;
            dict.readEntry("stopTime", stopTime_);

            Info<< type()
                << " : Stop requested, Time = " << stopTime_ << endl;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Expected 'stopAt' as none, now, or stopTime"
                << " found '" << when << "'"
                << exit(FatalIOError);
        }
    }

    return true;
}


bool Foam::functionObjects::adiosWrite::open(const Time& t)
{
    // Create output directory if initially non-existent
    static bool checkdir = true;

    if (checkdir && UPstream::master() && !isDir(dataDir_))
    {
        mkDir(dataDir_);
    }
    checkdir = false;

    if (adiosCoreWrite::open(dataDir_/t.timeName()))
    {
        return true;
    }

    return false;
}


bool Foam::functionObjects::adiosWrite::execute()
{
    return restart();
}


bool Foam::functionObjects::adiosWrite::write()
{
    if (restartIndex_ == time().timeIndex())
    {
        return true;
    }

    Info<< type()
        << " : Writing ADIOS data for time "
        << time().timeName() << " index "
        << time().timeIndex() << endl;

    writeData();

    switch (stopAt_)
    {
        case stopNow:
        {
            if (time().stopAt(Time::saNoWriteNow))
            {
                Info<< "USER REQUESTED ABORT (timeIndex="
                    << time_.timeIndex()
                    << "): stop without writing data"
                    << endl;
            }
            break;
        }

        case stopTime:
        {
            if (time().timeOutputValue() >= stopTime_)
            {
                if (time_.stopAt(Time::saNoWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop without writing data"
                        << endl;
                }
            }
            break;
        }

        default:
            break;

    }

    return true;
}


// ************************************************************************* //
