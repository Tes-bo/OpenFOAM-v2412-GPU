/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "fieldAverage.H"
#include "volFields.H"
#include "fieldAverageItem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldAverage, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldAverage::initialize()
{
    for (fieldAverageItem& item : faItems_)
    {
        // Note: not clearing data needed for restart
        item.clear(obr(), false);
    }

    Log << type() << " " << name() << ":" << nl;

    // Add mean fields to the field lists
    for (fieldAverageItem& item : faItems_)
    {
        addMeanField<scalar>(item);
        addMeanField<vector>(item);
        addMeanField<sphericalTensor>(item);
        addMeanField<symmTensor>(item);
        addMeanField<tensor>(item);
    }

    // Add prime-squared mean fields to the field lists
    for (fieldAverageItem& item : faItems_)
    {
        addPrime2MeanField<scalar, scalar>(item);
        addPrime2MeanField<vector, symmTensor>(item);
    }

    // Add window fields to the field lists
    for (const fieldAverageItem& item : faItems_)
    {
        restoreWindowFields<scalar>(item);
        restoreWindowFields<vector>(item);
        restoreWindowFields<sphericalTensor>(item);
        restoreWindowFields<symmTensor>(item);
        restoreWindowFields<tensor>(item);
    }


    for (const fieldAverageItem& item : faItems_)
    {
        if (!item.active())
        {
            WarningInFunction
                << "Field " << item.fieldName()
                << " not found in database for averaging";
        }
    }

    // Ensure first averaging works unconditionally
    prevTimeIndex_ = -1;

    Log << endl;
    initialised_ = true;
}


void Foam::functionObjects::fieldAverage::restart()
{
    Log << "    Restarting averaging at time "
        << obr().time().timeOutputValue()
        << nl << endl;

    for (fieldAverageItem& item : faItems_)
    {
        item.clear(obr(), true);
    }

    initialize();
}


void Foam::functionObjects::fieldAverage::calcAverages()
{
    if (!initialised_)
    {
        initialize();
    }

    const label currentTimeIndex = obr().time().timeIndex();
    const scalar currentTime = obr().time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }
    else
    {
        prevTimeIndex_ = currentTimeIndex;
    }

    bool doRestart = false;
    if (periodicRestart_)
    {
        const scalar deltaT = obr().time().deltaTValue();
        const scalar nextTime = restartPeriod_*periodIndex_ + 0.5*deltaT;

        if (currentTime > nextTime)
        {
            doRestart = true;
            ++periodIndex_;
        }
    }

    if (currentTime >= restartTime_)
    {
        doRestart = true;      // Restart is overdue.
        restartTime_ = GREAT;  // Avoid triggering again
    }

    if (doRestart)
    {
        restart();
    }

    Log << type() << " " << name() << " write:" << nl
        << "    Calculating averages" << nl;

    forAll(faItems_, fieldi)
    {
        faItems_[fieldi].evolve(obr());
    }

    storeWindowFields<scalar>();
    storeWindowFields<vector>();
    storeWindowFields<sphericalTensor>();
    storeWindowFields<symmTensor>();
    storeWindowFields<tensor>();

    addMeanSqrToPrime2Mean<scalar, scalar>();
    addMeanSqrToPrime2Mean<vector, symmTensor>();

    calculateMeanFields<scalar>();
    calculateMeanFields<vector>();
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculatePrime2MeanFields<scalar, scalar>();
    calculatePrime2MeanFields<vector, symmTensor>();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAverages() const
{
    Log << "    Writing average fields" << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAveragingProperties()
{
    for (const fieldAverageItem& item : faItems_)
    {
        dictionary propsDict;
        item.writeState(propsDict);
        setProperty(item.fieldName(), propsDict);
    }
}


void Foam::functionObjects::fieldAverage::readAveragingProperties()
{
    if (restartOnRestart_ || restartOnOutput_)
    {
        Info<< "    Starting averaging at time "
            << obr().time().timeOutputValue()
            << nl;
    }
    else
    {
        Info<< "    Restarting averaging for fields:" << nl;


        for (fieldAverageItem& item : faItems_)
        {
            const word& fieldName = item.fieldName();
            if (foundProperty(fieldName))
            {
                dictionary fieldDict;
                getDict(fieldName, fieldDict);
                item.readState(fieldDict);

                if (item.allowRestart())
                {
                    scalar userTotalTime =
                        obr().time().timeToUserTime(item.totalTime());

                    Info<< "        " << fieldName
                        << ": iters = " << item.totalIter()
                        << " time = " << userTotalTime << nl;
                }
                else
                {
                    item.clear(obr(), true);

                    Info<< "        " << fieldName
                        << ": starting averaging at time "
                        << obr().time().timeOutputValue() << endl;
                }
            }
            else
            {
                Info<< "        " << fieldName
                    << ": starting averaging at time "
                    << obr().time().timeOutputValue() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::fieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    initialised_(false),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(GREAT),
    restartTime_(GREAT),
    faItems_(),
    periodIndex_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Make certain that the values are consistent with the defaults:
    initialised_ = false;
    restartOnRestart_ = false;
    restartOnOutput_ = false;
    periodicRestart_ = false;
    restartPeriod_ = GREAT;
    restartTime_ = GREAT;

    Info<< type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput",  restartOnOutput_);
    dict.readIfPresent("periodicRestart",  periodicRestart_);

    List<fieldAverageItem> faItems0;
    dict.readEntry("fields", faItems0);

    DynamicList<fieldAverageItem> faItems(faItems0.size());

    wordHashSet names;
    for (auto& item : faItems0)
    {
        if (names.insert(item.fieldName()))
        {
            item.setMeanFieldName(scopedName(item.meanFieldName()));
            item.setPrime2MeanFieldName(scopedName(item.prime2MeanFieldName()));
            faItems.push_back(item);
        }
        else
        {
            WarningInFunction
                << "Duplicate entry found: " << item.fieldName()
                << " (ignored)" << endl;
        }
    }

    faItems_.transfer(faItems);

    const scalar currentTime = obr().time().value();

    if (periodicRestart_)
    {
        scalar userRestartPeriod = dict.get<scalar>("restartPeriod");
        restartPeriod_ = obr().time().userTimeToTime(userRestartPeriod);

        if (restartPeriod_ > 0)
        {
            // Determine the appropriate interval for the next restart
            periodIndex_ = 1;
            while (currentTime > restartPeriod_*periodIndex_)
            {
                ++periodIndex_;
            }

            Info<< "    Restart period " << userRestartPeriod
                << " - next restart at " << (userRestartPeriod*periodIndex_)
                << nl << endl;
        }
        else
        {
            periodicRestart_ = false;

            Info<< "    Restart period " << userRestartPeriod
                << " - ignored"
                << nl << endl;
        }
    }

    scalar userRestartTime = 0;
    if (dict.readIfPresent("restartTime", userRestartTime))
    {
        restartTime_ = obr().time().userTimeToTime(userRestartTime);

        if (currentTime > restartTime_)
        {
            // The restart time is already in the past - ignore
            restartTime_ = GREAT;
        }
        else
        {
            Info<< "    Restart scheduled at time " << userRestartTime
                << nl << endl;
        }
    }

    readAveragingProperties();

    Info<< endl;

    return true;
}


bool Foam::functionObjects::fieldAverage::execute()
{
    calcAverages();

    return true;
}


bool Foam::functionObjects::fieldAverage::write()
{
    writeAverages();
    writeAveragingProperties();

    if (restartOnOutput_)
    {
        restart();
    }

    return true;
}


// ************************************************************************* //
