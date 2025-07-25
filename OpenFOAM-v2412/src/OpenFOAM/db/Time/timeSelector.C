/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "timeSelector.H"
#include "ListOps.H"
#include "argList.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeSelector::timeSelector(const std::string& str)
:
    ranges_(str)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeSelector::contains(const scalar value) const
{
    return ranges_.contains(value);
}


bool Foam::timeSelector::contains(const instant& t) const
{
    return ranges_.contains(t.value());
}


Foam::List<bool> Foam::timeSelector::selected(const instantList& times) const
{
    List<bool> selectTimes(times.size());

    // Check ranges, avoid false positive on constant/
    forAll(times, timei)
    {
        selectTimes[timei] =
        (
            times[timei].name() != "constant" && contains(times[timei])
        );
    }

    // Check specific values
    for (const scalarRange& range : ranges_)
    {
        if (range.single())
        {
            const scalar target = range.value();

            const label nearestIndex =
                TimePaths::findClosestTimeIndex(times, target);

            // Note could also test if the index is too far away.
            // Eg, for times (0 10 20 30 40) selecting 100 will currently
            // return the closest time (40), but perhaps we should limit that
            // to the last deltaT?

            if (nearestIndex >= 0)
            {
                selectTimes[nearestIndex] = true;
            }
        }
    }

    return selectTimes;
}


Foam::instantList Foam::timeSelector::select(const instantList& times) const
{
    return subset(selected(times), times);
}


void Foam::timeSelector::inplaceSelect(instantList& times) const
{
    inplaceSubset(selected(times), times);
}


void Foam::timeSelector::addOptions
(
    const bool constant,
    const bool withZero
)
{
    if (constant)
    {
        argList::addBoolOption
        (
            "constant",
            "Include 'constant/' dir in the times list"
        );
    }
    if (withZero)
    {
        argList::addBoolOption
        (
            "withZero",
            "Include '0/' dir in the times list"
        );
    }
    argList::addBoolOption
    (
        "noZero",
        string("Exclude '0/' dir from the times list")
      + (
            withZero
          ? ", has precedence over the -withZero option"
          : ""
        )
    );
    argList::addBoolOption
    (
        "latestTime",
        "Select the latest time"
    );
    argList::addOption
    (
        "time",
        "ranges",
        "List of ranges. Eg, ':10,20 40:70 1000:', 'none', etc"
    );
}


void Foam::timeSelector::addOptions_singleTime()
{
    argList::addBoolOption
    (
        "constant",
        "Include 'constant/' dir in the times"
    );
    argList::addBoolOption
    (
        "noZero",
        "Exclude '0/' dir from the times (currently ignored)"
    );
    argList::addBoolOption
    (
        "latestTime",
        "Select the latest time"
    );
    argList::addOption
    (
        "time",
        "value",
        "Select the nearest time to the specified value"
    );
}


Foam::instantList Foam::timeSelector::select
(
    const instantList& times,
    const argList& args,
    const word& constantName
)
{
    if (times.size())
    {
        List<bool> selectTimes(times.size(), true);

        label constantIdx = -1;
        label zeroIdx = -1;
        label latestIdx = -1;

        // Determine locations of constant/ and 0/ directories
        forAll(times, timei)
        {
            if (times[timei].name() == constantName)
            {
                constantIdx = timei;
            }
            else if (times[timei].value() == 0)
            {
                zeroIdx = timei;
            }

            if (constantIdx >= 0 && zeroIdx >= 0)
            {
                break;
            }
        }

        // Determine latestTime selection (if any)
        // This must appear before the -time option processing
        if (args.found("latestTime"))
        {
            selectTimes = false;
            latestIdx = times.size() - 1;

            // Avoid false match on constant/
            if (latestIdx == constantIdx)
            {
                latestIdx = -1;
            }
        }

        if (args.found("time"))
        {
            // Can match 0/, but can never match constant/
            selectTimes = timeSelector(args["time"]).selected(times);
        }

        // Add in latestTime (if selected)
        if (latestIdx >= 0)
        {
            selectTimes[latestIdx] = true;
        }

        if (constantIdx >= 0)
        {
            // Only add constant/ if specifically requested
            selectTimes[constantIdx] = args.found("constant");
        }

        // Special treatment for 0/
        if (zeroIdx >= 0)
        {
            if (args.found("noZero"))
            {
                // Exclude 0/ if specifically requested
                selectTimes[zeroIdx] = false;
            }
            else if (argList::validOptions.found("withZero"))
            {
                // With -withZero enabled, drop 0/ unless specifically requested
                selectTimes[zeroIdx] = args.found("withZero");
            }
        }

        return subset(selectTimes, times);
    }

    return times;
}


Foam::instantList Foam::timeSelector::select0
(
    Time& runTime,
    const argList& args
)
{
    instantList times
    (
        timeSelector::select
        (
            runTime.times(),
            args,
            runTime.constant()
        )
    );

    if (times.empty())
    {
        WarningInFunction
            << "No time specified or available, selecting 'constant'"
            << endl;

        times.push_back(instant(0, runTime.constant()));
    }

    runTime.setTime(times.front(), 0);

    return times;
}


Foam::instantList Foam::timeSelector::selectIfPresent
(
    Time& runTime,
    const argList& args
)
{
    if
    (
        args.found("latestTime")
     || args.found("time")
     || args.found("constant")
     || args.found("noZero")
     || args.found("withZero")
    )
    {
        return select0(runTime, args);
    }

    // No timeSelector option specified. Do not change runTime.
    return instantList(one{}, instant(runTime.value(), runTime.timeName()));
}


bool Foam::timeSelector::setTimeIfPresent
(
    Time& runTime,
    const argList& args,
    const bool forceInitial
)
{
    label timei = -1;
    instantList times;

    if
    (
        forceInitial
     || args.found("constant")
     || args.found("latestTime")
     || args.found("time")
        // Currently ignoring -noZero, -withZero
    )
    {
        // Get times list
        times = runTime.times();
    }

    if (times.size())
    {
        // Start from first time (eg, for -constant or forced)
        timei = 0;

        // Determine latestTime selection (if any)
        // This must appear before the -time option processing
        if (args.found("latestTime"))
        {
            timei = times.size() - 1;
        }
        else if (args.found("time"))
        {
            const scalar target = args.get<scalar>("time");

            timei = TimePaths::findClosestTimeIndex(times, target);
        }


        // Avoid "constant" unless specifically requested with -constant,
        // and the -constant option is actually an expected option

        if
        (
            (timei >= 0 && timei < times.size()-1)
         && times[timei].name() == "constant"
         && (argList::validOptions.found("constant") && !args.found("constant"))
        )
        {
            ++timei;
        }
    }


    if (timei >= 0 && timei < times.size())
    {
        // Specified a timeSelector option, or forceInitial.
        // Set the runTime accordingly.

        runTime.setTime(times[timei], timei);
        return true;
    }
    else
    {
        // No timeSelector option specified. Do not change runTime.
        return false;
    }
}


// ************************************************************************* //
