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

#include "timeControl.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::timeControl::timeControls
>
Foam::timeControl::controlNames_
({
    { timeControl::ocNone, "none" },
    { timeControl::ocAlways, "always" },
    { timeControl::ocTimeStep, "timeStep" },
    { timeControl::ocWriteTime, "writeTime" },
    { timeControl::ocWriteTime, "outputTime" },
    { timeControl::ocRunTime, "runTime" },
    { timeControl::ocAdjustableRunTime, "adjustable" },
    { timeControl::ocAdjustableRunTime, "adjustableRunTime" },
    { timeControl::ocClockTime, "clockTime" },
    { timeControl::ocCpuTime, "cpuTime" },
    { timeControl::ocOnEnd, "onEnd" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeControl::timeControl
(
    const Time& runTime,
    const word& prefix
)
:
    time_(runTime),
    prefix_(prefix),
    timeControl_(ocAlways),
    intInterval_(0),
    interval_(0),
    executionIndex_(0)
{}


Foam::timeControl::timeControl
(
    const Time& runTime,
    const dictionary& dict,
    const word& prefix
)
:
    timeControl(runTime, prefix)
{
    read(dict);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::timeControl::entriesPresent
(
    const dictionary& dict,
    const word& prefix
)
{
    const word controlName(prefix + "Control");

    return dict.found(controlName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeControl::clear()
{
    timeControl_ = ocAlways;
    intInterval_ = 0;
    interval_ = 0;
    executionIndex_ = 0;
}


void Foam::timeControl::read(const dictionary& dict)
{
    // Default is timeStep
    timeControl_ = ocTimeStep;
    intInterval_ = 0;
    interval_ = 0;

    word controlName(prefix_ + "Control");
    word intervalName(prefix_ + "Interval");

    if (prefix_ == "write")
    {
        // TBD: Could have   timeControl_ = ocWriteTime;

        if (dict.found("outputControl"))
        {
            // Accept deprecated 'outputControl' instead of 'writeControl'

            // Change to the old names for this option
            controlName = "outputControl";
            intervalName = "outputInterval";

            IOWarningInFunction(dict)
                << "Found deprecated 'outputControl'" << nl
                << "    Use 'writeControl' with 'writeInterval'"
                << endl;
            error::warnAboutAge("keyword", 1606);
        }
    }


    timeControl_ = controlNames_.getOrDefault(controlName, dict, timeControl_);

    switch (timeControl_)
    {
        case ocTimeStep:
        case ocWriteTime:
        {
            intInterval_ = dict.getOrDefault<label>(intervalName, 0);
            interval_ = intInterval_;  // Mirrored as scalar (for output)
            break;
        }

        case ocClockTime:
        case ocRunTime:
        case ocCpuTime:
        case ocAdjustableRunTime:
        {
            const scalar userTime = dict.get<scalar>(intervalName);
            interval_ = time_.userTimeToTime(userTime);
            break;
        }

        default:
        {
            break;
        }
    }
}


bool Foam::timeControl::execute()
{
    switch (timeControl_)
    {
        case ocNone:
        {
            return false;
            break;
        }

        case ocAlways:
        {
            return true;
            break;
        }

        case ocTimeStep:
        {
            return
            (
                (intInterval_ <= 1)
             || !(time_.timeIndex() % intInterval_)
            );
            break;
        }

        case ocWriteTime:
        {
            if (time_.writeTime())
            {
                ++executionIndex_;
                return
                (
                    (intInterval_ <= 1)
                 || !(executionIndex_ % intInterval_)
                );
            }
            break;
        }

        case ocRunTime:
        case ocAdjustableRunTime:
        {
            label executionIndex = label
            (
                (
                    (time_.value() - time_.startTime().value())
                  + 0.5*time_.deltaTValue()
                )
               /interval_
            );

            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocCpuTime:
        {
            label executionIndex = label
            (
                returnReduce(time_.elapsedCpuTime(), maxOp<double>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocClockTime:
        {
            label executionIndex = label
            (
                returnReduce(time_.elapsedClockTime(), maxOp<double>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case ocOnEnd:
        {
            scalar endTime = time_.endTime().value() - 0.5*time_.deltaTValue();
            return time_.value() > endTime;
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined time control: "
                << controlNames_[timeControl_] << nl
                << abort(FatalError);
            break;
        }
    }

    return false;
}


// ************************************************************************* //
