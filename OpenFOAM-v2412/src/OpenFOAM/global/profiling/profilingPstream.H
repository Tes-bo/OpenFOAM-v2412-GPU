/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::profilingPstream

Description
    Timers and values for simple (simplistic) mpi-profiling.
    The entire class behaves as a singleton.

SourceFiles
    profilingPstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_profilingPstream_H
#define Foam_profilingPstream_H

#include "cpuTime.H"
#include "FixedList.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class profilingPstream Declaration
\*---------------------------------------------------------------------------*/

class profilingPstream
{
public:

    // Public Types

        //- The enumerated timing categories (for times and counts arrays)
        enum timingType : unsigned
        {
            ALL_TO_ALL = 0,
            BROADCAST,
            PROBE,
            REDUCE,
            GATHER,         // gather (or recv)
            SCATTER,        // scatter (or send)
            REQUEST,
            WAIT,
            OTHER,
            nCategories     // Dimensioning size
        };

        //- Fixed-size container for timing values
        typedef FixedList<double, timingType::nCategories> timingList;

        //- Fixed-size container for timing counts
        typedef FixedList<uint64_t, timingType::nCategories> countList;


private:

    // Private Static Data

        //- The timer to use
        static std::unique_ptr<cpuTime> timer_;

        //- Is timer in a suspend state?
        static bool suspend_;

        //- The accumulated values for various timing categories
        static timingList times_;

        //- The timing frequency for various timing categories
        static countList counts_;


public:

    // Static Member Functions

    // Management

        //- True if timer is active (ie, enabled and not suspended)
        static bool active() noexcept { return !suspend_ && timer_; }

        //- Create timer for measuring communication or un-suspend existing
        static void enable();

        //- Remove timer for measuring communication activity.
        //- Does not affect times/counts.
        static void disable() noexcept;

        //- Reset times/counts. Does not affect the timer itself
        static void reset();

        //- Suspend use of timer. Return old status
        static bool suspend() noexcept
        {
            bool old(suspend_);
            suspend_ = bool(timer_);
            return old;
        }

        //- Resume use of timer (if previously active)
        static void resume() noexcept
        {
            suspend_ = false;
        }


    // Timing/Counts

        //- Access to the timing information
        static timingList& times() noexcept { return times_; }

        //- Access to the timing counts
        static countList& counts() noexcept { return counts_; }

        //- Access to the timing information for given timing category
        static double times(const timingType idx)
        {
            return times_[idx];
        }

        //- Access to the count for given timing category
        static uint64_t counts(const timingType idx)
        {
            return counts_[idx];
        }

        //- The total of times
        static double elapsedTime();

        //- Update timer prior to measurement
        static void beginTiming()
        {
            if (!suspend_ && timer_)
            {
                timer_->resetCpuTimeIncrement();
            }
        }

        //- Add time increment
        static void addTime(const timingType idx)
        {
            if (!suspend_ && timer_)
            {
                times_[idx] += timer_->cpuTimeIncrement();
                ++counts_[idx];
            }
        }

        //- Add time increment to \em broadcast time
        static void addBroadcastTime()
        {
            addTime(timingType::BROADCAST);
        }

        //- Add time increment to \em reduce time
        static void addReduceTime()
        {
            addTime(timingType::REDUCE);
        }

        //- Add time increment to \em probe time
        static void addProbeTime()
        {
            addTime(timingType::PROBE);
        }

        //- Add time increment to \em request time
        static void addRequestTime()
        {
            addTime(timingType::REQUEST);
        }

        //- Add time increment to \em wait time
        static void addWaitTime()
        {
            addTime(timingType::WAIT);
        }

        //- Add time increment to \em gather time
        static void addGatherTime()
        {
            addTime(timingType::GATHER);
        }

        //- Add time increment to \em scatter time
        static void addScatterTime()
        {
            addTime(timingType::SCATTER);
        }

        //- Add time increment to \em allToAll time
        static void addAllToAllTime()
        {
            addTime(timingType::ALL_TO_ALL);
        }

        //- Add time increment to \em other time
        static void addOtherTime()
        {
            addTime(timingType::OTHER);
        }


    // Output

        //- Report current information. Uses parallel communication!
        static void report(const int reportLevel = 0);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
