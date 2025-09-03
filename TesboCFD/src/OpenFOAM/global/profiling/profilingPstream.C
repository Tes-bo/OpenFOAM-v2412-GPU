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

\*---------------------------------------------------------------------------*/

#include "profilingPstream.H"
#include "List.H"
#include "Tuple2.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

std::unique_ptr<Foam::cpuTime> Foam::profilingPstream::timer_(nullptr);

bool Foam::profilingPstream::suspend_(false);

Foam::profilingPstream::timingList Foam::profilingPstream::times_(double(0));
Foam::profilingPstream::countList Foam::profilingPstream::counts_(uint64_t(0));


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::profilingPstream::enable()
{
    if (!timer_)
    {
        timer_.reset(new cpuTime);
        times_ = double(0);
        counts_ = uint64_t(0);
    }
    suspend_ = false;
}


void Foam::profilingPstream::disable() noexcept
{
    timer_.reset(nullptr);
    suspend_ = false;
}


void Foam::profilingPstream::reset()
{
    times_ = double(0);
    counts_ = uint64_t(0);
}


double Foam::profilingPstream::elapsedTime()
{
    double total = 0;
    for (const double val : times_)
    {
        total += val;
    }

    return total;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Loop over all values (with striding) and extract the value at given index
template<class Type>
inline static void extractValues
(
    UList<Type>& result,
    const int index,
    const UList<Type>& allValues
)
{
    if (result.empty())
    {
        return;
    }

    const label numProc = result.size();
    const Type* values = allValues.cbegin();
    const label stride = allValues.size() / numProc;

    if (!values || !stride)
    {
        result = Type(0);
        return;
    }

    for (label proci = 0; proci < numProc; ++proci, values += stride)
    {
        result[proci] = values[index];
    }
}


// Loop over all values (with striding) and extract combined value
// using the given unary function
template<class Type, class Extract>
inline static void extractValues
(
    UList<Type>& result,
    const UList<Type>& allValues,
    const Extract& extract
)
{
    if (result.empty())
    {
        return;
    }

    const label numProc = result.size();
    const Type* values = allValues.cbegin();
    const label stride = allValues.size() / numProc;

    if (!values || !stride)
    {
        result = Type(0);
        return;
    }

    for (label proci = 0; proci < numProc; ++proci, values += stride)
    {
        result[proci] = extract(values);
    }
}


inline static void printTimingDetail(const UList<double>& values)
{
    const label numProc = values.size();

    if (numProc)
    {
        Info<< indent << "    times   " << numProc << '(';

        for (label proci = 0; proci < numProc; ++proci)
        {
            if (proci) Info<< ' ';
            Info<< values[proci];
        }

        Info<< ')' << nl;
    }
}


inline static void printTimingDetail(const UList<uint64_t>& values)
{
    const label numProc = values.size();

    if (numProc)
    {
        // Output via std::ostream to avoid conversion to Foam::label
        // that Ostream performs

        auto& os = Info.stdStream();

        Info<< indent << "    counts  " << numProc << '(';

        for (label proci = 0; proci < numProc; ++proci)
        {
            if (proci) os << ' ';
            os << values[proci];
        }

        Info<< ')' << nl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingPstream::report(const int reportLevel)
{
    const label numProc = (UPstream::parRun() ? UPstream::nProcs() : 1);

    if (numProc < 2)
    {
        return;
    }

    // Use mpiGather on all values and perform the combinations
    // and statistics locally. This reduces the overall number of MPI
    // calls. For detailed output we need this information anyhow.

    // NB: profilingPstream uses a FixedList for timings(), counts()
    // so sizes are guaranteed to be consistent and identical everywhere.

    List<double> allTimes;
    List<uint64_t> allCounts;

    // Avoid disturbing any information
    const bool oldSuspend = suspend();

    {
        // The timings
        const auto& procValues = times_;

        if (UPstream::master())
        {
            allTimes.resize(numProc * procValues.size());
        }

        UPstream::mpiGather
        (
            procValues.cdata_bytes(),   // Send
            allTimes.data_bytes(),      // Recv
            procValues.size_bytes(),    // Num send/recv data per rank
            UPstream::commWorld()
        );
    }

    if (reportLevel > 1)
    {
        // The counts
        const auto& procValues = counts_;

        if (UPstream::master())
        {
            allCounts.resize(numProc * procValues.size());
        }

        UPstream::mpiGather
        (
            procValues.cdata_bytes(),   // Send
            allCounts.data_bytes(),     // Recv
            procValues.size_bytes(),    // Num send/recv data per rank
            UPstream::commWorld()
        );
    }

    // Resume if not previously suspended
    if (!oldSuspend)
    {
        resume();
    }


    // (Time, Processor) for each of: min/max/sum(avg)
    typedef FixedList<Tuple2<double, int>, 3> statData;

    // Extract min/max/average
    auto calcStats = [](const UList<double>& data) -> statData
    {
        statData stats;
        stats = Tuple2<double, int>((data.empty() ? 0 : data[0]), 0);

        const label np = data.size();
        for (label proci = 1; proci < np; ++proci)
        {
            Tuple2<double, int> tup(data[proci], proci);

            // 0: min, 1: max, 2: total(avg)
            if (stats[0].first() > tup.first()) stats[0] = tup;
            if (stats[1].first() < tup.first()) stats[1] = tup;
            stats[2].first() += tup.first();
        }

        // From total -> average value
        if (np) { stats[2].first() /= np; }

        return stats;
    };


    const auto printTimingStats =
        [&](Ostream& os, const char* tag, const statData& stats)
        {
            os  << indent << tag << ": avg = " << stats[2].first()
                << ", min = " << stats[0].first()
                << " (proc " << stats[0].second() << ')'
                << ", max = " << stats[1].first()
                << " (proc " << stats[1].second() << ')'
                << nl;
        };


    if (UPstream::master())
    {
        Info<< "profiling(parallel):" << nl
            << incrIndent;

        statData stats;
        List<double> extractedTimes(numProc);
        List<uint64_t> extractedCounts;

        if (reportLevel > 1)
        {
            extractedCounts.resize(numProc);
        }

        // Total times
        {
            extractValues
            (
                extractedTimes,
                allTimes,
                [=](const double values[])
                {
                    double total = 0;
                    for (unsigned i = 0; i < timingType::nCategories; ++i)
                    {
                        total += values[i];
                    }
                    return total;
                }
            );
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "total     ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
        }

        // all-all
        {
            const int index = int(timingType::ALL_TO_ALL);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "all-all   ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        // broadcast
        {
            const int index = int(timingType::BROADCAST);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "broadcast ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        // probe
        {
            const int index = int(timingType::PROBE);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "probe     ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        // Reduce/scatter times
        {
            // const int index = int(timingType::REDUCE);

            extractValues
            (
                extractedTimes,
                allTimes,
                [=](const double values[])
                {
                    return
                    (
                        values[timingType::REDUCE]
                      + values[timingType::GATHER]
                      + values[timingType::SCATTER]
                    );
                }
            );
            extractValues
            (
                extractedCounts,
                allCounts,
                [=](const uint64_t values[])
                {
                    return
                    (
                        values[timingType::REDUCE]
                      + values[timingType::GATHER]
                      + values[timingType::SCATTER]
                    );
                }
            );
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "reduce    ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        // Recv/send times
        #if 0  // FUTURE?
        {
            // const int index = int(timingType::RECV);

            extractValues
            (
                extractedTimes,
                allTimes,
                [=](const double values[])
                {
                    return
                    (
                        values[timingType::RECV]
                      + values[timingType::SEND]
                    );
                }
            );
            extractValues
            (
                extractedCounts,
                allCounts,
                [=](const uint64_t values[])
                {
                    return
                    (
                        values[timingType::RECV]
                      + values[timingType::SEND]
                    );
                }
            );
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "send/recv ", stats);
            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }
        #endif

        // request
        {
            const int index = int(timingType::REQUEST);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "request   ", stats);

            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        // wait
        {
            const int index = int(timingType::WAIT);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "wait      ", stats);

            if (reportLevel > 0) printTimingDetail(extractedTimes);
            if (reportLevel > 1) printTimingDetail(extractedCounts);
        }

        Info<< decrIndent;
    }
}


// ************************************************************************* //
