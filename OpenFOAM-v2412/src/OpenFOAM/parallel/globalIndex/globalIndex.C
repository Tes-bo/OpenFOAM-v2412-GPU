/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "globalIndex.H"
#include "Pair.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::globalIndex::reportOverflowAndExit
(
    const label idx,
    const label prevOffset,
    const label count
)
{
    if (idx < 0)
    {
        // No overflow tagged
        return;
    }

    FatalErrorInFunction
        << "Overflow : sum of sizes exceeds labelMax ("
        << labelMax << ") after index " << idx;

    if (prevOffset >= 0 && count >= 0)
    {
        FatalError
            << " while trying to add (" << count
            << ") to offset (" << prevOffset << ")";
    }

    FatalError
        << nl
        << "Please recompile with larger datatype for label." << nl
        << exit(FatalError);
}


Foam::labelRange
Foam::globalIndex::calcRange
(
    const label localSize,
    const label comm,
    const bool checkOverflow
)
{
    // Range with 0-offset initially
    labelRange myRange(0, localSize);

    if (!UPstream::is_parallel(comm))
    {
        return myRange;
    }

    const label myProci = UPstream::myProcNo(comm);
    const labelList counts = UPstream::allGatherValues(localSize, comm);

    if (checkOverflow)
    {
        const label len = counts.size();

        label start = 0;

        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];

            if (i == myProci)
            {
                myRange.start() = start;
            }

            const label prev = start;
            start += count;

            if (start < prev)
            {
                reportOverflowAndExit(i, prev, count);
            }
        }
    }
    else
    {
        // std::accumulate
        // (
        //     counts.cbegin(),
        //     counts.cbegin(myProci),
        //     label(0)
        // );

        label start = 0;

        for (label i = 0; i < myProci; ++i)
        {
            start += counts[i];
        }
        myRange.start() = start;
    }

    return myRange;
}


Foam::label
Foam::globalIndex::calcOffset
(
    const label localSize,
    const label comm,
    const bool checkOverflow
)
{
    // Placeholder value
    label myOffset = 0;

    if (!UPstream::is_parallel(comm))
    {
        return myOffset;
    }

    const label myProci = UPstream::myProcNo(comm);
    const labelList counts = UPstream::allGatherValues(localSize, comm);

    if (checkOverflow)
    {
        const label len = counts.size();

        label start = 0;

        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];
            if (i == myProci)
            {
                myOffset = start;
            }

            const label prev = start;
            start += count;

            if (start < prev)
            {
                reportOverflowAndExit(i, prev, count);
            }
        }
    }
    else
    {
        // std::accumulate
        // (
        //     counts.cbegin(),
        //     counts.cbegin(myProci),
        //     label(0)
        // );

        label start = 0;

        for (label i = 0; i < myProci; ++i)
        {
            start += counts[i];
        }
        myOffset = start;
    }

    return myOffset;
}


Foam::labelList
Foam::globalIndex::calcOffsets
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    labelList values;

    const label len = counts.size();

    if (len)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];
            values[i] = start;
            start += count;

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, values[i], count);
            }
        }
        values[len] = start;
    }

    return values;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::calcRanges
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    List<labelRange> values;

    const label len = counts.size();

    if (len)
    {
        values.resize(len);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];
            values[i].reset(start, count);
            start += count;

            if
            (
                checkOverflow
             && (start < values[i].start())
             && (i < len-1)  // Do not check the one beyond the end range
            )
            {
                reportOverflowAndExit(i, values[i].start(), count);
            }
        }
    }

    return values;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex(Istream& is)
{
    is >> offsets_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::CompactListList<Foam::label>
Foam::globalIndex::bin
(
    const labelUList& offsets,
    const labelUList& globalIds,
    labelList& order,
    DynamicList<label>& validBins
)
{
    Foam::sortedOrder(globalIds, order);
    validBins.clear();

    CompactListList<label> bins;

    if (!globalIds.empty())
    {
        labelList& binOffsets = bins.offsets();
        binOffsets.resize(offsets.size(), Foam::zero{});

        labelList& binValues = bins.values();
        binValues = UIndirectList<label>(globalIds, order);

        const label id = binValues[0];
        label proci = findLower(offsets, id+1);

        validBins.push_back(proci);
        label binSize = 1;

        for (label i = 1; i < order.size(); i++)
        {
            const label id = binValues[i];

            if (id < offsets[proci+1])
            {
                ++binSize;
            }
            else
            {
                // Not local. Reset proci
                const label oldProci = proci;
                proci = findLower(offsets, id+1);

                // Set offsets
                for (label j = oldProci+1; j < proci; ++j)
                {
                    binOffsets[j] = binOffsets[oldProci]+binSize;
                }
                binOffsets[proci] = i;
                validBins.push_back(proci);
                binSize = 1;
            }
        }

        for (label j = proci+1; j < binOffsets.size(); ++j)
        {
            binOffsets[j] = binOffsets[proci]+binSize;
        }
    }

    return bins;
}


void Foam::globalIndex::resize(const label n)
{
    if (n < 1)
    {
        offsets_.clear();
    }
    else
    {
        offsets_.resize(n+1, end_value());
    }
}


void Foam::globalIndex::reset
(
    const label localSize,
    const label comm,
    const bool parallel
)
{
    const label len = UPstream::nProcs(comm);

    if (len)
    {
        labelList counts;

        if (parallel && UPstream::parRun())  // or UPstream::is_parallel(comm)
        {
            counts = UPstream::allGatherValues(localSize, comm);
        }
        else
        {
            // Non-parallel branch: use localSize on-proc, zero elsewhere
            // TBD: check for (proci >= 0) ?
            const auto proci = UPstream::myProcNo(comm);

            counts.resize(len, Foam::zero{});
            counts[proci] = localSize;
        }

        reset(counts, true);  // checkOverflow = true
    }
    else
    {
        // Nothing to do
        offsets_.clear();
    }
}


void Foam::globalIndex::reset
(
    const labelUList& counts,
    const bool checkOverflow
)
{
    const label len = counts.size();

    if (len)
    {
        offsets_.resize_nocopy(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];
            offsets_[i] = start;
            start += count;

            if (checkOverflow && start < offsets_[i])
            {
                reportOverflowAndExit(i, offsets_[i], count);
            }
        }
        offsets_[len] = start;
    }
    else
    {
        offsets_.clear();
    }
}


void Foam::globalIndex::setLocalSize(const label proci, const label len)
{
    if (proci >= 0 && proci+1 < offsets_.size() && len >= 0)
    {
        const label delta = (len - (offsets_[proci+1] - offsets_[proci]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = proci+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


Foam::labelList Foam::globalIndex::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci] = offsets_[proci+1] - offsets_[proci];
    }

    return values;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::ranges() const
{
    List<labelRange> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci].reset
        (
            offsets_[proci],
            (offsets_[proci+1] - offsets_[proci])
        );
    }

    return values;
}


Foam::label Foam::globalIndex::maxNonLocalSize(const label proci) const
{
    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return 0;
    }

    label maxLen = 0;

    for (label i=0; i < len; ++i)
    {
        if (i != proci)
        {
            const label count = (offsets_[i+1] - offsets_[i]);
            maxLen = max(maxLen, count);
        }
    }

    return maxLen;
}


Foam::labelRange Foam::globalIndex::front() const
{
    return
    (
        (offsets_.size() < 2)
      ? labelRange()
      : labelRange(Pair<label>(offsets_[0], offsets_[1]))
    );
}


Foam::labelRange Foam::globalIndex::back() const
{
    return
    (
        (offsets_.size() < 2)
      ? labelRange()
      : labelRange
        (
            Pair<label>
            (
                offsets_[offsets_.size()-2],
                offsets_[offsets_.size()-1]
            )
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, globalIndex& gi)
{
    return is >> gi.offsets();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const globalIndex& gi)
{
    return os << gi.offsets();
}


// ************************************************************************* //
