/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "BitOps.H"
#include "bitSet.H"
#include "HashSet.H"
#include "List.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * * BitOps  * * * * * * * * * * * * * * * * //

// See bitSet::setMany for original implementation
void Foam::BitOps::set(List<bool>& bools, const labelUList& locations)
{
    // Check the max expected value first
    const auto max = std::max_element(locations.begin(), locations.end());
    const label len = (max != locations.end() ? (1 + *max) : 0);

    if (len > bools.size())
    {
        bools.resize(len, false);
    }

    for (label i : locations)
    {
        if (i >= 0)
        {
            bools[i] = true;
        }
    }
}


// See bitSet::set(labelRange) for original implementation
void Foam::BitOps::set(List<bool>& bools, const labelRange& range)
{
    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    // Range is invalid (zero-sized or entirely negative) - noop
    if (slice.empty())
    {
        return;
    }

    // Check maximum extent of the range.
    // The end_value() method is the exclusive end-value,
    // which corresponds to our potential new length.
    // - resize now to avoid allocations within the loop

    if (slice.end_value() >= bools.size())
    {
        bools.resize(slice.end_value(), false);
    }

    for (const label i : slice)
    {
        bools.set(i);
    }
}


// See bitSet::set(labelRange) for original implementation
void Foam::BitOps::set(labelHashSet& hashset, const labelRange& range)
{
    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    for (const label i : slice)
    {
        hashset.set(i);
    }
}


void Foam::BitOps::set(bitSet& bitset, const labelRange& range)
{
    bitset.set(range);
}


void Foam::BitOps::unset(List<bool>& bools, const labelUList& locations)
{
    for (const label i : locations)
    {
        bools.unset(i);
    }
}


// See bitSet::unset(labelRange) for original implementation
void Foam::BitOps::unset(List<bool>& bools, const labelRange& range)
{
    for (const label i : range)
    {
        bools.unset(i);
    }
}


void Foam::BitOps::unset(labelHashSet& hashset, const labelRange& range)
{
    for (const label i : range)
    {
        hashset.unset(i);
    }
}


void Foam::BitOps::unset(bitSet& bitset, const labelRange& range)
{
    bitset.unset(range);
}


Foam::List<bool> Foam::BitOps::select
(
    const label n,
    const labelUList& locations
)
{
    List<bool> bools(n, false);

    BitOps::set(bools, locations);

    return bools;
}


Foam::List<bool> Foam::BitOps::select(const labelUList& locations)
{
    List<bool> bools;

    BitOps::set(bools, locations);

    return bools;
}


// Note: code is like ListOps findIndices() and/or bitSet toc()
Foam::List<Foam::label> Foam::BitOps::toc(const UList<bool>& bools)
{
    const label len = bools.size();

    // Pass 1: count occurrences
    label count = 0;

    for (const bool b : bools)
    {
        if (b) ++count;
    }

    labelList indices(count);

    // Pass 2: fill content
    if (count)
    {
        const label total(count);
        count = 0;

        for (label i = 0; i < len; ++i)
        {
            if (bools[i])
            {
                indices[count] = i;
                if (++count == total)  // Terminate early
                {
                    break;
                }
            }
        }
    }

    return indices;
}


Foam::List<Foam::label> Foam::BitOps::sortedToc(const UList<bool>& bools)
{
    return BitOps::toc(bools);
}


Foam::List<Foam::label> Foam::BitOps::toc(const bitSet& bitset)
{
    return bitset.toc();
}


Foam::List<Foam::label> Foam::BitOps::sortedToc(const bitSet& bitset)
{
    return bitset.sortedToc();
}


Foam::List<Foam::label> Foam::BitOps::toc(const labelHashSet& hashset)
{
    return hashset.sortedToc();
}


Foam::List<Foam::label> Foam::BitOps::sortedToc(const labelHashSet& hashset)
{
    return hashset.sortedToc();
}


// * * * * * * * * * * * * * * * * BitSetOps * * * * * * * * * * * * * * * * //

Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const labelHashSet& locations,
    const bool on
)
{
    bitSet output(n, !on);

    for (const label idx : locations)
    {
        // Restrict the input size
        if (idx < n)
        {
            output.set(idx, on);
        }
    }

    return output;
}


Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const labelUList& locations,
    const bool on
)
{
    bitSet output(n, !on);

    for (const label idx : locations)
    {
        // Restrict the input size
        if (idx < n)
        {
            output.set(idx, on);
        }
    }

    return output;
}


Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const label select,
    const labelUList& values,
    const bool on
)
{
    bitSet output(n, !on);

    // Restrict the input size
    const label len = std::min(n, values.size());

    for (label idx = 0; idx < len; ++idx)
    {
        if (select == values[idx])
        {
            output.set(idx, on);
        }
    }

    return output;
}


// ************************************************************************* //
