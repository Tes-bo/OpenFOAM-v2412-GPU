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

#include "PackedList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelUList& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(addr[i]));
    }
}


template<unsigned Width>
template<class Addr>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const IndirectListBase<label, Addr>& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(addr[i]));
    }
}


template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelRange& range
)
:
    PackedList<Width>(range.size())
{
    label pos = range.start();
    const label len = range.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(pos));
        ++pos;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
bool Foam::PackedList<Width>::uniform() const
{
    // Trivial cases
    if (empty())
    {
        return false;
    }
    else if (size() == 1)
    {
        return true;
    }

    // The value of the first element for testing
    const unsigned int val = get(0);

    bool identical = true;

    if (!val)
    {
        // No bits set: just check there are no non-zero blocks
        // - like bitSet::none()
        identical = (-1 == first_block());
    }
    else if (val == PackedList<Width>::max_value)
    {
        // All bits set: just check there are no zero blocks
        // - like bitSet::all()
        identical = (-1 == first_not_block());
    }
    else
    {
        const label nblocks = num_blocks(size());

        if (nblocks > 1)
        {
            // Fill value for complete blocks
            const unsigned int blockval = repeated_value(val);

            // Check each complete block (nblocks-1)
            for (label blocki = 0; identical && blocki < (nblocks-1); ++blocki)
            {
                identical = (blocks_[blocki] == blockval);
            }
        }

        // Partial block: check manually
        for
        (
            label elemi = elem_per_block*(nblocks-1);
            identical && elemi < size();
            ++elemi
        )
        {
            identical = (val == get(elemi));
        }
    }

    return identical;
}


template<unsigned Width>
bool Foam::PackedList<Width>::equal(const PackedList<Width>& other) const
{
    if (size() != other.size())
    {
        return false;
    }

    const label nblocks = num_blocks(size());
    const auto& rhs = other.blocks_;

    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        if (blocks_[blocki] != rhs[blocki])
        {
            return false;
        }
    }

    return true;
}


template<unsigned Width>
Foam::labelList Foam::PackedList<Width>::values() const
{
    return this->unpack<label>();
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack() const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold result"
    );

    List<IntType> output(size());

    if (empty())
    {
        return output;
    }
    else if (uniform())
    {
        output = static_cast<IntType>(get(0));
        return output;
    }

    // NON-UNIFORM and len > 0

    label outi = 0;

    // Process n-1 complete blocks
    const label nblocks = num_blocks(size());

    for (label blocki=0; blocki < nblocks-1; ++blocki)
    {
        unsigned int blockval = blocks_[blocki];

        for (unsigned nget = elem_per_block; nget; --nget, ++outi)
        {
            output[outi] = IntType(blockval & PackedList<Width>::max_value);
            blockval >>= Width;
        }
    }

    // Any partial blocks
    for (/*nil*/; outi < size(); ++outi)
    {
        output[outi] = get(outi);
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelRange& range) const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for unpack output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold unpack output."
    );


    // Could be more efficient but messier with block-wise access.
    // - automatically handles any invalid positions

    auto pos = range.start();

    List<IntType> output(range.size());

    for (IntType& out : output)
    {
        out = IntType(get(pos));
        ++pos;
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelUList& locations) const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for unpack output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold unpack output."
    );


    label pos = 0;

    List<IntType> output(locations.size());

    for (IntType& out : output)
    {
        out = IntType(get(locations[pos]));
        ++pos;
    }

    return output;
}


// ************************************************************************* //
