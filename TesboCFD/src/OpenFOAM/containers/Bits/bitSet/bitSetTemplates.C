/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "FixedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned N>
Foam::bitSet::bitSet(const label n, const FixedList<label, N>& locations)
:
    bitSet(n)
{

    setMany(locations.begin(), locations.end());
}


template<unsigned N>
Foam::bitSet::bitSet(const FixedList<label, N>& locations)
:
    bitSet()
{

    setMany(locations.begin(), locations.end());
}


template<class Addr>
Foam::bitSet::bitSet
(
    const bitSet& bitset,
    const IndirectListBase<label, Addr>& addr
)
:
    bitSet(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, bitset.get(addr[i]));
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class InputIter>
Foam::label Foam::bitSet::setMany(InputIter first, InputIter last)
{
    // Check the max expected value first
    const auto max = std::max_element(first, last);
    const label len = (max != last ? (1 + *max) : 0);

    label changed = 0;

    if (len > 0)
    {
        reserve(len);

        for (; first != last; ++first)
        {
            if (set(*first))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class InputIter>
Foam::label Foam::bitSet::unset(InputIter first, InputIter last)
{
    label changed = 0;

    for (; first != last; ++first)
    {
        if (unset(*first))
        {
            ++changed;
        }
    }

    return changed;
}


template<unsigned N>
Foam::label Foam::bitSet::set(const FixedList<label, N>& locations)
{
    return setMany(locations.begin(), locations.end());
}


template<unsigned N>
Foam::label Foam::bitSet::unset(const FixedList<label, N>& locations)
{
    return unset(locations.begin(), locations.end());
}


// ************************************************************************* //
