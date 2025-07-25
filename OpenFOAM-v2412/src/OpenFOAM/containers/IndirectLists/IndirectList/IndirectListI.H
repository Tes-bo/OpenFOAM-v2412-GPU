/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    const labelUList& addr
)
:
    IndirectListAddressing<labelList>(addr),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    labelList&& addr
)
:
    IndirectListAddressing<labelList>(std::move(addr)),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    const Foam::zero
)
:
    IndirectListAddressing<labelList>(labelList()),  // zero-size addressing
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(const IndirectList<T>& list)
:
    // Copy addressing
    IndirectListAddressing<labelList>(list.addressing()),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(IndirectList<T>&& list)
:
    // Move addressing
    IndirectListAddressing<labelList>(std::move(list.addressing())),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(const UIndirectList<T>& list)
:
    // Copy addressing
    IndirectListAddressing<labelList>(list.addressing()),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
template<class UnaryCondition>
Foam::IndirectList<T> Foam::IndirectList<T>::subset
(
    const UList<T>& values,
    const UnaryCondition& select,
    const bool invert
)
{
    const label len = values.size();

    IndirectList<T> result(values, Foam::zero{});
    labelList& addr = result.addressing();

    addr.resize_nocopy(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        // Test on position
        if (select(i) ? !invert : invert)
        {
            addr[count] = i;
            ++count;
        }
    }
    addr.resize(count);

    return result;
}


template<class T>
template<class UnaryPredicate>
Foam::IndirectList<T> Foam::IndirectList<T>::subset_if
(
    const UList<T>& values,
    const UnaryPredicate& pred,
    const bool invert
)
{
    const label len = values.size();

    IndirectList<T> result(values, Foam::zero{});
    labelList& addr = result.addressing();

    addr.resize_nocopy(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        // Test on value
        if (pred(values[i]) ? !invert : invert)
        {
            addr[count] = i;
            ++count;
        }
    }
    addr.resize(count);

    return result;
}


template<class T>
Foam::IndirectList<T> Foam::IndirectList<T>::uniq
(
    const UList<T>& values,
    const bool sorted
)
{
    const label len = values.size();

    IndirectList<T> result(values, Foam::zero{});
    labelList& order = result.addressing();

    // Start from sorted order
    Foam::sortedOrder(values, order);

    if (len > 1)
    {
        label count = 0;
        for (label i = 1; i < len; ++i)
        {
            // Eliminate duplicates
            if (values[order[count]] != values[order[i]])
            {
                ++count;
                order[count] = order[i];
            }
        }
        ++count;
        order.resize(count);

        // Recover the original input order
        if (!sorted)
        {
            Foam::sort(order);
        }
    }

    return result;
}


// ************************************************************************* //
