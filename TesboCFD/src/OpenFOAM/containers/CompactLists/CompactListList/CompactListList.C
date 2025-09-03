/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "CompactListList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
void Foam::CompactListList<T>::reportOverflowAndExit
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


template<class T>
template<class ListListType>
Foam::CompactListList<T> Foam::CompactListList<T>::pack_impl
(
    const ListListType& lists,
    const bool checkOverflow
)
{
    CompactListList<T> compact;

    auto& newOffsets = compact.offsets_;
    auto& newValues = compact.values_;

    label total = 0;
    const label len = lists.size();

    if (len)
    {
        newOffsets.resize(len+1, Foam::zero{});

        for (label i = 0; i < len; ++i)
        {
            const label count = lists[i].size();

            newOffsets[i] = total;
            total += count;

            if (checkOverflow && total < newOffsets[i])
            {
                reportOverflowAndExit(i, newOffsets[i], count);
            }
        }
        newOffsets[len] = total;
    }

    if (total)
    {
        // Make a deepCopy of data
        newValues.resize(total);

        auto iter = newValues.begin();

        // NB: operator[] for sub-list read access (eg, an indirect list)
        // cannot replace with std::copy

        for (const auto& list : lists)
        {
            const label sublen = list.size();

            for (label i = 0; i < sublen; (void)++i, (void)++iter)
            {
                *iter = list[i];
            }
        }
    }

    return compact;
}


template<class T>
template<class SubListType>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const UList<SubListType>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::pack_impl<UList<SubListType>>
    (
        lists,
        checkOverflow
    );
}


template<class T>
template<class SubListType, class Addr>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const IndirectListBase<SubListType, Addr>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::pack_impl<IndirectListBase<SubListType, Addr>>
    (
        lists,
        checkOverflow
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
Foam::label Foam::CompactListList<T>::resize_offsets
(
    const labelUList& listSizes,
    const bool checkOverflow
)
{
    const label len = listSizes.size();
    label total = 0;

    if (len)
    {
        offsets_.resize(len+1);

        for (label i = 0; i < len; ++i)
        {
            const label count = listSizes[i];

            offsets_[i] = total;
            total += count;

            if (checkOverflow && total < offsets_[i])
            {
                reportOverflowAndExit(i, offsets_[i], count);
            }
        }

        offsets_[len] = total;
    }
    else
    {
        clear();
    }
    return total;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(const labelUList& listSizes)
{
    #ifdef FULLDEBUG
    const label total = resize_offsets(listSizes, true);
    #else
    const label total = resize_offsets(listSizes, false);
    #endif
    values_.resize(total);
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const labelUList& listSizes,
    const T& val
)
{
    #ifdef FULLDEBUG
    const label total = resize_offsets(listSizes, true);
    #else
    const label total = resize_offsets(listSizes, false);
    #endif
    values_.resize(total, val);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::CompactListList<T>::maxNonLocalSize(const label rowi) const
{
    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return 0;
    }

    label maxLen = 0;

    for (label i=0; i < len; ++i)
    {
        if (i != rowi)
        {
            const label count = (offsets_[i+1] - offsets_[i]);
            maxLen = max(maxLen, count);
        }
    }

    return maxLen;
}


template<class T>
std::streamsize Foam::CompactListList<T>::byteSize() const
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return this->size_bytes();
}


template<class T>
Foam::labelRange Foam::CompactListList<T>::range(const label i) const
{
    return labelRange(offsets_[i], offsets_[i+1] - offsets_[i]);
}


template<class T>
Foam::List<Foam::labelRange>
Foam::CompactListList<T>::ranges() const
{
    List<labelRange> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i=0; i < len; ++i)
    {
        values[i].reset(offsets_[i], (offsets_[i+1] - offsets_[i]));
    }

    return values;
}


template<class T>
void Foam::CompactListList<T>::resize(const labelUList& listSizes)
{
    #ifdef FULLDEBUG
    const label total = resize_offsets(listSizes, true);
    #else
    const label total = resize_offsets(listSizes, false);
    #endif
    values_.resize(total);
}


template<class T>
void Foam::CompactListList<T>::resize_nocopy(const labelUList& listSizes)
{
    #ifdef FULLDEBUG
    const label total = resize_offsets(listSizes, true);
    #else
    const label total = resize_offsets(listSizes, false);
    #endif
    values_.resize_nocopy(total);
}


template<class T>
void Foam::CompactListList<T>::setLocalSize(const label rowi, const label len)
{
    if (rowi >= 0 && rowi+1 < offsets_.size() && len >= 0)
    {
        const label delta = (len - (offsets_[rowi+1] - offsets_[rowi]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = rowi+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


template<class T>
Foam::labelList Foam::CompactListList<T>::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i=0; i < len; ++i)
    {
        values[i] = offsets_[i+1] - offsets_[i];
    }

    return values;
}


template<class T>
void Foam::CompactListList<T>::swap
(
    CompactListList<T>& other
)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    offsets_.swap(other.offsets_);
    values_.swap(other.values_);
}


template<class T>
void Foam::CompactListList<T>::transfer
(
    CompactListList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    offsets_.transfer(list.offsets_);
    values_.transfer(list.values_);
}


template<class T>
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const label pos,
    label len
) const
{
    if (pos >= 0 && pos < this->size())
    {
        // Change sub-length to (one-past) end position
        // len == -1 (like std::string::npos) - search until end

        if (len > 0) len += pos;
        if (len < 0 || len > this->size())
        {
            len = this->size();
        }

        for (label i = pos; i < len; ++i)
        {
            *d_iter = SubListType(this->localList(i));
            ++d_iter;
        }
    }

    return d_iter;
}


template<class T>
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const labelRange& range
) const
{
    return this->copy_unpack<SubListType>(d_iter, range.start(), range.size());
}


template<class T>
template<class SubListType, class OutputIter>
OutputIter Foam::CompactListList<T>::copy_unpack
(
    OutputIter d_iter,
    const labelUList& indices
) const
{
    for (label i : indices)
    {
        *d_iter = SubListType(this->localList(i));
        ++d_iter;
    }

    return d_iter;
}


// Could also support copy_unpack() with IndirectListBase, as required...
// or the caller can also just use copy_unpack with len = 1 and the
// desired position


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack() const
{
    List<SubListType> lists(size());

    this->copy_unpack<SubListType>(lists.begin());

    return lists;
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack(const labelRange& range) const
{
    List<SubListType> lists(range.size());

    this->copy_unpack<SubListType>(lists.begin(), range.start(), range.size());

    return lists;
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack(const labelUList& indices) const
{
    List<SubListType> lists(indices.size());

    this->copy_unpack<SubListType>(lists.begin(), indices);

    return lists;
}



// ************************************************************************* //
