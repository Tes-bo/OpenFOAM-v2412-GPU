/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

inline Foam::labelRanges::labelRanges(const label initialCapacity)
:
    ranges_(initialCapacity)
{}


inline Foam::labelRanges::labelRanges(const UList<labelRange>& list)
:
    ranges_(list)
{}


inline Foam::labelRanges::labelRanges(List<labelRange>&& list)
:
    ranges_(std::move(list))
{}


template<int AnySizeMin>
inline Foam::labelRanges::labelRanges
(
    DynamicList<labelRange, AnySizeMin>&& list
)
:
    ranges_(std::move(list))
{}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline constexpr Foam::labelRanges::const_iterator::
const_iterator
(
    const UList<labelRange>* list,
    const label idx,
    const label subIdx
) noexcept
:
    list_(list),
    index_(idx),
    subIndex_(subIdx)
{}


inline Foam::label Foam::labelRanges::const_iterator::
operator*() const
{
    return (*list_)[index_][subIndex_];
}


inline Foam::labelRanges::const_iterator&
Foam::labelRanges::const_iterator::
operator++()
{
    // TBD: trap (index_ >= list_->size()) and make a no-op?
    if (++subIndex_ >= (*list_)[index_].size())
    {
        // Move to the next range
        ++index_;
        subIndex_ = 0;
    }

    return *this;
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::const_iterator::
operator++(int)
{
    const_iterator old(*this);
    this->operator++();
    return old;
}


inline constexpr bool
Foam::labelRanges::const_iterator::
operator==
(
    const const_iterator& iter
) const noexcept
{
    return
    (
        index_    == iter.index_
     && subIndex_ == iter.subIndex_
    );
}


inline constexpr bool
Foam::labelRanges::const_iterator::
operator!=
(
    const const_iterator& iter
) const noexcept
{
    return !(*this == iter);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::labelRanges::totalSize() const noexcept
{
    label total = 0;
    for (const labelRange& range : ranges_)
    {
        if (range.size() > 0)  // Ignore negative size (paranoid)
        {
            total += range.size();
        }
    }
    return total;
}


inline bool Foam::labelRanges::contains(const label value) const noexcept
{
    for (const labelRange& range : ranges_)
    {
        if (range.contains(value))
        {
            return true;
        }
    }

    return false;
}


template<class... Args>
inline Foam::labelRange& Foam::labelRanges::emplace_back(Args&&... args)
{
    return ranges_.emplace_back(args...);
}


inline void Foam::labelRanges::sort()
{
    Foam::sort(ranges_);
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::labelRanges::const_iterator
Foam::labelRanges::cbegin() const noexcept
{
    return const_iterator(&ranges_);
}


inline const Foam::labelRanges::const_iterator
Foam::labelRanges::cend() const noexcept
{
    return const_iterator(&ranges_, ranges_.size());
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::begin() const noexcept
{
    return const_iterator(&ranges_);
}


inline const Foam::labelRanges::const_iterator
Foam::labelRanges::end() const noexcept
{
    return const_iterator(&ranges_, ranges_.size());
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::cbegin(const label i) const
{
    if (i <= 0) return this->cbegin();

    label idx = 0;
    label subIdx = i;

    for (const labelRange& range : ranges_)
    {
        if (subIdx < range.size())
        {
            return const_iterator(&ranges_, idx, subIdx);
        }
        else
        {
            ++idx;
            if (range.size() > 0)  // Ignore negative size (paranoid)
            {
                subIdx -= range.size();
            }
        }
    }

    return this->cend();
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::begin(const label i) const
{
    return this->cbegin(i);
}


// ************************************************************************* //
