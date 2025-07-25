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

#include "UList.H"
#include "contiguous.H"
#include "labelRange.H"

#include <random>

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
Foam::labelRange
Foam::UList<T>::validateRange(const labelRange& requestedRange) const
{
    const labelRange range(requestedRange.subset0(this->size()));

    #ifdef FULLDEBUG
    this->checkStart(range.start());
    this->checkSize(range.start() + range.size());
    #endif

    return range;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::moveFirst(const label i)
{
    checkIndex(i);

    for (label lower = 0; lower < i; ++lower)
    {
        Foam::Swap(this->operator[](lower), this->operator[](i));
    }
}


template<class T>
void Foam::UList<T>::moveLast(const label i)
{
    checkIndex(i);

    for (label upper = size()-1; upper > i; --upper)
    {
        Foam::Swap(this->operator[](i), this->operator[](upper));
    }
}


template<class T>
void Foam::UList<T>::swapFirst(const label i)
{
    checkIndex(i);

    if (i > 0)
    {
        Foam::Swap(this->operator[](0), this->operator[](i));
    }
}


template<class T>
void Foam::UList<T>::swapLast(const label i)
{
    checkIndex(i);

    const label upper = size()-1;

    if (i < upper)
    {
        Foam::Swap(this->operator[](i), this->operator[](upper));
    }
}


template<class T>
void Foam::UList<T>::deepCopy(const UList<T>& list)
{
    if (this->size_ != list.size_)
    {
        FatalErrorInFunction
            << "Lists have different sizes: "
            << this->size_ << " != " << list.size() << nl
            << abort(FatalError);
    }
    else if (this->size_ > 0)
    {
        // Can dispatch with
        // - std::execution::parallel_unsequenced_policy
        // - std::execution::unsequenced_policy
        std::copy(list.cbegin(), list.cend(), this->v_);
    }
}


template<class T>
template<class Addr>
void Foam::UList<T>::deepCopy(const IndirectListBase<T, Addr>& list)
{
    if (this->size_ != list.size())
    {
        FatalErrorInFunction
            << "Lists have different sizes: "
            << this->size_ << " != " << list.size() << nl
            << abort(FatalError);
    }
    else if (this->size_)
    {
        // Copy the indirect list contents

        // NB: operator[] for list read access (eg, an indirect list)
        // cannot replace with std::copy

        const label len = this->size_;

        auto iter = this->v_;

        for (label i = 0; i < len; (void)++i, (void)++iter)
        {
            *iter = list[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// This is non-inlined to allow template specializations
template<class T>
void Foam::UList<T>::operator=(const Foam::zero)
{
    this->fill_uniform(Foam::zero{});
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
std::streamsize Foam::UList<T>::byteSize() const
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
Foam::label Foam::UList<T>::find(const T& val) const
{
    const auto iter = std::find(this->cbegin(), this->cend(), val);
    return (iter != this->cend() ? label(iter - this->cbegin()) : label(-1));
}


template<class T>
Foam::label Foam::UList<T>::find(const T& val, label pos, label len) const
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

        const auto iter = std::find
        (
            (this->cbegin() + pos),
            (this->cbegin() + len),
            val
        );

        if (iter != (this->cbegin() + len))
        {
            return label(iter - this->cbegin());
        }
    }

    return -1;
}


template<class T>
Foam::label Foam::UList<T>::rfind(const T& val, label pos) const
{
    // pos == -1 (like std::string::npos) - search from end

    if (pos < 0 || pos >= this->size())
    {
        pos = this->size()-1;
    }

    while (pos >= 0)
    {
        if (this->v_[pos] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
bool Foam::UList<T>::operator==(const UList<T>& list) const
{
    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    return
    (
        (this->size_ == list.size_)
     && std::equal(this->cbegin(), this->cend(), list.cbegin())
    );
}


template<class T>
bool Foam::UList<T>::operator!=(const UList<T>& list) const
{
    return !operator==(list);
}


template<class T>
bool Foam::UList<T>::operator<(const UList<T>& list) const
{
    // Can dispatch with
    // - std::execution::parallel_unsequenced_policy
    // - std::execution::unsequenced_policy
    return std::lexicographical_compare
    (
        this->cbegin(), this->cend(),
        list.cbegin(), list.cend()
    );
}


template<class T>
bool Foam::UList<T>::operator>(const UList<T>& list) const
{
    return list.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator<=(const UList<T>& list) const
{
    return !list.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator>=(const UList<T>& list) const
{
    return !operator<(list);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::sort(UList<T>& list)
{
    // With which std::execution:: ?
    std::sort(list.begin(), list.end());
}


template<class T, class Compare>
void Foam::sort(UList<T>& list, const Compare& comp)
{
    // With which std::execution:: ?
    std::sort(list.begin(), list.end(), comp);
}


template<class T>
void Foam::stableSort(UList<T>& list)
{
    // With which std::execution:: ?
    std::stable_sort(list.begin(), list.end());
}


template<class T, class Compare>
void Foam::stableSort(UList<T>& list, const Compare& comp)
{
    // With which std::execution:: ?
    std::stable_sort(list.begin(), list.end(), comp);
}


template<class T>
void Foam::shuffle(UList<T>& list)
{
    std::shuffle(list.begin(), list.end(), std::default_random_engine());
}


// ************************************************************************* //
