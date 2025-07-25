/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

template<class IntType>
inline constexpr Foam::IntRange<IntType>::IntRange() noexcept
:
    start_(0),
    size_(0)
{}


template<class IntType>
inline constexpr Foam::IntRange<IntType>::IntRange
(
    const IntType len
) noexcept
:
    start_(0),
    size_(len)
{}


template<class IntType>
inline constexpr Foam::IntRange<IntType>::IntRange
(
    const IntType beg,
    const IntType len
) noexcept
:
    start_(beg),
    size_(len)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
inline IntType Foam::IntRange<IntType>::begin_value() const noexcept
{
    return (start_);
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::end_value() const noexcept
{
    return (start_ + size_);
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::rbegin_value() const noexcept
{
    return (start_ + (size_ - 1));
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::rend_value() const noexcept
{
    return (start_ - 1);
}


// * * * * * * * * * * * * * * Forward Iterators * * * * * * * * * * * * * * //

template<class IntType>
inline constexpr IntType
Foam::IntRange<IntType>::const_iterator::
operator[](const IntType n) const noexcept
{
    return (value_ + n);
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator&
Foam::IntRange<IntType>::const_iterator::
operator++() noexcept
{
    ++value_;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::const_iterator::
operator++(int) noexcept
{
    const_iterator old(*this);
    ++value_;
    return old;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator&
Foam::IntRange<IntType>::const_iterator::
operator--() noexcept
{
    --value_;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::const_iterator::
operator--(int) noexcept
{
    const_iterator old(*this);
    --value_;
    return old;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator&
Foam::IntRange<IntType>::const_iterator::
operator+=(const IntType n) noexcept
{
    value_ += n;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator&
Foam::IntRange<IntType>::const_iterator::
operator-=(const IntType n) noexcept
{
    value_ -= n;
    return *this;
}


template<class IntType>
inline constexpr typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::const_iterator::
operator+(const IntType n) const noexcept
{
    return const_iterator(value_ + n);
}


template<class IntType>
inline constexpr typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::const_iterator::
operator-(const IntType n) const noexcept
{
    return const_iterator(value_ - n);
}


template<class IntType>
inline constexpr IntType
Foam::IntRange<IntType>::const_iterator::
operator-(const const_iterator& iter) const noexcept
{
    return (value_ - iter.value_);
}


template<class IntType>
inline constexpr bool
Foam::IntRange<IntType>::const_iterator::
operator==(const const_iterator& iter) const noexcept
{
    return (value_ == iter.value_);
}


template<class IntType>
inline constexpr bool
Foam::IntRange<IntType>::const_iterator::
operator<(const const_iterator& iter) const noexcept
{
    return (value_ < iter.value_);
}


// * * * * * * * * * * * * * * Reverse Iterators * * * * * * * * * * * * * * //

template<class IntType>
inline constexpr IntType
Foam::IntRange<IntType>::const_reverse_iterator::
operator[](const IntType n) const noexcept
{
    return (value_ - n);
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator&
Foam::IntRange<IntType>::const_reverse_iterator::
operator++() noexcept
{
    --value_;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::const_reverse_iterator::
operator++(int) noexcept
{
    const_reverse_iterator old(*this);
    --value_;
    return old;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator&
Foam::IntRange<IntType>::const_reverse_iterator::
operator--() noexcept
{
    ++value_;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::const_reverse_iterator::
operator--(int) noexcept
{
    const_reverse_iterator old(*this);
    ++value_;
    return old;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator&
Foam::IntRange<IntType>::const_reverse_iterator::
operator+=(const IntType n) noexcept
{
    value_ -= n;
    return *this;
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator&
Foam::IntRange<IntType>::const_reverse_iterator::
operator-=(const IntType n) noexcept
{
    value_ += n;
    return *this;
}


template<class IntType>
inline constexpr typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::const_reverse_iterator::
operator+(const IntType n) const noexcept
{
    return const_reverse_iterator(value_ - n);
}


template<class IntType>
inline constexpr typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::const_reverse_iterator::
operator-(const IntType n) const noexcept
{
    return const_reverse_iterator(value_ + n);
}


template<class IntType>
inline constexpr IntType
Foam::IntRange<IntType>::const_reverse_iterator::
operator-(const const_reverse_iterator& iter) const noexcept
{
    return (iter.value_ - value_);
}


template<class IntType>
inline constexpr bool
Foam::IntRange<IntType>::const_reverse_iterator::
operator==(const const_reverse_iterator& iter) const noexcept
{
    return (value_ == iter.value_);
}


template<class IntType>
inline constexpr bool
Foam::IntRange<IntType>::const_reverse_iterator::
operator<(const const_reverse_iterator& iter) const noexcept
{
    return (value_ > iter.value_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::at(const IntType i) const
{
    return const_iterator(start_ + ((i < 0 || i > size_) ? size_ : i));
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::begin() const noexcept
{
    return const_iterator(begin_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::cbegin() const noexcept
{
    return const_iterator(begin_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::end() const noexcept
{
    return const_iterator(end_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_iterator
Foam::IntRange<IntType>::cend() const noexcept
{
    return const_iterator(end_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::rbegin() const noexcept
{
    return const_reverse_iterator(rbegin_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::crbegin() const noexcept
{
    return const_reverse_iterator(rbegin_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::rend() const noexcept
{
    return const_reverse_iterator(rend_value());
}


template<class IntType>
inline typename Foam::IntRange<IntType>::const_reverse_iterator
Foam::IntRange<IntType>::crend() const noexcept
{
    return const_reverse_iterator(rend_value());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
inline void Foam::IntRange<IntType>::reset() noexcept
{
    start_ = size_ = 0;
}


template<class IntType>
inline void Foam::IntRange<IntType>::reset
(
    const IntType beg,
    const IntType len
) noexcept
{
    start_ = beg;
    size_  = len;
}


template<class IntType>
inline void Foam::IntRange<IntType>::setStart(const IntType i) noexcept
{
    start_ = i;
}


template<class IntType>
inline void Foam::IntRange<IntType>::setSize(const IntType n) noexcept
{
    size_ = n;
}


template<class IntType>
inline void Foam::IntRange<IntType>::resize(const IntType n) noexcept
{
    size_ = n;
}


template<class IntType>
inline void Foam::IntRange<IntType>::clampSize() noexcept
{
    if (size_ < 0) size_ = 0;
}


template<class IntType>
inline bool Foam::IntRange<IntType>::contains(const IntType value)
const noexcept
{
    // Avoid overflow risk:
    //   (value < (start_ + size_)) --> ((value - start_) < size_)

    return (size_ && start_ <= value && (value - start_) < size_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class IntType>
inline constexpr IntType Foam::IntRange<IntType>::
operator[](const IntType i) const noexcept
{
    return start_ + i;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator++() noexcept
{
    return ++size_;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator++(int) noexcept
{
    return size_++;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator--() noexcept
{
    --size_;
    clampSize();
    return size_;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator--(int) noexcept
{
    const IntType old(size_);
    --size_;
    clampSize();
    return old;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator+=(const IntType n) noexcept
{
    size_ += n;
    return size_;
}


template<class IntType>
inline IntType Foam::IntRange<IntType>::
operator-=(const IntType n) noexcept
{
    size_ -= n;
    clampSize();
    return size_;
}


// ************************************************************************* //
