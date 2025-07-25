/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::CircularBuffer<T>::toGlobal(label i) const
{
    const label len = this->size();

    if (!len)
    {
        // Bounds error
        return -1;
    }
    else if (i < 0)
    {
        // Wrap any number of times
        while (i < 0) i += len;
    }
    else
    {
        // Wrap any number of times
        while (i >= len) i -= len;
    }

    i += begin_;

    if (i >= storage_.size())
    {
        i -= storage_.size();
    }

    return i;
}


template<class T>
inline Foam::label Foam::CircularBuffer<T>::size_one() const noexcept
{
    return
    (
        (end_ >= begin_)
      ? (end_ - begin_)
      : (storage_.size() - begin_)
    );
}


template<class T>
inline Foam::label Foam::CircularBuffer<T>::size_two() const noexcept
{
    return
    (
        (end_ && end_ < begin_)
      ? end_
      : static_cast<label>(0)
    );
}


template<class T>
template<class OtherListType>
inline void Foam::CircularBuffer<T>::copyList(const OtherListType& rhs)
{
    this->clear();

    const label len = rhs.size();

    if (len)
    {
        reserve(len + 1);

        // Never overfilled, simply write at end_ (one-past position)

        // - after clear(), begin_ and end_ are both 0

        for (label i = 0; i < len; ++i)
        {
            storage_[end_] = rhs[i];  // copy element
            ++end_;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::CircularBuffer<T>::CircularBuffer() noexcept
:
    storage_(),
    begin_(0),
    end_(0)
{}


template<class T>
inline Foam::CircularBuffer<T>::CircularBuffer(const label len)
:
    storage_(max(min_size(), len + 1)),
    begin_(0),
    end_(0)
{}


template<class T>
inline Foam::CircularBuffer<T>::CircularBuffer
(
    const CircularBuffer<T>& list
)
:
    storage_(list.storage_),
    begin_(list.begin_),
    end_(list.end_)
{}


template<class T>
inline Foam::CircularBuffer<T>::CircularBuffer
(
    CircularBuffer<T>&& list
)
:
    storage_(std::move(list.storage_)),
    begin_(list.begin_),
    end_(list.end_)
{
    list.begin_ = 0;
    list.end_ = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::CircularBuffer<T>::capacity() const noexcept
{
    // or storage_.capacity();
    return storage_.size();
}


template<class T>
inline bool Foam::CircularBuffer<T>::empty() const noexcept
{
    return storage_.empty() || (begin_ == end_);
}


template<class T>
inline Foam::label Foam::CircularBuffer<T>::size() const noexcept
{
    const label diff(end_ - begin_);

    if (diff < 0)
    {
        return (storage_.size() + diff);
    }

    return diff;
}


template<class T>
inline Foam::label Foam::CircularBuffer<T>::space() const noexcept
{
    return (storage_.size() - size());
}


template<class T>
inline Foam::labelRange Foam::CircularBuffer<T>::range_one() const noexcept
{
    return
    (
        (begin_ == end_)
      ? labelRange()
      : labelRange(begin_, this->size_one())
    );
}


template<class T>
inline Foam::labelRange Foam::CircularBuffer<T>::range_two() const noexcept
{
    return labelRange(0, this->size_two());
}


template<class T>
inline void Foam::CircularBuffer<T>::clear() noexcept
{
    begin_ = end_ = 0;
}


template<class T>
inline void Foam::CircularBuffer<T>::clearStorage()
{
    storage_.clear();
    begin_ = end_ = 0;
}


template<class T>
inline void Foam::CircularBuffer<T>::swap(CircularBuffer<T>& other)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressing
    storage_.swap(other.storage_);
    std::swap(begin_, other.begin_);
    std::swap(end_, other.end_);
}


template<class T>
inline void Foam::CircularBuffer<T>::reserve(const label len)
{
    this->doReserve(false, len);
}


template<class T>
inline void Foam::CircularBuffer<T>::reserve_nocopy(const label len)
{
    this->doReserve(true, len);
}


template<class T>
bool Foam::CircularBuffer<T>::contains(const T& val) const
{
    return (this->array_one().contains(val) || this->array_two().contains(val));
}


template<class T>
inline bool Foam::CircularBuffer<T>::contains(const T& val, label pos) const
{
    return (this->find(val, pos) >= 0);
}


template<class T>
inline T& Foam::CircularBuffer<T>::front()
{
    if (empty())
    {
        FatalErrorInFunction << "Buffer is empty" << abort(FatalError);
    }

    return storage_[begin_];
}


template<class T>
inline const T& Foam::CircularBuffer<T>::front() const
{
    if (empty())
    {
        FatalErrorInFunction << "Buffer is empty" << abort(FatalError);
    }

    return storage_[begin_];
}


template<class T>
inline T& Foam::CircularBuffer<T>::back()
{
    if (empty())
    {
        FatalErrorInFunction << "Buffer is empty" << abort(FatalError);
    }

    return storage_.rcValue(end_);
}


template<class T>
inline const T& Foam::CircularBuffer<T>::back() const
{
    if (empty())
    {
        FatalErrorInFunction << "Buffer is empty" << abort(FatalError);
    }

    return storage_.rcValue(end_);
}


template<class T>
inline void Foam::CircularBuffer<T>::push_front(const T& val)
{
    reserve(size() + 2);

    // Never overfilled. Move begin and write

    begin_ = storage_.rcIndex(begin_);
    storage_[begin_] = val;  // copy assign element
}


template<class T>
inline void Foam::CircularBuffer<T>::push_front(T&& val)
{
    reserve(size() + 2);

    // Never overfilled. Move begin and write

    begin_ = storage_.rcIndex(begin_);
    storage_[begin_] = std::move(val);  // move assign element
}


template<class T>
template<class... Args>
inline T& Foam::CircularBuffer<T>::emplace_front(Args&&... args)
{
    reserve(size() + 2);

    // Never overfilled. Move begin and write

    begin_ = storage_.rcIndex(begin_);
    storage_[begin_] = T(std::forward<Args>(args)...);

    return storage_[begin_];
}


template<class T>
inline void Foam::CircularBuffer<T>::push_back(const T& val)
{
    reserve(size() + 2);

    // Never overfilled, simply write at end_ (one-past position)

    storage_[end_] = val;  // copy assign element
    end_ = storage_.fcIndex(end_);
}


template<class T>
inline void Foam::CircularBuffer<T>::push_back(T&& val)
{
    reserve(size() + 2);

    // Never overfilled, simply write at end_ (one-past position)

    storage_[end_] = std::move(val);  // move assign element
    end_ = storage_.fcIndex(end_);
}


template<class T>
template<class... Args>
inline T& Foam::CircularBuffer<T>::emplace_back(Args&&... args)
{
    reserve(size() + 2);

    // Never overfilled, simply write at end_ (one-past position)

    const label backIndex = end_;
    storage_[end_] = T(std::forward<Args>(args)...);
    end_ = storage_.fcIndex(end_);

    return storage_[backIndex];
}


template<class T>
inline void Foam::CircularBuffer<T>::pop_front(label n)
{
    if (n >= size())
    {
        begin_ = end_;
    }
    else
    {
        while (n-- > 0)
        {
            begin_ = storage_.fcIndex(begin_);
        }
    }
}


template<class T>
inline void Foam::CircularBuffer<T>::pop_back(label n)
{
    if (n >= size())
    {
        end_ = begin_;
    }
    else
    {
        while (n-- > 0)
        {
            end_ = storage_.rcIndex(end_);
        }
    }
}


template<class T>
inline Foam::label Foam::CircularBuffer<T>::push_uniq(const T& val)
{
    if (this->contains(val))
    {
        return 0;
    }
    else
    {
        this->push_back(val);
        return 1;  // Increased list length by one
    }
}


template<class T>
inline void Foam::CircularBuffer<T>::push_back(const UList<T>& rhs)
{
    const label len = rhs.size();

    if (len)
    {
        reserve(size() + len + 1);

        // Never overfilled, simply write at end_ (one-past position)

        for (label i = 0; i < len; ++i)
        {
            storage_[end_] = rhs[i];  // copy element
            end_ = storage_.fcIndex(end_);
        }
    }
}


template<class T>
template<class Addr>
inline void Foam::CircularBuffer<T>::push_back
(
    const IndirectListBase<T, Addr>& rhs
)
{
    const label len = rhs.size();

    if (len)
    {
        reserve(size() + len + 1);

        // Never overfilled, simply write at end_ (one-past position)

        for (label i = 0; i < len; ++i)
        {
            storage_[end_] = rhs[i];  // copy element
            end_ = storage_.fcIndex(end_);
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::CircularBuffer<T>::operator[](label i)
{
    const label idx = this->toGlobal(i);
    return storage_[idx];
}


template<class T>
inline const T& Foam::CircularBuffer<T>::operator[](label i) const
{
    const label idx = this->toGlobal(i);
    return storage_[idx];
}


template<class T>
inline void Foam::CircularBuffer<T>::operator=(const CircularBuffer<T>& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->clear();

    const auto list1 = rhs.array_one();
    const auto list2 = rhs.array_two();
    const label len = (list1.size() + list2.size());

    if (len)
    {
        reserve(len + 1);

        // Never overfilled, simply write at end_ (one-past position)

        // - after clear(), begin_ and end_ are both 0

        for (const T& val : list1)
        {
            storage_[end_] = val;
            ++end_;
        }

        for (const T& val : list2)
        {
            storage_[end_] = val;
            ++end_;
        }
    }
}


template<class T>
inline void Foam::CircularBuffer<T>::operator=(CircularBuffer<T>&& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->clearStorage();
    this->swap(rhs);
}


template<class T>
inline void Foam::CircularBuffer<T>::operator=(const T& val)
{
    this->array_one() = val;
    this->array_two() = val;
}


template<class T>
inline void Foam::CircularBuffer<T>::operator=(const Foam::zero)
{
    this->array_one() = Foam::zero{};
    this->array_two() = Foam::zero{};
}


template<class T>
inline void Foam::CircularBuffer<T>::operator=(const UList<T>& rhs)
{
    this->copyList(rhs);
}


template<class T>
template<class AnyAddr>
inline void Foam::CircularBuffer<T>::operator=
(
    const IndirectListBase<T, AnyAddr>& rhs
)
{
    this->copyList(rhs);
}


// ************************************************************************* //
