/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicList<T, SizeMin>::doAssignDynList
(
    const ListType& list
)
{
    const label len = list.size();

    if (capacity_ < len)
    {
        // Needs more space for the copy operation
        List<T>::setAddressableSize(capacity_);  // Use entire space
        List<T>::resize_nocopy(len);
        capacity_ = List<T>::size();
    }

    // Perform copy into addressable portion
    List<T>::setAddressableSize(len);
    List<T>::operator=(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::doCapacity
(
    const bool nocopy,
    const label newCapacity
)
{
    if (newCapacity == capacity_)
    {
        return;
    }

    // Addressable length, possibly truncated by new capacity
    const label currLen = min(List<T>::size(), newCapacity);

    // Corner case...
    if (List<T>::size() == newCapacity)
    {
        // Adjust addressable size to trigger proper resizing.
        // Using (old size+1) is safe since it does not affect the 'overlap'
        // of old and new addressable regions, but incurs fewew copy
        // operations than extending to use the current capacity would.
        List<T>::setAddressableSize(currLen+1);
    }

    if (nocopy)
    {
        List<T>::resize_nocopy(newCapacity);
    }
    else
    {
        List<T>::resize(newCapacity);
    }

    capacity_ = List<T>::size();
    List<T>::setAddressableSize(currLen);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::doReserve
(
    const bool nocopy,
    const label len
)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = List<T>::size();

        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(len, label(2*capacity_)));

        if (nocopy)
        {
            List<T>::resize_nocopy(capacity_);
        }
        else
        {
            List<T>::resize(capacity_);
        }
        List<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::doResize
(
    const bool nocopy,
    const label len
)
{
    this->doReserve(nocopy, len);
    List<T>::setAddressableSize(len);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::DynamicList<T, SizeMin>::DynamicList() noexcept
:
    List<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList(const label initialCapacity)
:
    List<T>(),
    capacity_(0)
{
    reserve_nocopy(initialCapacity);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    const T& val
)
:
    List<T>(len, val),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    const Foam::zero
)
:
    List<T>(len, Foam::zero{}),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, SizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, AnySizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& list,
    const labelUList& indices
)
:
    List<T>(list, indices),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<unsigned N>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const FixedList<T, N>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    std::initializer_list<T> lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const IndirectListBase<T, Addr>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, SizeMin>&& list
) noexcept
:
    List<T>(std::move(static_cast<List<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, AnySizeMin>&& list
) noexcept
:
    List<T>(std::move(static_cast<List<T>&>(list))),
    capacity_(list.capacity())
{
    list.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    List<T>&& list
) noexcept
:
    List<T>(std::move(list)),
    capacity_(List<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline std::streamsize
Foam::DynamicList<T, SizeMin>::capacity_bytes() const noexcept
{
    return std::streamsize(capacity_)*sizeof(T);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity
(
    const label len
)
{
    this->doCapacity(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity_nocopy
(
    const label len
)
{
    this->doCapacity(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity_unsafe
(
    const label len
) noexcept
{
    capacity_ = len;
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve
(
    const label len
)
{
    this->doReserve(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve_nocopy
(
    const label len
)
{
    this->doReserve(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label len
)
{
    this->doResize(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize_fill
(
    const label len,
    const T& val
)
{
    this->doResize(true, len);  // nocopy = true
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize_nocopy
(
    const label len
)
{
    this->doResize(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label len,
    const T& val
)
{
    const label oldLen = List<T>::size();
    resize(len);

    // Fill newly exposed with constant value
    if (oldLen < List<T>::size())
    {
        std::fill
        (
            this->begin(oldLen), this->end(), val
        );
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clear() noexcept
{
    List<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::shrink_to_fit()
{
    const label currLen = List<T>::size();
    if (currLen < capacity_)
    {
        // Adjust addressable size to trigger proper resizing
        List<T>::setAddressableSize(currLen+1);
        List<T>::resize(currLen);
        capacity_ = List<T>::size();
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::shrink_unsafe()
{
    if (List<T>::empty())
    {
        // Delete storage if empty
        List<T>::clear();
    }
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::shrink()
{
    this->shrink_to_fit();
    return *this;
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::swap(List<T>& list)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Remove unused storage
    this->shrink_to_fit();

    // Swap storage and addressable size
    UList<T>::swap(list);

    // Update capacity
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::swap
(
    DynamicList<T, AnySizeMin>& other
) noexcept
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&other)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer(List<T>& list)
{
    List<T>::transfer(list);
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    DynamicList<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    // Take over storage as-is (without shrink)
    capacity_ = list.capacity();

    List<T>::transfer(static_cast<List<T>&>(list));
    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::DynamicList<T, SizeMin>::emplace_back(Args&&... args)
{
    // This could/should be better with inplace construction
    // (as per std::vector), but currently lacking the methods for that
    // so resize and move assign

    const label idx = List<T>::size();
    resize(idx + 1);

    // move assign element
    UList<T>::operator[](idx) = T(std::forward<Args>(args)...);
    return UList<T>::operator[](idx);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const T& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    UList<T>::operator[](idx) = val;  // copy element
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    T&& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    UList<T>::operator[](idx) = std::move(val);  // move assign element
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const UList<T>& list
)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
template<unsigned N>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const FixedList<T, N>& list
)
{
    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    std::initializer_list<T> list
)
{
    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    const IndirectListBase<T, Addr>& list
)
{
    // Note: push_back will still work even if the indirect list
    // actually references *this, since its source elements will not
    // overlap the new destinations.

    const label idx = this->size();
    const label n = list.size();

    resize(list + n);

    auto iter = this->begin(idx);

    for (label i = 0; i < n; (void)++i, (void)++iter)
    {
        *iter = list[i];  // copy element
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    List<T>&& list
)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    std::move(list.begin(), list.end(), this->begin(idx));

    list.clear();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::push_back
(
    DynamicList<T, AnySizeMin>&& list
)
{
    push_back(std::move(static_cast<List<T>&>(list)));
    list.clearStorage();  // Deletion, capacity=0 etc.
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::push_uniq(const T& val)
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


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::pop_back(label n)
{
    if (n >= this->size())
    {
        clear();
    }
    else if (n > 0)
    {
        resize(this->size() - n);
    }
}


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = List<T>::size() - 1;

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = List<T>::operator[](idx);

    List<T>::setAddressableSize(idx);

    return val;
}


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove
(
    const label idx,
    const bool fast
)
{
    if (fast)
    {
        // Simply swap idx <=> last
        this->swapLast(idx);
    }
    else
    {
        // Move element to the end and move everything down
        this->moveLast(idx);
    }

    // Element to remove is now at the end
    return this->remove();
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::remove
(
    const labelRange& range
)
{
    return this->removeElements(this->validateRange(range));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::remove
(
    std::initializer_list<label> start_size
)
{
    return this->removeElements(this->validateRange(start_size));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::subset
(
    const labelRange& range
)
{
    return this->subsetElements(this->validateRange(range));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::subset
(
    std::initializer_list<label> start_size
)
{
    return this->subsetElements(this->validateRange(start_size));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynamicList<T, SizeMin>::operator()
(
    const label i
)
{
    if (i >= List<T>::size())
    {
        resize(i + 1);
    }

    return UList<T>::operator[](i);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const Foam::zero
)
{
    UList<T>::operator=(Foam::zero{});
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const UList<T>& lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
template<unsigned N>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const FixedList<T, N>& lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, SizeMin>& lst
)
{
    if (this == &lst)
    {
        return;  // Self-assignment is a no-op
    }

    doAssignDynList(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    doAssignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    std::initializer_list<T> lst
)
{
    doAssignDynList(lst);
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const IndirectListBase<T, Addr>& lst
)
{
    // NOTE: Self-assignment needs special handling
    /// if
    /// (
    ///     static_cast<const UList<T>*>(this)
    ///  == static_cast<const UList<T>*>(&list.values())
    /// )

    doAssignDynList(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    List<T>&& lst
)
{
    clear();
    transfer(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, SizeMin>&& lst
)
{
    if (this == &lst)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    transfer(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, AnySizeMin>&& list
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    transfer(list);
}


// ************************************************************************* //
