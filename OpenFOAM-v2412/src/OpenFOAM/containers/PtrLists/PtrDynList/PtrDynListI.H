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

#include "autoPtr.H"
#include "refPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::PtrDynList<T, SizeMin>::PtrDynList() noexcept
:
    PtrList<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(const label len)
:
    PtrList<T>(),
    capacity_(0)
{
    reserve(len);
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    const PtrDynList<T, SizeMin>& list
)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, SizeMin>&& list
)
:
    PtrList<T>(std::move(list)),
    capacity_(list.capacity())
{
    // FUTURE:
    // list.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept

    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, AnySizeMin>&& list
)
:
    PtrList<T>(std::move(list)),
    capacity_(list.capacity())
{
    // FUTURE:
    // list.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept

    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrList<T>&& list
) noexcept
:
    PtrList<T>(std::move(list)),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(UList<T*>& list)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reserve(const label len)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = PtrList<T>::size();

        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(len, label(2*capacity_)));

        PtrList<T>::resize(capacity_);
        PtrList<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize(const label newLen)
{
    auto& ptrs = this->ptrs_;

    const label oldLen = ptrs.size();

    if (capacity_ < newLen)
    {
        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(newLen, label(2*capacity_)));

        PtrList<T>::resize(capacity_);
    }
    else if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i = newLen; i < oldLen; ++i)
        {
            delete ptrs[i];
            ptrs[i] = nullptr;
        }
    }

    // Adjust addressed size
    PtrList<T>::setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize_null(const label newLen)
{
    if (capacity_ < newLen)
    {
        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(newLen, label(2*capacity_)));

        PtrList<T>::resize_null(capacity_);
    }
    else
    {
        PtrList<T>::free();  // Free (and nullify) old pointers
    }

    // Adjust addressed size
    PtrList<T>::setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clear()
{
    PtrList<T>::free();  // Free (and nullify) old pointers
    PtrList<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clearStorage()
{
    PtrList<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::shrink_to_fit()
{
    const label currLen = PtrList<T>::size();
    if (currLen < capacity_)
    {
        // Adjust addressable size to trigger proper resizing
        PtrList<T>::setAddressableSize(currLen+1);
        PtrList<T>::resize(currLen);
        capacity_ = PtrList<T>::size();
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::shrink_unsafe()
{
    if (PtrList<T>::empty())
    {
        // Delete empty list
        PtrList<T>::clear();
    }
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::squeezeNull()
{
    const label newLen = UPtrList<T>::squeezeNull();
    PtrList<T>::setAddressableSize(newLen);
    return newLen;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::swap(PtrList<T>& list)
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&list)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Remove unused storage
    this->shrink_to_fit();

    // Swap storage and addressable size
    UPtrList<T>::swap(list);

    // Update capacity
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::swap
(
    PtrDynList<T, AnySizeMin>& other
) noexcept
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&other)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UPtrList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::transfer(PtrList<T>& list)
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&list)
    )
    {
        return;  // Self assignment is a no-op
    }

    PtrList<T>::transfer(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::transfer
(
    PtrDynList<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&list)
    )
    {
        return;  // Self assignment is a no-op
    }

    // Take over storage as-is (without shrink)
    capacity_ = list.capacity();

    PtrList<T>::transfer(static_cast<PtrList<T>&>(list));
    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace_back(Args&&... args)
{
    T* ptr = new T(std::forward<Args>(args)...);
    this->push_back(ptr);
    return *ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(T* ptr)
{
    const label idx = this->size();
    resize(idx + 1);
    this->ptrs_[idx] = ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(std::unique_ptr<T>&& ptr)
{
    this->push_back(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(autoPtr<T>&& ptr)
{
    this->push_back(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(const refPtr<T>& ptr)
{
    this->push_back(ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(const tmp<T>& ptr)
{
    this->push_back(ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back(PtrList<T>&& other)
{
    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i = 0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // Take pointer ownership
    }

    other.clear();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::push_back
(
    PtrDynList<T, AnySizeMin>&& other
)
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&other)
    )
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i = 0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // Take pointer ownership
    }

    other.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::pop_back(label n)
{
    if (n >= this->size())
    {
        this->clear();
    }
    else if (n > 0)
    {
        this->resize(this->size() - n);
    }
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace_set
(
    const label i,
    Args&&... args
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }
    return PtrList<T>::emplace_set(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::emplace
(
    const label i,
    Args&&... args
)
{
    return this->emplace_set(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::PtrDynList<T, SizeMin>::try_emplace
(
    const label i,
    Args&&... args
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }
    return PtrList<T>::try_emplace(i, std::forward<Args>(args)...);
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    T* ptr
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }

    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    std::unique_ptr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    autoPtr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const refPtr<T>& ptr
)
{
    return this->set(i, ptr.ptr());  // release or clone

}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const tmp<T>& ptr
)
{
    return this->set(i, ptr.ptr());  // release or clone
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reorder(const labelUList& oldToNew)
{
    // Shrinking first is a bit annoying, but saves needing a special version.
    this->shrink_to_fit();
    PtrList<T>::reorder(oldToNew);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, SizeMin>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const PtrList<T>*>(this)
     == static_cast<const PtrList<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrList<T>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, SizeMin>&& list
)
{
    this->transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, AnySizeMin>&& list
)
{
    this->transfer(list);
}


// ************************************************************************* //
