/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::PtrList<T>::PtrList() noexcept
:
    UPtrList<T>()
{}


template<class T>
inline Foam::PtrList<T>::PtrList(const label len)
:
    UPtrList<T>(len)
{}


template<class T>
inline Foam::PtrList<T>::PtrList(const PtrList<T>& list)
:
    UPtrList<T>(list.ptrs_.clone())
{}


template<class T>
inline Foam::PtrList<T>::PtrList(PtrList<T>&& list) noexcept
:
    UPtrList<T>(std::move(list))
{}


template<class T>
inline Foam::PtrList<T>::PtrList(UList<T*>& list)
:
    UPtrList<T>(list)
{
    // Took ownership of the pointers
    list = static_cast<T*>(nullptr);
}


template<class T>
template<class CloneArg>
inline Foam::PtrList<T>::PtrList
(
    const PtrList<T>& list,
    const CloneArg& cloneArg
)
:
    UPtrList<T>(list.clone(cloneArg)())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::clear()
{
    (this->ptrs_).free();  // Free (and nullify) old pointers
    UPtrList<T>::clear();
}


template<class T>
inline void Foam::PtrList<T>::free()
{
    (this->ptrs_).free();  // Free and nullify old pointers
}


template<class T>
inline void Foam::PtrList<T>::resize_null(const label newLen)
{
    (this->ptrs_).free();  // Free (and nullify) old pointers
    UPtrList<T>::resize_null(newLen);
}


template<class T>
template<class... Args>
inline T& Foam::PtrList<T>::emplace_back(Args&&... args)
{
    T* ptr = new T(std::forward<Args>(args)...);
    UPtrList<T>::push_back(ptr);
    return *ptr;
}


template<class T>
inline void Foam::PtrList<T>::push_back(T* ptr)
{
    UPtrList<T>::push_back(ptr);
}


template<class T>
inline void Foam::PtrList<T>::push_back(std::unique_ptr<T>&& ptr)
{
    UPtrList<T>::push_back(ptr.release());
}


template<class T>
inline void Foam::PtrList<T>::push_back(autoPtr<T>&& ptr)
{
    UPtrList<T>::push_back(ptr.release());
}


template<class T>
inline void Foam::PtrList<T>::push_back(const refPtr<T>& ptr)
{
    UPtrList<T>::push_back(ptr.ptr());  // release or clone
}


template<class T>
inline void Foam::PtrList<T>::push_back(const tmp<T>& ptr)
{
    UPtrList<T>::push_back(ptr.ptr());  // release or clone
}


template<class T>
inline void Foam::PtrList<T>::push_back(PtrList<T>&& other)
{
    if (this == &other)
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

    other.clear();
}


template<class T>
template<class... Args>
inline T& Foam::PtrList<T>::emplace_set(const label i, Args&&... args)
{
    (void) this->release(i);  // delete old entry
    T* ptr = new T(std::forward<Args>(args)...);
    (void) UPtrList<T>::set(i, ptr);
    return *ptr;
}


template<class T>
template<class... Args>
inline T& Foam::PtrList<T>::emplace(const label i, Args&&... args)
{
    return this->emplace_set(i, std::forward<Args>(args)...);
}


template<class T>
template<class... Args>
inline T& Foam::PtrList<T>::try_emplace(const label i, Args&&... args)
{
    T* ptr = UPtrList<T>::get(i);
    if (ptr)
    {
        return *ptr;
    }
    return this->emplace_set(i, std::forward<Args>(args)...);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
    // UPtrList::set returns a nullptr if trying to set with the same
    // pointer (again). This prevents the autoPtr from managing the
    // memory (avoids possible double deletion).

    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    std::unique_ptr<T>&& ptr
)
{
    return set(i, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    autoPtr<T>&& ptr
)
{
    return set(i, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const refPtr<T>& ptr
)
{
    return set(i, ptr.ptr());  // release or clone
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const tmp<T>& ptr
)
{
    return set(i, ptr.ptr());  // release or clone
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::release(const label i)
{
    if (i < 0 || i >= this->size())
    {
        return nullptr;
    }

    return autoPtr<T>(UPtrList<T>::set(i, nullptr));
}


template<class T>
inline void Foam::PtrList<T>::transfer(PtrList<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    (this->ptrs_).free();  // Free and nullify old pointers
    UPtrList<T>::transfer(list);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::operator=(PtrList<T>&& list)
{
    this->transfer(list);
}


// ************************************************************************* //
