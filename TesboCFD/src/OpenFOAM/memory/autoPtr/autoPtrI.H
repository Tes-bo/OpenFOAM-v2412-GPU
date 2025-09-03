/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T* Foam::autoPtr<T>::release() noexcept
{
    T* p = ptr_;
    ptr_ = nullptr;
    return p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(T* p) noexcept
{
    delete ptr_;
    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(autoPtr<T>&& other) noexcept
{
    // Could also make Fatal with FULLDEBUG
    if (&other == this)
    {
        return;  // No self-assignment
    }

    reset(other.release());
}


template<class T>
template<class... Args>
inline T& Foam::autoPtr<T>::emplace(Args&&... args)
{
    delete ptr_;  // delete old entry
    ptr_ = new T(std::forward<Args>(args)...);
    return *ptr_;
}


template<class T>
inline void Foam::autoPtr<T>::swap(autoPtr<T>& other) noexcept
{
    // Swap is just copy/assign for pointer and enum types
    // Self-swap is effectively ignored
    T* p = ptr_;
    ptr_ = other.ptr_;
    other.ptr_ = p;
}


template<class T>
template<class... Args>
inline Foam::autoPtr<T> Foam::autoPtr<T>::clone(Args&&... args) const
{
    if (ptr_)
    {
        return autoPtr<T>(ptr_->clone(std::forward<Args>(args)...).ptr());
    }

    return autoPtr<T>();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::autoPtr<T>::operator*()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "unallocated autoPtr of type " << typeid(T).name()
            << abort(FatalError);
    }
    return *ptr_;
}


template<class T>
inline const T& Foam::autoPtr<T>::operator*() const
{
    return const_cast<autoPtr<T>*>(this)->operator*();
}


template<class T>
inline T* Foam::autoPtr<T>::operator->()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "unallocated autoPtr of type " << typeid(T).name()
            << abort(FatalError);
    }
    return ptr_;
}


template<class T>
inline const T* Foam::autoPtr<T>::operator->() const
{
    return const_cast<autoPtr<T>*>(this)->operator->();
}


template<class T>
inline T& Foam::autoPtr<T>::operator()()
{
    return operator*();
}


template<class T>
inline const T& Foam::autoPtr<T>::operator()() const
{
    return operator*();
}


// ************************************************************************* //
