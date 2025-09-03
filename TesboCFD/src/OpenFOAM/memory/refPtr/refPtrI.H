/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
inline Foam::word Foam::refPtr<T>::typeName()
{
    return Foam::word("refPtr<" + std::string(typeid(T).name()) + '>', false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::refPtr<T>::refPtr() noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(std::nullptr_t) noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(T* p) noexcept
:
    ptr_(p),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::refPtr<T>::refPtr(refPtr<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::refPtr<T>::refPtr(const refPtr<T>& rhs)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (is_pointer())
    {
        if (ptr_)
        {
            rhs.type_ = REF;  // (shallow copy)
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy/move of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::refPtr<T>::refPtr(const refPtr<T>& rhs, bool reuse)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (is_pointer())
    {
        if (ptr_)
        {
            if (reuse)
            {
                rhs.ptr_ = nullptr;
                rhs.type_ = PTR;
            }
            else
            {
                rhs.type_ = REF;  // (shallow copy)
            }
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy/move of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::refPtr<T>::refPtr(std::unique_ptr<T>&& rhs) noexcept
:
    refPtr<T>(rhs.release())
{}


template<class T>
inline Foam::refPtr<T>::refPtr(autoPtr<T>&& rhs) noexcept
:
    refPtr<T>(rhs.release())
{}


template<class T>
inline Foam::refPtr<T>::refPtr(const tmp<T>& rhs, bool reuse)
:
    refPtr<T>()
{
    reset(const_cast<tmp<T>&>(rhs), reuse);
}


template<class T>
inline Foam::refPtr<T>::refPtr(const tmp<T>& rhs)
:
    refPtr<T>()
{
    reset(const_cast<tmp<T>&>(rhs), false);
}


template<class T>
inline Foam::refPtr<T>::refPtr(tmp<T>&& rhs)
:
    refPtr<T>()
{
    reset(const_cast<tmp<T>&>(rhs), true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::refPtr<T>::~refPtr() noexcept
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::refPtr<T>::is_const() const noexcept
{
    return (type_ == CREF);
}


template<class T>
inline bool Foam::refPtr<T>::is_pointer() const noexcept
{
    //OR: return (type_ == PTR);
    return (type_ < REF_Types);
}


template<class T>
inline bool Foam::refPtr<T>::is_reference() const noexcept
{
    //OR: return (type_ == CREF || type_ == REF);
    return (type_ > REF_Types);
}


template<class T>
inline bool Foam::refPtr<T>::movable() const noexcept
{
    return (ptr_ && is_pointer());
}


template<class T>
inline const T& Foam::refPtr<T>::cref() const
{
    if (!ptr_ && is_pointer())
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;  // const reference
}


template<class T>
inline T& Foam::refPtr<T>::ref() const
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempted non-const reference to const object: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_ && is_pointer())
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;  // non-const reference
}


template<class T>
inline Foam::refPtr<T> Foam::refPtr<T>::shallowClone() const noexcept
{
    refPtr<T> dup;

    if (ptr_)
    {
        dup.ptr_ = ptr_;
        dup.type_ = (is_const() ? CREF : REF);
    }

    return dup;
}


template<class T>
inline T* Foam::refPtr<T>::release() noexcept
{
    if (is_pointer())
    {
        T* p = ptr_;
        ptr_ = nullptr;
        return p;
    }

    return nullptr;
}


template<class T>
inline T* Foam::refPtr<T>::ptr() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    if (is_pointer())
    {
        // Release pointer
        T* p = ptr_;
        ptr_ = nullptr;

        return p;
    }

    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::refPtr<T>::clear() const noexcept
{
    if (is_pointer())
    {
        delete ptr_;
        ptr_ = nullptr;
    }
}


template<class T>
inline void Foam::refPtr<T>::reset(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::refPtr<T>::reset(refPtr<T>&& other) noexcept
{
    // Could also make Fatal with FULLDEBUG
    if (&other == this)
    {
        return;  // No self-assignment
    }

    clear();
    ptr_ = other.ptr_;
    type_ = other.type_;

    other.ptr_ = nullptr;
    other.type_ = PTR;
}


template<class T>
inline void Foam::refPtr<T>::reset(tmp<T>& other, bool reuse)
{
    if (other.get())
    {
        if (reuse && other.is_pointer())
        {
            // Acquire pointer.
            // Fatal if pointer is not unique (avoids potential leaks)
            reset(other.ptr());
        }
        else if (other.is_const())
        {
            cref(other.get());
        }
        else
        {
            ref(other.get());
        }
    }
    else
    {
        clear();
    }
}


template<class T>
template<class... Args>
inline T& Foam::refPtr<T>::emplace(Args&&... args)
{
    clear();  // delete old entry
    ptr_ = new T(std::forward<Args>(args)...);
    type_ = PTR;
    return *ptr_;
}


template<class T>
inline void Foam::refPtr<T>::cref(const refPtr<T>& other) noexcept
{
    // Could also make Fatal with FULLDEBUG
    if (&other == this)
    {
        return;  // No self-assignment
    }

    clear();
    ptr_ = other.ptr_;
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::refPtr<T>::cref(const T* p) noexcept
{
    clear();
    ptr_ = const_cast<T*>(p);
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::ref(T& obj) noexcept
{
    clear();
    ptr_ = &obj;
    type_ = REF;
}


template<class T>
inline void Foam::refPtr<T>::ref(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = (ptr_ ? REF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::swap(refPtr<T>& other) noexcept
{
    // Swap is just copy/assign for pointer and enum types
    // Self-swap is effectively ignored
    T* p = ptr_;
    ptr_ = other.ptr_;
    other.ptr_ = p;

    refType t = type_;
    type_ = other.type_;
    other.type_ = t;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::refPtr<T>::operator*() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }
    return *ptr_;
}


template<class T>
inline T& Foam::refPtr<T>::operator*()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T* Foam::refPtr<T>::operator->() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline T* Foam::refPtr<T>::operator->()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::refPtr<T>::operator=(const refPtr<T>& other)
{
    // Could also make Fatal with FULLDEBUG
    if (&other == this)
    {
        return;  // No self-assignment
    }

    clear();

    if (other.is_pointer())
    {
        ptr_ = other.ptr_;
        type_ = PTR;

        other.ptr_ = nullptr;
        other.type_ = PTR;

        if (!ptr_)
        {
            FatalErrorInFunction
                << "Attempted assignment of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempted assignment of an object reference of type "
            << typeid(T).name()
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::refPtr<T>::operator=(tmp<T>&& other)
{
    reset(other, true);  // reuse
}


template<class T>
inline Foam::refPtr<T>::operator tmp<T>()
{
    if (is_pointer())
    {
        return tmp<T>(ptr());
    }
    else
    {
        return tmp<T>(cref());
    }
}


// ************************************************************************* //
