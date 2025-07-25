/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
inline Foam::word Foam::tmp<T>::typeName()
{
    return Foam::word("tmp<" + std::string(typeid(T).name()) + '>', false);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::tmp<T>::checkUseCount() const
{
    if (ptr_ && ptr_->refCount::use_count() > 1)
    {
        FatalErrorInFunction
            << "Attempt to create more than "
            << (ptr_->refCount::use_count() + 1)
            << " tmp's referring to the same object of type tmp<"
            << typeid(T).name() << '>'
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::tmp<T>::tmp() noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::tmp<T>::tmp(std::nullptr_t) noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline Foam::tmp<T>::tmp(T* p)
:
    ptr_(p),
    type_(PTR)
{
    if (ptr_ && !ptr_->refCount::unique())
    {
        FatalErrorInFunction
            << "Attempted construction of a "
            << this->typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline constexpr Foam::tmp<T>::tmp(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::tmp<T>::tmp(tmp<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& rhs)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (is_pointer())
    {
        if (ptr_)
        {
            ptr_->refCount::operator++();
            this->checkUseCount();
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
inline Foam::tmp<T>::tmp(const tmp<T>& rhs, bool reuse)
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
                ptr_->refCount::operator++();
                this->checkUseCount();
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
inline Foam::tmp<T>::tmp(autoPtr<T>&& rhs) noexcept
:
    tmp<T>(rhs.release())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::tmp<T>::~tmp() noexcept
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::is_const() const noexcept
{
    return (type_ == CREF);
}


template<class T>
inline bool Foam::tmp<T>::is_pointer() const noexcept
{
    //OR: return (type_ == PTR || type_ == CACHE_PTR);
    return (type_ < REF_Types);
}


template<class T>
inline bool Foam::tmp<T>::is_reference() const noexcept
{
    //OR: return (type_ == CREF || type_ == REF);
    return (type_ > REF_Types);
}


template<class T>
inline bool Foam::tmp<T>::movable() const noexcept
{
    return (ptr_ && type_ == PTR && ptr_->refCount::unique());
}


template<class T>
inline const T& Foam::tmp<T>::cref() const
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
inline T& Foam::tmp<T>::ref() const
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
inline T* Foam::tmp<T>::ptr() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    if (type_ == PTR)
    {
        if (!ptr_->refCount::unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type "
                << this->typeName()
                << abort(FatalError);
        }

        // Release pointer
        T* p = ptr_;
        ptr_ = nullptr;

        return p;
    }

    // CACHE_PTR (immovable) is cloned, as per a reference
    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::tmp<T>::clear() const noexcept
{
    if (ptr_ && is_pointer())
    {
        if (ptr_->refCount::unique())
        {
            delete ptr_;
        }
        else
        {
            ptr_->refCount::operator--();
        }
        ptr_ = nullptr;
    }
}


template<class T>
inline void Foam::tmp<T>::protect(bool on) noexcept
{
    if (on)
    {
        // ON: from PTR -> CACHE_PTR, but not nullptr
        if (ptr_ && type_ == PTR)
        {
            type_ = CACHE_PTR;
        }
    }
    else
    {
        // OFF: from CACHE_PTR -> PTR
        if (type_ == CACHE_PTR)
        {
            type_ = PTR;
        }
    }
}


template<class T>
inline void Foam::tmp<T>::reset(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::tmp<T>::reset(tmp<T>&& other) noexcept
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
template<class... Args>
inline T& Foam::tmp<T>::emplace(Args&&... args)
{
    clear();  // delete old entry
    ptr_ = new T(std::forward<Args>(args)...);
    type_ = PTR;
    return *ptr_;
}


template<class T>
inline void Foam::tmp<T>::cref(const tmp<T>& other) noexcept
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
inline void Foam::tmp<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::tmp<T>::cref(const T* p) noexcept
{
    clear();
    ptr_ = const_cast<T*>(p);
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::tmp<T>::ref(T& obj) noexcept
{
    clear();
    ptr_ = &obj;
    type_ = REF;
}


template<class T>
inline void Foam::tmp<T>::ref(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = (ptr_ ? REF : PTR);
}


template<class T>
inline void Foam::tmp<T>::swap(tmp<T>& other) noexcept
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
inline const T* Foam::tmp<T>::operator->() const
{
    if (!ptr_ && is_pointer())
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_ && is_pointer())
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& other)
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
        type_ = other.type_;

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
inline void Foam::tmp<T>::operator=(tmp<T>&& other) noexcept
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
inline void Foam::tmp<T>::operator=(T* p)
{
    if (!p)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!p->refCount::unique())
    {
        FatalErrorInFunction
            << "Attempted assignment of a "
            << this->typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    reset(p);
}


// ************************************************************************* //
