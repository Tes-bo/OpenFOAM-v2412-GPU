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

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline const T* Foam::HashPtrTable<T, Key, Hash>::test(const Key& key) const
{
    // Like lookup() with a nullptr
    const const_iterator iter(this->cfind(key));
    if (iter.good())
    {
        return iter.val();
    }
    return nullptr;
}


template<class T, class Key, class Hash>
inline const T* Foam::HashPtrTable<T, Key, Hash>::get(const Key& key) const
{
    // Like lookup() with a nullptr
    const const_iterator iter(this->cfind(key));
    if (iter.good())
    {
        return iter.val();
    }
    return nullptr;
}


template<class T, class Key, class Hash>
template<class... Args>
inline bool Foam::HashPtrTable<T, Key, Hash>::emplace
(
    const Key& key,
    Args&&... args
)
{
    // Use insertion semantics
    return
    (
        !parent_type::contains(key)
     && this->parent_type::set(key, new T(std::forward<Args>(args)...))
    );
}


template<class T, class Key, class Hash>
template<class... Args>
inline T& Foam::HashPtrTable<T, Key, Hash>::emplaceImpl
(
    const bool overwrite,
    const Key& key,
    Args&&... args
)
{
    T* ptr = nullptr;

    // Replace existing entry unconditionally (overwrite = true)
    // or only if nullptr (overwrite = false)
    iterator iter(this->find(key));
    if (iter.good())
    {
        ptr = iter.val();

        if (overwrite || !ptr)
        {
            // Overwrite existing or replace nullptr
            delete ptr;
            ptr = new T(std::forward<Args>(args)...);
            iter.val() = ptr;
        }
    }
    else
    {
        ptr = new T(std::forward<Args>(args)...);
        this->parent_type::set(key, ptr);
    }

    return *ptr;
}


template<class T, class Key, class Hash>
template<class... Args>
inline T& Foam::HashPtrTable<T, Key, Hash>::emplace_set
(
    const Key& key,
    Args&&... args
)
{
    // overwrite = true
    return this->emplaceImpl(true, key, std::forward<Args>(args)...);
}


template<class T, class Key, class Hash>
template<class... Args>
inline T& Foam::HashPtrTable<T, Key, Hash>::try_emplace
(
    const Key& key,
    Args&&... args
)
{
    // overwrite = false
    return this->emplaceImpl(false, key, std::forward<Args>(args)...);
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    std::unique_ptr<T>&& ptr
)
{
    if (parent_type::insert(key, ptr.get()))
    {
        ptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::insert
(
    const Key& key,
    autoPtr<T>&& ptr
)
{
    if (parent_type::insert(key, ptr.get()))
    {
        ptr.release(); // Now owned by HashPtrTable
        return true;
    }

    return false;
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    T* ptr
)
{
    // Replace existing entry unconditionally
    iterator iter(this->find(key));
    if (iter.good())
    {
        if (iter.val() != ptr)
        {
            delete iter.val();
            iter.val() = ptr;
        }
        return true;
    }

    // Or add new entry
    return this->parent_type::set(key, ptr);
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    std::unique_ptr<T>&& ptr
)
{
    return this->set(key, ptr.release());
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    autoPtr<T>&& ptr
)
{
    return this->set(key, ptr.release());
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    const refPtr<T>& ptr
)
{
    return this->set(key, ptr.ptr());  // release or clone
}


template<class T, class Key, class Hash>
inline bool Foam::HashPtrTable<T, Key, Hash>::set
(
    const Key& key,
    const tmp<T>& ptr
)
{
    return this->set(key, ptr.ptr());  // release or clone
}


// ************************************************************************* //
