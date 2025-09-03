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

#include "error.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::HashTable<T, Key, Hash>::hashKeyIndex(const Key& key) const
{
    // capacity is always a power of two - this is the modulus
    return Hash()(key) & (capacity_ - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::at(const Key& key)
{
    iterator iter(this->find(key));

    if (!iter.good())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return iter.val();
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::at(const Key& key) const
{
    const const_iterator iter(this->cfind(key));

    if (!iter.good())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return iter.val();
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::contains(const Key& key) const
{
    if (size_)
    {
        return Iterator<true>(this, key).good();
    }

    return false;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::find
(
    const Key& key
)
{
    if (size_)
    {
        return iterator(Iterator<false>(this, key));
    }

    return iterator();
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    return this->cfind(key);
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cfind
(
    const Key& key
) const
{
    if (size_)
    {
        return const_iterator(Iterator<true>(this, key));
    }

    return const_iterator();
}


template<class T, class Key, class Hash>
template<class... Args>
inline bool Foam::HashTable<T, Key, Hash>::emplace
(
    const Key& key,
    Args&&... args
)
{
    return this->setEntry(false, key, std::forward<Args>(args)...);
}


template<class T, class Key, class Hash>
template<class... Args>
inline bool Foam::HashTable<T, Key, Hash>::emplace_set
(
    const Key& key,
    Args&&... args
)
{
    return this->setEntry(true, key, std::forward<Args>(args)...);
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& obj
)
{
    return this->setEntry(false, key, obj);
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    T&& obj
)
{
    return this->setEntry(false, key, std::forward<T>(obj));
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& obj
)
{
    return this->setEntry(true, key, obj);   // Overwrite
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    T&& obj
)
{
    return this->setEntry(true, key, std::forward<T>(obj));  // Overwrite
}


// FUTURE?
// If 'key' already exists, assign \c std::forward<T>(obj) to the
// entry. If it does not exist, uses insert()
//
// template<class T, class Key, class Hash>
// inline void Foam::HashTable<T, Key, Hash>::insert_or_assign
// (
//     const Key& key,
//     T&& obj
// )
// {
//     iterator iter(this->find(key));
//
//     if (iter.good())
//     {
//         iter.val() = std::forward<T>(obj);
//     }
//     else
//     {
//         this->setEntry(true, key, std::forward<T>(obj));
//         // iter = this->find(key);
//     }
//     // return iter.val();
// }


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::lookup
(
    const Key& key,
    const T& deflt
) const
{
    const const_iterator iter(this->cfind(key));
    return iter.good() ? iter.val() : deflt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key)
{
    return this->at(key);
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    return this->at(key);
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()(const Key& key)
{
    iterator iter(this->find(key));

    if (iter.good())
    {
        return iter.val();
    }

    this->emplace(key);
    return find(key).val();
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()
(
    const Key& key,
    const T& deflt
)
{
    iterator iter(this->find(key));

    if (iter.good())
    {
        return iter.val();
    }

    this->insert(key, deflt);
    return find(key).val();
}


// ************************************************************************* //
