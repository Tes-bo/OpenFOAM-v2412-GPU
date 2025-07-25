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

#ifndef HashSet_C
#define HashSet_C

#include "HashSet.H"
#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashSet<Key, Hash>::assignMany
(
    const label nItems,
    InputIter first,
    InputIter last
)
{
    this->clear();
    this->reserve(nItems);

    return insert(first, last);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Key, class Hash>
template<unsigned N>
Foam::HashSet<Key, Hash>::HashSet(const FixedList<Key, N>& list)
:
    parent_type(2*list.size())
{
    insert(list.begin(), list.end());
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>::HashSet(const UList<Key>& list)
:
    parent_type(2*list.size())
{
    insert(list.begin(), list.end());
}


template<class Key, class Hash>
template<class Addr>
Foam::HashSet<Key, Hash>::HashSet(const IndirectListBase<Key, Addr>& list)
:
    parent_type(2*list.size())
{
    insert(list.begin(), list.end());
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>::HashSet(std::initializer_list<Key> list)
:
    parent_type(2*list.size())
{
    insert(list.begin(), list.end());
}


template<class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::HashSet<Key, Hash>::HashSet
(
    const HashTable<AnyType, Key, AnyHash>& tbl
)
:
    parent_type(2*tbl.size())
{
    for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
    {
        this->insert(iter.key());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashSet<Key, Hash>::insert
(
    InputIter first,
    InputIter last
)
{
    label changed = 0;
    while (first != last)
    {
        if (insert(*first))
        {
            ++changed;
        }
        ++first;
    }
    return changed;
}


template<class Key, class Hash>
inline Foam::label Foam::HashSet<Key, Hash>::insert
(
    std::initializer_list<Key> list
)
{
    return insert(list.begin(), list.end());
}


template<class Key, class Hash>
template<unsigned N>
inline Foam::label Foam::HashSet<Key, Hash>::insert
(
    const FixedList<Key, N>& list
)
{
    return insert(list.begin(), list.end());
}


template<class Key, class Hash>
inline Foam::label Foam::HashSet<Key, Hash>::insert
(
    const UList<Key>& list
)
{
    return insert(list.begin(), list.end());
}


template<class Key, class Hash>
template<class Addr>
inline  Foam::label Foam::HashSet<Key, Hash>::insert
(
    const IndirectListBase<Key, Addr>& list
)
{
    return insert(list.begin(), list.end());
}


template<class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashSet<Key, Hash>::unset
(
    InputIter first,
    InputIter last
)
{
    return this->parent_type::erase(first, last);
}


template<class Key, class Hash>
inline Foam::label Foam::HashSet<Key, Hash>::unset
(
    std::initializer_list<Key> list
)
{
    return unset(list.begin(), list.end());
}


template<class Key, class Hash>
template<unsigned N>
inline Foam::label Foam::HashSet<Key, Hash>::unset
(
    const FixedList<Key, N>& list
)
{
    return unset(list.begin(), list.end());
}


template<class Key, class Hash>
inline Foam::label Foam::HashSet<Key, Hash>::unset
(
    const UList<Key>& list
)
{
    return unset(list.begin(), list.end());
}


template<class Key, class Hash>
template<class Addr>
inline Foam::label Foam::HashSet<Key, Hash>::unset
(
    const IndirectListBase<Key, Addr>& list
)
{
    return unset(list.begin(), list.end());
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::merge(HashSet<Key, Hash>& source)
{
    parent_type::merge(source);
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::merge(HashSet<Key, Hash>&& source)
{
    parent_type::merge(source);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Key, class Hash>
inline bool Foam::HashSet<Key, Hash>::operator()(const Key& key) const noexcept
{
    return this->contains(key);
}


template<class Key, class Hash>
inline bool Foam::HashSet<Key, Hash>::operator[](const Key& key) const noexcept
{
    return this->contains(key);
}


template<class Key, class Hash>
template<unsigned N>
void Foam::HashSet<Key, Hash>::operator=(const FixedList<Key, N>& rhs)
{
    assignMany(rhs.size(), rhs.begin(), rhs.end());
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator=(const UList<Key>& rhs)
{
    assignMany(rhs.size(), rhs.begin(), rhs.end());
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator=(std::initializer_list<Key> rhs)
{
    assignMany(rhs.size(), rhs.begin(), rhs.end());
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator==(const HashSet<Key, Hash>& rhs) const
{
    // Trivial checks first
    if (this == &rhs)
    {
        return true;
    }
    else if (this->size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (!this->contains(iter.key()))
        {
            return false;
        }
    }

    return true;
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator!=(const HashSet<Key, Hash>& rhs) const
{
    return !operator==(rhs);
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>&
Foam::HashSet<Key, Hash>::operator|=(const HashSet<Key, Hash>& rhs)
{
    // Add rhs elements into lhs - avoid no-ops
    if (this != &rhs)
    {
        for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
        {
            this->insert(iter.key());
        }
    }

    return *this;
}


template<class Key, class Hash>
inline Foam::HashSet<Key, Hash>&
Foam::HashSet<Key, Hash>::operator&=(const HashSet<Key, Hash>& rhs)
{
    this->parent_type::retain(rhs);
    return *this;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>&
Foam::HashSet<Key, Hash>::operator^=(const HashSet<Key, Hash>& rhs)
{
    // Add missed rhs elements, remove duplicate elements
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (this->contains(iter.key()))
        {
            this->erase(iter.key());
        }
        else
        {
            this->insert(iter.key());
        }
    }

    return *this;
}


template<class Key, class Hash>
inline Foam::HashSet<Key, Hash>&
Foam::HashSet<Key, Hash>::operator+=(const HashSet<Key, Hash>& rhs)
{
    return this->operator|=(rhs);
}


template<class Key, class Hash>
inline Foam::HashSet<Key, Hash>&
Foam::HashSet<Key, Hash>::operator-=(const HashSet<Key, Hash>& rhs)
{
    this->parent_type::erase(rhs);

    return *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Key, class Hash>
Foam::Ostream& Foam::operator<<(Ostream& os, const HashSet<Key, Hash>& rhs)
{
    return rhs.writeKeys(os, Detail::ListPolicy::short_length<Key>::value);
}


/* * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * */

template<class Key, class Hash>
Foam::HashSet<Key, Hash> Foam::operator|
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
)
{
    HashSet<Key, Hash> result(a);
    result |= b;
    return result;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash> Foam::operator&
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
)
{
    HashSet<Key, Hash> result(a.capacity());

    for (const Key& k : a)
    {
        if (b.contains(k))
        {
            result.insert(k);
        }
    }

    return result;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash> Foam::operator^
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
)
{
    HashSet<Key, Hash> result(a);
    result ^= b;
    return result;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash> Foam::operator-
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
)
{
    HashSet<Key, Hash> result(a.capacity());

    for (const Key& k : a)
    {
        if (!b.contains(k))
        {
            result.insert(k);
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::iterator
Foam::HashSet<Key, Hash>::begin()
{
    return iterator
    (
        static_cast<parent_type&>(*this).begin()
    );
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::begin() const
{
    return const_iterator
    (
        static_cast<const parent_type&>(*this).begin()
    );
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::cbegin() const
{
    return const_iterator
    (
        static_cast<const parent_type&>(*this).cbegin()
    );
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::iterator
Foam::HashSet<Key, Hash>::end() noexcept
{
    return iterator();
}


template<class Key, class Hash>
inline typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::end() const noexcept
{
    return const_iterator();
}


template<class Key, class Hash>
inline constexpr typename Foam::HashSet<Key, Hash>::const_iterator
Foam::HashSet<Key, Hash>::cend() const noexcept
{
    return const_iterator();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
