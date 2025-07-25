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

#include "HashPtrTable.H"
#include "autoPtr.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable
(
    const HashPtrTable<T, Key, Hash>& rhs
)
:
    parent_type(rhs.capacity())
{
    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        const Key& k = iter.key();
        const T* ptr = iter.val();

        if (ptr)
        {
            this->set(k, new T(*ptr));
        }
        else
        {
            this->set(k, nullptr);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::~HashPtrTable()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::release(iterator& iter)
{
    if (iter.good())
    {
        autoPtr<T> old(iter.val());
        iter.val() = nullptr;
        return old;
    }

    return nullptr;
}


template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::release(const Key& key)
{
    iterator iter(this->find(key));
    return this->release(iter);
}


template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::remove(iterator& iter)
{
    if (iter.good())
    {
        autoPtr<T> old(iter.val());
        this->parent_type::erase(iter);
        return old;
    }

    return nullptr;
}


template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::remove(const Key& key)
{
    iterator iter(this->find(key));
    return this->remove(iter);
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(iterator& iter)
{
    if (iter.good())
    {
        T* ptr = iter.val();

        if (this->parent_type::erase(iter))
        {
            delete ptr;
            return true;
        }
    }

    return false;
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(const Key& key)
{
    iterator iter(this->find(key));
    return this->erase(iter);
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::clear()
{
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        delete iter.val();
    }

    this->parent_type::clear();
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::merge
(
    HashPtrTable<T, Key, Hash>& source
)
{
    // Use parent merge (of raw pointer entries)
    // and by-pass any pointer deletions etc.
    parent_type::merge(static_cast<parent_type&>(source));
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::merge
(
    HashPtrTable<T, Key, Hash>&& source
)
{
    // Use parent merge (of raw pointer entries)
    // and by-pass any pointer deletions etc.
    parent_type::merge(static_cast<parent_type&>(source));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::operator=
(
    const HashPtrTable<T, Key, Hash>& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->clear();

    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        const Key& k = iter.key();
        const T* ptr = iter.val();

        if (ptr)
        {
            this->set(k, new T(*ptr));
        }
        else
        {
            this->set(k, nullptr);
        }
    }
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::operator=
(
    HashPtrTable<T, Key, Hash>&& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->clear();
    this->transfer(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashPtrTableIO.C"

// ************************************************************************* //
