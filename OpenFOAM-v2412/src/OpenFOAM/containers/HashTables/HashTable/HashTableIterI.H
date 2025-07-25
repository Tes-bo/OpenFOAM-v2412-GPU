/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * iterator base * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<bool Const>
inline constexpr
Foam::HashTable<T, Key, Hash>::Iterator<Const>::Iterator() noexcept
:
    entry_(nullptr),
    container_(nullptr),
    index_(0)
{}


template<class T, class Key, class Hash>
template<bool Const>
inline Foam::HashTable<T, Key, Hash>::Iterator<Const>::Iterator
(
    table_type* tbl
)
:
    entry_(nullptr),
    container_(tbl),
    index_(0)
{
    if (container_ && container_->size_)
    {
        // Locate the first non-nullptr table entry
        while
        (
            !(entry_ = container_->table_[index_])
         && ++index_ < container_->capacity_
        )
        {}

        if (index_ >= container_->capacity_)
        {
            // Nothing found - make it an end iterator
            entry_ = nullptr;
            index_ = 0;
        }
    }
}


//
// Any changes here may need changes in iterator_erase(), iterator_extract()
// methods too
//
template<class T, class Key, class Hash>
template<bool Const>
inline void
Foam::HashTable<T, Key, Hash>::Iterator<Const>::increment()
{
    if (index_ < 0)
    {
        // Negative index is a special value from erase
        //
        // Saved as (-index-1), retrieved as (-(index-1)) but to with an
        // extra (-1) to compensate for the ++ in the following while loop
        index_ = -(index_+1) - 1;
    }
    else if (index_ < container_->capacity_ && entry_ && entry_->next_)
    {
        // Move to next element on the linked-list
        entry_ = entry_->next_;
        return;
    }

    // Move to the next non-nullptr table entry
    while
    (
        ++index_ < container_->capacity_
     && !(entry_ = container_->table_[index_])
    )
    {}

    if (index_ >= container_->capacity_)
    {
        // Nothing found - make it an end iterator
        entry_ = nullptr;
        index_ = 0;
    }
}


template<class T, class Key, class Hash>
template<bool Const>
inline Foam::Ostream& Foam::HashTable<T, Key, Hash>::Iterator<Const>::print
(
    Ostream& os
) const
{
    if (entry_)
    {
        entry_->print(os);
    }
    return os;
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator&
Foam::HashTable<T, Key, Hash>::iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::iterator::operator++(int)
{
    iterator iter(*this);
    this->increment();
    return iter;
}


// * * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::const_iterator::operator++()
{
    this->increment();
    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::const_iterator::operator++(int)
{
    const_iterator iter(*this);
    this->increment();
    return iter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::begin()
{
    return iterator(Iterator<false>(this));
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::begin() const
{
    return const_iterator(Iterator<true>(this));
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cbegin() const
{
    return const_iterator(Iterator<true>(this));
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::end() noexcept
{
    return iterator();
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::end() const noexcept
{
    return const_iterator();
}


template<class T, class Key, class Hash>
inline constexpr typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cend() const noexcept
{
    return const_iterator();
}


// ************************************************************************* //
