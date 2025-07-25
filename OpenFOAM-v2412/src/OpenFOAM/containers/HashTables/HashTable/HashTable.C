/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#ifndef Foam_HashTable_C
#define Foam_HashTable_C

#include "HashTable.H"
#include "List.H"
#include "FixedList.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable() noexcept
:
    HashTableCore(),
    size_(0),
    capacity_(0),
    table_(nullptr)
{}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(const Foam::zero) noexcept
:
    HashTable<T, Key, Hash>()
{}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(const label initialCapacity)
:
    HashTable<T, Key, Hash>()
{
    if (initialCapacity > 0)
    {
        // Like resize() but no initial content to transfer
        capacity_ = HashTableCore::canonicalSize(initialCapacity);
        table_ = new node_type*[capacity_];
        std::fill_n(table_, capacity_, static_cast<node_type*>(nullptr));
    }
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(const HashTable<T, Key, Hash>& ht)
:
    HashTable<T, Key, Hash>(2*ht.size())
{
    for (const_iterator iter = ht.cbegin(); iter != ht.cend(); ++iter)
    {
        insert(iter.key(), iter.val());
    }
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(HashTable<T, Key, Hash>&& rhs) noexcept
:
    HashTableCore(),
    size_(rhs.size_),
    capacity_(rhs.capacity_),
    table_(rhs.table_)
{
    // Stole all contents
    rhs.size_ = 0;
    rhs.capacity_ = 0;
    rhs.table_ = nullptr;
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable
(
    std::initializer_list<std::pair<Key, T>> list,
    const bool overwrite
)
:
    HashTable<T, Key, Hash>()
{
    reserve(list.size());

    for (const auto& keyval : list)
    {
        this->setEntry(overwrite, keyval.first, keyval.second);
    }
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable
(
    const UList<Key>& keys,
    const UList<T>& values,
    const bool overwrite
)
:
    HashTable<T, Key, Hash>()
{
    const label len = std::min(keys.size(), values.size());

    reserve(len);

    for (label i = 0; i < len; ++i)
    {
        this->setEntry(overwrite, keys[i], values[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::~HashTable()
{
    // Remove all entries from table
    clear();

    // Remove the table itself
    capacity_ = 0;
    delete[] table_;
    table_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::toc() const
{
    List<Key> list(size_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        list[count++] = iter.key();
    }

    return list;
}


template<class T, class Key, class Hash>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::sortedToc() const
{
    List<Key> list(this->toc());
    Foam::sort(list);

    return list;
}


template<class T, class Key, class Hash>
template<class Compare>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::sortedToc
(
    const Compare& comp
) const
{
    List<Key> list(this->toc());
    Foam::sort(list, comp);

    return list;
}


template<class T, class Key, class Hash>
Foam::UPtrList<const typename Foam::HashTable<T, Key, Hash>::node_type>
Foam::HashTable<T, Key, Hash>::csorted() const
{
    UPtrList<const node_type> result(size_);

    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        result.set(count++, iter.node());
    }

    Foam::sort(result);

    return result;
}


template<class T, class Key, class Hash>
Foam::UPtrList<typename Foam::HashTable<T, Key, Hash>::node_type>
Foam::HashTable<T, Key, Hash>::sorted()
{
    UPtrList<node_type> result(size_);

    label count = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        result.set(count++, iter.node());
    }

    Foam::sort(result);

    return result;
}



template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocKeys
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    List<Key> list(size_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key()) ? !invert : invert))
        {
            list[count++] = iter.key();
        }
    }

    list.resize(count);
    Foam::sort(list);

    return list;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocValues
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    List<Key> list(size_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.val()) ? !invert : invert))
        {
            list[count++] = iter.key();
        }
    }

    list.resize(count);
    Foam::sort(list);

    return list;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::List<Key> Foam::HashTable<T, Key, Hash>::tocEntries
(
    const BinaryPredicate& pred,
    const bool invert
) const
{
    List<Key> list(size_);
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key(), iter.val()) ? !invert : invert))
        {
            list[count++] = iter.key();
        }
    }

    list.resize(count);
    Foam::sort(list);

    return list;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countKeys
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countValues
(
    const UnaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.val()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::countEntries
(
    const BinaryPredicate& pred,
    const bool invert
) const
{
    label count = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if ((pred(iter.key(), iter.val()) ? !invert : invert))
        {
            ++count;
        }
    }

    return count;
}


template<class T, class Key, class Hash>
template<class... Args>
bool Foam::HashTable<T, Key, Hash>::setEntry
(
    const bool overwrite,
    const Key& key,
    Args&&... args
)
{
    if (!capacity_)
    {
        setCapacity(128);  // Impose an initial capacity
    }

    const label index = hashKeyIndex(key);

    node_type* curr = nullptr;
    node_type* prev = nullptr;

    for (node_type* ep = table_[index]; ep; ep = ep->next_)
    {
        if (key == ep->key())
        {
            curr = ep;
            break;
        }
        prev = ep;
    }

    if (!curr)
    {
        // Not found, insert it at the head
        table_[index] =
            new node_type(table_[index], key, std::forward<Args>(args)...);

        ++size_;
        if (0.8*capacity_ < size_)  // Resize after 0.8 load factor
        {
            if (capacity_ < maxTableSize) setCapacity(2*capacity_);
        }
    }
    else if (overwrite)
    {
        // Overwrite current entry (Perl convention).

        // Can skip if the value is not stored anyhow (Eg, HashSet)
        // - this avoids a useless delete/new
        if (!node_type::stores_value())
        {
            return true;
        }

        node_type* ep = curr->next_;  // next in the linked list

        // In some cases the delete/new could be avoided in favour of move
        // assignment, but cannot be certain that all objects support this
        // or that it behaves the same as a copy construct.

        delete curr;
        ep = new node_type(ep, key, std::forward<Args>(args)...);

        // Replace current element - within list or insert at the head
        if (prev)
        {
            prev->next_ = ep;
        }
        else
        {
            table_[index] = ep;
        }
    }
    else
    {
        // Not overwriting existing entry
        return false;
    }

    return true;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::insert_node(node_type* entry)
{
    if (!entry) return;

    if (!capacity_)
    {
        setCapacity(128);  // Impose an initial capacity
    }

    const label index = hashKeyIndex(entry->key());

    node_type* curr = nullptr;
    //node_type* prev = nullptr;

    for (node_type* ep = table_[index]; ep; ep = ep->next_)
    {
        if (entry->key() == ep->key())
        {
            curr = ep;
            break;
        }
        //prev = ep;
    }

    if (!curr)
    {
        // Not found, insert it at the head
        table_[index] = entry;

        ++size_;
        if (0.8*capacity_ < size_)  // Resize after 80% fill factor
        {
            if (capacity_ < maxTableSize) setCapacity(2*capacity_);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Not inserting " << entry->key() << ": already in table\n"
            << exit(FatalError);
    }
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::erase(const iterator& iter)
{
    // NOTE: we use (const iterator&) here, but treat its contents as mutable.
    //
    // The parameter should be (iterator&), but then the compiler doesn't find
    // it correctly and tries to call as (iterator) instead.

    return iterator_erase(const_cast<iterator&>(iter));
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::erase(const Key& key)
{
    if (size_)
    {
        iterator iter(find(key));
        if (iter.good()) return iterator_erase(iter);
    }
    return false;
}


template<class T, class Key, class Hash>
template<class InputIter>
inline Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    InputIter first,
    InputIter last
)
{
    label changed = 0;

    for
    (
        const label nTotal = this->size();
        changed < nTotal && first != last; // terminate early
        ++first
    )
    {
        if (this->erase(*first))
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    std::initializer_list<Key> keys
)
{
    return erase(keys.begin(), keys.end());
}


template<class T, class Key, class Hash>
template<unsigned N>
inline Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    const FixedList<Key, N>& keys
)
{
    return erase(keys.cbegin(), keys.cend());
}


template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    const UList<Key>& keys
)
{
    return erase(keys.cbegin(), keys.cend());
}


template<class T, class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::label Foam::HashTable<T, Key, Hash>::erase
(
    const HashTable<AnyType, Key, AnyHash>& other
)
{
    const label nTotal = this->size();
    label changed = 0;

    if (other.size() <= nTotal)
    {
        // The other is smaller/same-size, use its keys for removal

        for
        (
            auto iter = other.cbegin();
            changed < nTotal && iter != other.cend(); // Terminate early
            ++iter
        )
        {
            if (erase(iter.key()))
            {
                ++changed;
            }
        }
    }
    else
    {
        // We are smaller: remove if found in the other hash
        for
        (
            iterator iter = begin();
            changed < nTotal && iter != end(); // Terminate early
            ++iter
        )
        {
            if (other.contains(iter.key()) && erase(iter))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::label Foam::HashTable<T, Key, Hash>::retain
(
    const HashTable<AnyType, Key, AnyHash>& other
)
{
    const label nTotal = this->size();
    label changed = 0;

    if (other.empty())
    {
        // Trivial case
        changed = nTotal;
        this->clear();
    }
    else
    {
        // Inverted logic: remove if *not* found in the other hash

        for (iterator iter = begin(); iter != end(); ++iter)
        {
            if (!other.contains(iter.key()) && erase(iter))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::setCapacity(label newCapacity)
{
    newCapacity = HashTableCore::canonicalSize(newCapacity);

    if (newCapacity == capacity_)
    {
        return;
    }

    if (!size_)
    {
        // Table is unpopulated - can already remove now
        capacity_ = 0;
        delete[] table_;
        table_ = nullptr;
    }

    if (!newCapacity)
    {
        // Special handling for capacity = 0.
        if (size_)
        {
            WarningInFunction
                << "HashTable contains " << size_
                << " elements, cannot set capacity to 0 buckets!" << nl;
        }
        return;
    }

    // Swap primary table entries: size_ is left untouched

    auto oldTable = table_;
    const label oldCapacity = capacity_;

    capacity_ = newCapacity;
    table_ = new node_type*[capacity_];
    std::fill_n(table_, capacity_, static_cast<node_type*>(nullptr));

    if (!oldTable)
    {
        return;
    }

    // Move to new table[] but with new chaining.

    for (label i = 0, pending = size_; pending && i < oldCapacity; ++i)
    {
        for (node_type* ep = oldTable[i]; ep; /*nil*/)
        {
            node_type* next = ep->next_;

            // Move to new location
            {
                const label newIdx = hashKeyIndex(ep->key());

                ep->next_ = table_[newIdx];  // add to head
                table_[newIdx] = ep;
            }

            ep = next;      // continue in the linked-list
            --pending;      // note any early completion
        }
        oldTable[i] = nullptr;
    }

    delete[] oldTable;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::resize(label newCapacity)
{
    setCapacity(newCapacity);
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::reserve(label numEntries)
{
    if (numEntries > size_)
    {
        // From number of entries to estimated capacity
        // - initial load factor of 0.5
        numEntries *= 2;
        if (numEntries > capacity_) setCapacity(numEntries);
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::clear()
{
    if (!table_)
    {
        capacity_ = 0;  // Paranoid
    }

    for (label i = 0, pending = size_; pending && i < capacity_; ++i)
    {
        for (node_type* ep = table_[i]; ep; /*nil*/)
        {
            node_type* next = ep->next_;

            delete ep;

            ep = next;  // continue in the linked-list
            --pending;  // note any early completion
        }
        table_[i] = nullptr;
    }
    size_ = 0;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::clearStorage()
{
    // Remove all entries from table
    clear();

    // Remove the table itself
    capacity_ = 0;
    delete[] table_;
    table_ = nullptr;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::swap(HashTable<T, Key, Hash>& rhs) noexcept
{
    if (this == &rhs)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(size_, rhs.size_);
    std::swap(capacity_, rhs.capacity_);
    std::swap(table_, rhs.table_);
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::transfer(HashTable<T, Key, Hash>& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    swap(rhs);
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterKeys
(
    const UnaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.key()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class UnaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterValues
(
    const UnaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.val()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
template<class BinaryPredicate>
Foam::label Foam::HashTable<T, Key, Hash>::filterEntries
(
    const BinaryPredicate& pred,
    const bool pruning
)
{
    label changed = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        // Matches? either prune (pruning) or keep (!pruning)
        if
        (
            (pred(iter.key(), iter.val()) ? pruning : !pruning)
         && erase(iter)
        )
        {
            ++changed;
        }
    }

    return changed;
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::merge(HashTable<T, Key, Hash>& source)
{
    // Self-merge implicitly a no-op

    if (node_type::stores_value())
    {
        // key/val pair
        for (iterator iter = source.begin(); iter != source.end(); ++iter)
        {
            if (!contains(iter.key()))
            {
                node_type* entry = source.iterator_extract(iter);
                this->insert_node(entry);
            }
        }
    }
    else
    {
        // key only, does not store a value.
        // Since the key is const, juggling the node does not work

        for (iterator iter = source.begin(); iter != source.end(); ++iter)
        {
            if (emplace(iter.key()))
            {
                source.erase(iter);
            }
        }
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::merge(HashTable<T, Key, Hash>&& source)
{
    merge(source);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    const HashTable<T, Key, Hash>& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    this->clear();
    this->reserve(rhs.size());

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        insert(iter.key(), iter.val());
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    std::initializer_list<std::pair<Key, T>> rhs
)
{
    this->clear();
    this->reserve(rhs.size());

    for (const auto& keyval : rhs)
    {
        set(keyval.first, keyval.second);
    }
}


template<class T, class Key, class Hash>
void Foam::HashTable<T, Key, Hash>::operator=
(
    HashTable<T, Key, Hash>&& rhs
)
{
    transfer(rhs);  // Includes self-assignment check
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::operator==
(
    const HashTable<T, Key, Hash>& rhs
) const
{
    // Sizes (number of keys) must match
    if (size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        const const_iterator other(this->cfind(iter.key()));

        if (!other.good() || other.val() != iter.val())
        {
            return false;
        }
    }

    return true;
}


template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::operator!=
(
    const HashTable<T, Key, Hash>& rhs
) const
{
    return !operator==(rhs);
}


template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>& Foam::HashTable<T, Key, Hash>::operator+=
(
    const HashTable<T, Key, Hash>& rhs
)
{
    // Avoid no-ops:
    if (rhs.size() && (this != &rhs))
    {
        if (this->size())
        {
            for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
            {
                insert(iter.key(), iter.val());
            }
        }
        else
        {
            (*this) = rhs;
        }
    }

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Iterators, Friend Operators

#include "HashTableIter.C"
#include "HashTableIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
