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

Class
    Foam::HashTable

Description
    A HashTable similar to \c std::unordered_map.

    The entries are considered \a unordered since their placement
    depends on the method used to generate the hash key index, the
    table capacity, insertion order etc. When the key order is
    important, use the sortedToc() method to obtain a list of sorted
    keys and use that for further access, or the csorted()/sorted() methods
    to obtain a UPtrList of entries to traverse in sorted order.

    Internally the table uses closed addressing into a flat storage space
    with collisions handled by linked-list chaining.

    - The max_load_factor is 0.8, but a load factor 0.5 is generally
      assumed when initially creating a HashTable (ie, use an capacity
      of twice the expected number elements).
    - When inserting into a table without an initial capacity,
      a default capacity (bucket count) of 128 is used.

    The end iterator of all hash-tables has a nullptr to the hash entry.
    Thus avoid separate allocation for each table and use a single one with
    a nullptr. The hash-table iterators always have an entry-pointer as the
    first member data, which allows reinterpret_cast from anything else with
    a nullptr as its first data member.
    The nullObject is such an item (with a nullptr data member).

Note
    For historical reasons, dereferencing the table iterator
    (eg, \a *iter) returns a reference to the stored object value
    rather than the stored key/val pair like std::unordered_map does.

    The HashTable iterator:
    \code
        forAllConstIters(table, iter)
        {
            Info<< "val:" << *iter << nl
                << "key:" << iter.key() << nl
                << "val:" << iter.val() << nl;
        }
    \endcode
    whereas for the \c std::unordered_map iterator:
    \code
        forAllConstIters(stdmap, iter)
        {
            Info<< "key/val:" << *iter << nl
                << "key:" << iter->first << nl
                << "val:" << iter->second << nl;
        }
    \endcode
    This difference is most evident when using range-for syntax.

SourceFiles
    HashTableI.H
    HashTableIterI.H
    HashTable.C
    HashTableIO.C
    HashTableIter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_HashTable_H
#define Foam_HashTable_H

#include "stdFoam.H"
#include "word.H"
#include "zero.H"
#include "Hash.H"
#include "HashTableDetail.H"
#include "HashTableCore.H"

#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class T> class List;
template<class T> class UList;
template<class T> class UPtrList;
template<class T, unsigned N> class FixedList;
template<class T, class Key, class Hash> class HashTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream&, HashTable<T, Key, Hash>&);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream&, const HashTable<T, Key, Hash>&);

/*---------------------------------------------------------------------------*\
                          Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=Foam::Hash<Key>>
class HashTable
:
    public HashTableCore
{
public:

    // Data Types

        //- The template instance used for this HashTable
        typedef HashTable<T, Key, Hash> this_type;

        //- A table entry (node) that encapsulates the key/val tuple
        //- with an additional linked-list entry for hash collisions
        typedef typename std::conditional
        <
            std::is_same<Foam::zero, typename std::remove_cv<T>::type>::value,
            Detail::HashTableSingle<Key>,
            Detail::HashTablePair<Key, T>
        >::type node_type;


    // STL type definitions

        //- The second template parameter, type of keys used.
        typedef Key key_type;

        //- The first template parameter, type of objects contained.
        typedef T mapped_type;

        //- Same as mapped_type for OpenFOAM HashTables
        //  Note that this is different than the std::map definition.
        typedef T value_type;

        //- The third template parameter, the hash index method.
        typedef Hash hasher;

        //- Pointer type for storing into value_type objects.
        //  This type is usually 'value_type*'.
        typedef T* pointer;

        //- Reference to the stored value_type.
        //  This type is usually 'value_type&'.
        typedef T& reference;

        //- Const pointer type for the stored value_type.
        typedef const T* const_pointer;

        //- Const reference to the stored value_type.
        typedef const T& const_reference;

        //- The type to represent the difference between two iterators
        typedef label difference_type;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;


        //- Forward iterator with non-const access
        class iterator;

        //- Forward iterator with const access
        class const_iterator;


private:

    // Private Data

        //- The number of elements currently stored in table
        label size_;

        //- Number of buckets allocated in table
        label capacity_;

        //- The table of primary nodes
        node_type** table_;


    // Private Member Functions

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key& key) const;

        //- Assign a new hash-entry to a possibly already existing key.
        //  \return True if the new entry was set.
        template<class... Args>
        bool setEntry(const bool overwrite, const Key& key, Args&&... args);

        //- Insert node (low-level). The node must not previously exist!
        void insert_node(node_type* entry);

        //- Low-level entry erasure using iterator internals.
        //  This invalidates the iterator until the next ++ operation.
        //  \return True if the corresponding entry existed and was removed
        bool iterator_erase(iterator& iter);

        //- Unlink node from iterator (low-level).
        node_type* iterator_extract(iterator& iter);

        //- Read hash table
        Istream& readTable(Istream& is);

        //- Write hash table
        Ostream& writeTable(Ostream& os) const;


public:

    // Constructors

        //- Default construct: empty without allocation (capacity=0)
        HashTable() noexcept;

        //- Construct empty without allocation (capacity=0)
        explicit HashTable(const Foam::zero) noexcept;

        //- Construct empty with initial table capacity
        explicit HashTable(const label initialCapacity);

        //- Construct from Istream
        HashTable(Istream& is);

        //- Copy construct
        HashTable(const this_type& ht);

        //- Move construct
        HashTable(this_type&& rhs) noexcept;

        //- Construct from key/value pairs in initializer list
        //  By default, uses insert not overwrite semantics for duplicates.
        HashTable
        (
            std::initializer_list<std::pair<Key, T>> list,
            const bool overwrite = false
        );

        //- Construct from key/value pairs
        //  By default, uses insert not overwrite semantics for duplicates.
        HashTable
        (
            const UList<Key>& keys,
            const UList<T>& values,
            const bool overwrite = false
        );


    //- Destructor
    ~HashTable();


    // Member Functions

    // Capacity

        //- True if the hash table is empty
        bool empty() const noexcept { return !size_; }

        //- The number of elements in table
        label size() const noexcept { return size_; }

        //- The size of the underlying table (the number of buckets)
        label capacity() const noexcept { return capacity_; }


    // Access

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline T& at(const Key& key);

        //- Find and return a hashed entry. FatalError if it does not exist.
        inline const T& at(const Key& key) const;

        //- True if hashed key is contained (found) in table
        inline bool contains(const Key& key) const;

        //- Find and return an iterator set at the hashed entry
        //  If not found iterator = end()
        inline iterator find(const Key& key);

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        inline const_iterator find(const Key& key) const;

        //- Find and return an const_iterator set at the hashed entry
        //  If not found iterator = end()
        inline const_iterator cfind(const Key& key) const;

        //- Return hashed entry if it exists, or return the given default
        inline const T& lookup(const Key& key, const T& deflt) const;


    // Table of contents

        //- The table of contents (the keys) in unsorted order.
        List<Key> toc() const;

        //- The table of contents (the keys) in sorted order
        List<Key> sortedToc() const;

        //- The table of contents (the keys) sorted according to the
        //- specified comparator
        template<class Compare>
        List<Key> sortedToc(const Compare& comp) const;

        //- The table of contents (the keys) selected according to the
        //- unary predicate applied to the \b keys.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class UnaryPredicate>
        List<Key> tocKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- The table of contents (the keys) selected according to the
        //- unary predicate applied to the \b values.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class UnaryPredicate>
        List<Key> tocValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- The table of contents (the keys) selected according to the
        //- binary predicate applied to the \b keys and \b values.
        //  \param invert changes the logic to select when the predicate
        //      is false
        //  \return sorted list of selected keys
        template<class BinaryPredicate>
        List<Key> tocEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


    // Sorted entries

        //- Const access to the hash-table contents in sorted order
        //- (sorted by keys).
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<const node_type> csorted() const;

        //- Non-const access to the hash-table contents in sorted order
        //- (sorted by keys).
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<node_type> sorted();


    // Counting

        //- Count the number of keys that satisfy the unary predicate
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class UnaryPredicate>
        label countKeys
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of values that satisfy the unary predicate
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class UnaryPredicate>
        label countValues
        (
            const UnaryPredicate& pred,
            const bool invert = false
        ) const;

        //- Count the number of entries that satisfy the binary predicate.
        //  \param invert changes the logic to select when the predicate
        //      is false
        template<class BinaryPredicate>
        label countEntries
        (
            const BinaryPredicate& pred,
            const bool invert = false
        ) const;


    // Edit

        //- Emplace insert a new entry, not overwriting existing entries.
        //  \return True if the entry did not previously exist in the table.
        template<class... Args>
        inline bool emplace(const Key& key, Args&&... args);

        //- Emplace set an entry, overwriting any existing entries.
        //  \return True, since it always overwrites any entries.
        template<class... Args>
        inline bool emplace_set(const Key& key, Args&&... args);

        //- Copy insert a new entry, not overwriting existing entries.
        //  \return True if the entry did not previously exist in the table.
        inline bool insert(const Key& key, const T& obj);

        //- Move insert a new entry, not overwriting existing entries.
        //  \return True if the entry did not previously exist in the table.
        inline bool insert(const Key& key, T&& obj);

        //- Copy assign a new entry, overwriting existing entries.
        //  \return True, since it always overwrites any entries.
        inline bool set(const Key& key, const T& obj);

        //- Move assign a new entry, overwriting existing entries.
        //  \return True, since it always overwrites any entries.
        inline bool set(const Key& key, T&& obj);

        //- Erase an entry specified by given iterator
        //  This invalidates the iterator until the next ++ operation.
        //
        //  Includes a safeguard against the end-iterator such that the
        //  following is safe:
        //  \code
        //      auto iter = table.find(unknownKey);
        //      table.erase(iter);
        //  \endcode
        //  which is what \code table.erase(unknownKey) \endcode does anyhow.
        //
        //  \return True if the corresponding entry existed and was removed
        bool erase(const iterator& iter);

        //- Erase an entry specified by the given key
        //  \return True if the entry existed and was removed
        bool erase(const Key& key);

        //- Remove table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items removed
        template<class AnyType, class AnyHash>
        label erase(const HashTable<AnyType, Key, AnyHash>& other);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        inline label erase(std::initializer_list<Key> keys);

        //- Remove multiple entries using an iterator range of keys
        template<class InputIter>
        inline label erase(InputIter first, InputIter last);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        template<unsigned N>
        inline label erase(const FixedList<Key, N>& keys);

        //- Remove table entries given by the listed keys
        //  \return The number of items removed
        inline label erase(const UList<Key>& keys);

        //- Retain table entries given by keys of the other hash-table.
        //
        //  The other hash-table must have the same type of key, but the
        //  type of values held and the hashing function are arbitrary.
        //
        //  \return The number of items changed (removed)
        template<class AnyType, class AnyHash>
        label retain(const HashTable<AnyType, Key, AnyHash>& other);

        //- Generalized means to filter table entries based on their keys.
        //  Keep (or optionally prune) entries with keys that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k);
        //  \endcode
        //
        //  For example,
        //  \code
        //  wordRes goodFields = ...;
        //  allFieldNames.filterKeys
        //  (
        //      [&goodFields](const word& k){ return goodFields.match(k); }
        //  );
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterKeys
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their values.
        //  Keep (or optionally prune) entries with values that satisfy
        //  the unary predicate, which has the following signature:
        //  \code
        //  bool operator()(const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterValues
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Generalized means to filter table entries based on their key/value.
        //  Keep (or optionally prune) entries with keys/values that satisfy
        //  the binary predicate, which has the following signature:
        //  \code
        //  bool operator()(const Key& k, const T& v);
        //  \endcode
        //
        //  \return The number of items changed (removed)
        template<class BinaryPredicate>
        label filterEntries
        (
            const BinaryPredicate& pred,
            const bool pruning = false
        );


        //- Remove all entries from table
        void clear();

        //- Remove all entries from table and the table itself.
        //  Equivalent to clear() followed by setCapacity(0)
        void clearStorage();

        //- Change the hash table capacity (number of buckets).
        //  Setting a zero capacity will only remove the underlying
        //  storage if the table is empty.
        void setCapacity(label newCapacity);

        //- Rehash the hash table with new number of buckets.
        //- Currently identical to setCapacity()
        void resize(label newCapacity);

        //- Reserve space for at least the specified number of elements
        //- (not the number of buckets) and regenerates the hash table.
        void reserve(label numEntries);

        //- Swap contents into this table
        void swap(HashTable<T, Key, Hash>& rhs) noexcept;

        //- Transfer contents into this table.
        void transfer(HashTable<T, Key, Hash>& rhs);

        //- Attempts to extract entries from source parameter and insert them
        //- into \c this, does not overwrite existing entries.
        //- The source will contains any items that could not be merged.
        void merge(HashTable<T, Key, Hash>& source);

        //- Attempts to extract entries from source parameter and insert them
        //- into \c this, does not overwrite existing entries.
        //- The source will contains any items that could not be merged.
        void merge(HashTable<T, Key, Hash>&& source);


    // Member Operators

        //- Find and return a hashed entry. FatalError if it does not exist.
        //  Same as at().
        inline T& operator[](const Key& key);

        //- Find and return a hashed entry. FatalError if it does not exist.
        //  Same as at().
        inline const T& operator[](const Key& key) const;

        //- Return existing entry or create a new entry.
        //  A newly created entry is created as a nameless T() and is thus
        //  value-initialized. For primitives, this will be zero.
        inline T& operator()(const Key& key);

        //- Return existing entry or insert a new entry.
        inline T& operator()(const Key& key, const T& deflt);

        //- Copy assign
        void operator=(const this_type& rhs);

        //- Copy assign from an initializer list
        //  Duplicate entries are handled by overwriting
        void operator=(std::initializer_list<std::pair<Key, T>> rhs);

        //- Move assign
        void operator=(this_type&& rhs);

        //- Equality. Tables are equal if all keys and values are equal,
        //- independent of order or underlying storage size.
        bool operator==(const this_type& rhs) const;

        //- The opposite of the equality operation.
        bool operator!=(const this_type& rhs) const;

        //- Add entries into this HashTable
        this_type& operator+=(const this_type& rhs);


protected:

    // Iterators and helpers

        //- Internally used base for iterator and const_iterator
        template<bool Const> class Iterator;

        //- Allow iterator access to HashTable internals.
        friend class Iterator<true>;

        //- Allow iterator access to HashTable internals.
        friend class Iterator<false>;

        //- The iterator base for HashTable (internal use only).
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        //  iterator and const_iterator have the same size, allowing
        //  us to reinterpret_cast between them (if desired)

        template<bool Const>
        class Iterator
        {
        public:

        // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type = this_type::difference_type;

            //- The HashTable container type
            using table_type = typename std::conditional
            <
                Const,
                const this_type,
                this_type
            >::type;

            //- The node-type being addressed
            using node_type = typename std::conditional
            <
                Const,
                const this_type::node_type,
                this_type::node_type
            >::type;

            //- The key type
            using key_type = this_type::key_type;

            //- The object type being addressed
            using mapped_type = typename std::conditional
            <
                Const,
                const this_type::mapped_type,
                this_type::mapped_type
            >::type;


        protected:

        // Protected Data

            //- The selected entry.
            //  MUST be the first member for easy comparison between iterators
            //  and to support reinterpret_cast from nullObject
            node_type* entry_;

            //- The hash-table container being iterated on.
            //  Uses pointer for default copy/assignment
            table_type* container_;

            //- Index within the hash-table data.
            //  A signed value, since iterator_erase() needs a negative value
            //  to mark the position.
            label index_;

            // Friendship with HashTable, for begin/find constructors
            friend class HashTable;


        // Protected Constructors

            //- Default construct. Also the same as the end iterator
            inline constexpr Iterator() noexcept;

            //- Construct from begin of hash-table
            inline explicit Iterator(table_type* tbl);

            //- Construct by finding key in hash table
            Iterator(table_type* tbl, const Key& key);


        // Protected Member Functions

            //- Increment to the next position
            inline void increment();

            //- Permit explicit cast to the other (const/non-const) iterator
            template<bool Any>
            explicit operator const Iterator<Any>&() const
            {
                return *reinterpret_cast<const Iterator<Any>*>(this);
            }


        public:

        // Member Functions

            //- True if iterator points to an entry
            //  This can be used directly instead of comparing to end()
            bool good() const noexcept { return entry_; }

            //- True if iterator points to an entry - same as good()
            bool found() const noexcept { return entry_; }

            //- The key associated with the iterator
            const Key& key() const { return entry_->key(); }

            //- Write the (key, val) pair
            inline Ostream& print(Ostream& os) const;


        // Member Operators

            //- True if iterator points to an entry
            //  This can be used directly instead of comparing to end()
            explicit operator bool() const noexcept { return entry_; }

            //- Compare hash-entry element pointers.
            //  Independent of const/non-const access
            template<bool Any>
            bool operator==(const Iterator<Any>& iter) const noexcept
            {
                return (entry_ == iter.entry_);
            }

            template<bool Any>
            bool operator!=(const Iterator<Any>& iter) const noexcept
            {
                return (entry_ != iter.entry_);
            }
        };

public:

        //- Forward iterator with non-const access
        class iterator
        :
            public Iterator<false>
        {
        public:

        // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = this_type::difference_type;

            using node_type   = this_type::node_type;
            using key_type    = this_type::key_type;
            using mapped_type = this_type::mapped_type;
            using value_type  = this_type::value_type;
            using pointer     = this_type::pointer;
            using reference   = this_type::reference;
            using const_pointer = this_type::const_pointer;
            using const_reference = this_type::const_reference;


        // Constructors

            //- Default construct (end iterator)
            iterator() = default;

            //- Copy construct from similar access type
            explicit iterator(const Iterator<false>& iter)
            :
                Iterator<false>(iter)
            {}


        // Member Functions/Operators

            //- Const access to the entry (node)
            const node_type* node() const noexcept
            {
                return Iterator<false>::entry_;
            }

            //- Non-const access to the entry (node)
            node_type* node() noexcept
            {
                return Iterator<false>::entry_;
            }

            //- Const access to referenced object (value)
            const_reference val() const
            {
                return Iterator<false>::entry_->cval();
            }

            //- Non-const access to referenced object (value)
            reference val()
            {
                return Iterator<false>::entry_->val();
            }

            //- Const access to referenced object (value)
            const_reference operator*() const  { return this->val(); }
            const_reference operator()() const { return this->val(); }

            //- Non-const access to referenced object (value)
            reference operator*()  { return this->val(); }
            reference operator()() { return this->val(); }

            inline iterator& operator++();
            inline iterator operator++(int);
        };


    // STL const_iterator

        //- Forward iterator with const access
        class const_iterator
        :
            public Iterator<true>
        {
        public:

        // Typedefs
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = this_type::difference_type;

            using node_type   = this_type::node_type;
            using key_type    = this_type::key_type;
            using mapped_type = const this_type::mapped_type;
            using value_type  = const this_type::value_type;
            using pointer     = this_type::const_pointer;
            using reference   = this_type::const_reference;


        // Generated Methods

            //- Default construct (end iterator)
            const_iterator() = default;

            //- Copy construct
            const_iterator(const const_iterator&) = default;

            //- Copy assignment
            const_iterator& operator=(const const_iterator&) = default;


        // Constructors

            //- Copy construct from any access type
            template<bool Any>
            const_iterator(const Iterator<Any>& iter)
            :
                Iterator<true>(static_cast<const Iterator<Any>&>(iter))
            {}

            //- Implicit conversion from dissimilar access type
            const_iterator(const iterator& iter)
            :
                const_iterator(reinterpret_cast<const const_iterator&>(iter))
            {}


        // Member Functions/Operators

            //- Const access to the entry (node)
            const node_type* node() const noexcept
            {
                return Iterator<true>::entry_;
            }

            //- Const access to referenced object (value)
            reference val() const
            {
                return Iterator<true>::entry_->cval();
            }

            //- Const access to referenced object (value)
            reference operator*() const  { return this->val(); }
            reference operator()() const { return this->val(); }

            inline const_iterator& operator++();
            inline const_iterator operator++(int);


        // Assignment

            // Allow assign from iterator to const_iterator
            const_iterator& operator=(const iterator& iter)
            {
                return this->operator=
                (
                    reinterpret_cast<const const_iterator&>(iter)
                );
            }
        };


    // Iterator (keys)

        //- An iterator wrapper for returning a reference to the key
        template<class Iter>
        class key_iterator_base
        :
            public Iter
        {
        public:

            using value_type = this_type::key_type;
            using pointer    = const Key*;
            using reference  = const Key&;

            //- Default construct (end iterator)
            constexpr key_iterator_base() noexcept
            :
                Iter()
            {}

            //- Copy construct with implicit conversion
            explicit key_iterator_base(const Iter& iter)
            :
                Iter(iter)
            {}

            //- Return the key
            reference operator*() const  { return this->key(); }
            reference operator()() const { return this->key(); }

            key_iterator_base& operator++()
            {
                this->increment();
                return *this;
            }

            key_iterator_base operator++(int)
            {
                key_iterator_base iter(*this);
                this->increment();
                return iter;
            }
        };


        //- Forward iterator returning the key
        using key_iterator = key_iterator_base<iterator>;

        //- Forward const iterator returning the key
        using const_key_iterator = key_iterator_base<const_iterator>;

        //- A const iterator begin/end pair for iterating over keys
        const_iterator_pair<const_key_iterator, this_type> keys() const
        {
            return const_iterator_pair<const_key_iterator, this_type>(*this);
        }


    // Iterator access

        //- iterator set to the beginning of the HashTable
        inline iterator begin();

        //- const_iterator set to the beginning of the HashTable
        inline const_iterator begin() const;

        //- const_iterator set to the beginning of the HashTable
        inline const_iterator cbegin() const;

        //- iterator to signal the end (for any HashTable)
        inline iterator end() noexcept;

        //- const_iterator to signal the end (for any HashTable)
        inline const_iterator end() const noexcept;

        //- const_iterator to signal the end (for any HashTable)
        inline constexpr const_iterator cend() const noexcept;


    // Reading/writing

        //- Print information
        Ostream& printInfo(Ostream& os) const;

        //- Write unordered keys (list), with line-breaks
        //- when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeKeys(Ostream& os, const label shortLen=0) const;


    // IOstream Operators

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            HashTable<T, Key, Hash>& tbl
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const HashTable<T, Key, Hash>& tbl
        );


    // Housekeeping

        //- Same as contains()
        bool found(const Key& key) const { return this->contains(key); }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const node_type> sorted() const { return this->csorted(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashTableI.H"
#include "HashTableIterI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTableC  // Excluded from token.H
#ifdef NoRepository
    #include "HashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
