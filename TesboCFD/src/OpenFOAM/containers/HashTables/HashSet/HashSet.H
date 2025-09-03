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

Class
    Foam::HashSet

Description
    A HashTable with keys but without contents that is similar to
    \c std::unordered_set.

    The entries are considered \a unordered since their placement
    depends on the method used to generate the hash key index, the
    table capacity, insertion order etc. When the key order is
    important, use the sortedToc() method to obtain a list of sorted
    keys and use that for further access.

Note
    The HashSet iterator dereferences to the key, so the following
    range-for works as expected:
    \code
        labelHashSet someLabels{10, 20, 30, 40, ...};
        for (const label i : someLabels)
        {
            Info<< "val:" << i << nl;
        }
    \endcode

Typedef
    Foam::wordHashSet

Description
    A HashSet with word keys and string hasher.

Typedef
    Foam::labelHashSet

Description
    A HashSet with label keys and label hasher.

\*---------------------------------------------------------------------------*/

#ifndef Foam_HashSet_H
#define Foam_HashSet_H

#include "HashTable.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class MinMax;
template<class Key, class Hash> class HashSet;

// Common hash-set types

//- A HashSet of words, uses string hasher.
typedef HashSet<word, Hash<word>> wordHashSet;

//- A HashSet of labels, uses label hasher.
typedef HashSet<label, Hash<label>> labelHashSet;


/*---------------------------------------------------------------------------*\
                           Class HashSet Declaration
\*---------------------------------------------------------------------------*/

template<class Key, class Hash=Foam::Hash<Key>>
class HashSet
:
    public HashTable<Foam::zero, Key, Hash>
{
    // Private Member Functions

        //- Assign from the input iterator range
        template<class InputIter>
        inline label assignMany
        (
            const label nItems,
            InputIter first,
            InputIter last
        );


public:

    //- The template instance used for this HashSet
    typedef HashSet<Key, Hash> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<Foam::zero, Key, Hash> parent_type;

    //- An iterator, returning reference to the key
    using iterator = typename parent_type::key_iterator;

    //- A const_iterator, returning reference to the key
    using const_iterator = typename parent_type::const_key_iterator;


    // Constructors

        //- Default construct: empty without allocation (capacity=0)
        HashSet() noexcept = default;

        //- Construct empty without allocation (capacity=0)
        explicit HashSet(const Foam::zero) noexcept
        :
            parent_type(Foam::zero{})
        {}

        //- Construct empty with initial table capacity
        explicit HashSet(const label initialCapacity)
        :
            parent_type(initialCapacity)
        {}

        //- Copy construct
        HashSet(const this_type& rhs)
        :
            parent_type(rhs)
        {}

        //- Move construct
        HashSet(this_type&& rhs) noexcept
        :
            parent_type(std::move(rhs))
        {}

        //- Construct from Istream with default initial table capacity
        explicit HashSet(Istream& is)
        :
            parent_type(is)
        {}

        //- Construct from FixedList of Key
        template<unsigned N>
        explicit HashSet(const FixedList<Key, N>& list);

        //- Construct from UList of Key
        explicit HashSet(const UList<Key>& list);

        //- Construct from an indirect list
        template<class Addr>
        explicit HashSet(const IndirectListBase<Key, Addr>& list);

        //- Construct from an initializer list of Key
        HashSet(std::initializer_list<Key> list);

        //- Construct from the keys of another HashTable,
        //- the type of values held is arbitrary.
        template<class AnyType, class AnyHash>
        explicit HashSet(const HashTable<AnyType, Key, AnyHash>& tbl);


    // Member Functions

        //- Same as contains() - return true if key exists in the set.
        //  Method name compatibility with bitSet and boolList.
        bool test(const Key& key) const
        {
            return this->contains(key);
        }


    // Edit

        //- Insert a new entry, not overwriting existing entries.
        //  \return True if the entry inserted, which means that it did
        //  not previously exist in the set.
        bool insert(const Key& key)
        {
            return this->parent_type::emplace(key);
        }

        //- Same as insert (no value to overwrite)
        bool set(const Key& key)
        {
            return this->parent_type::emplace(key);
        }

        //- Unset the specified key - same as erase
        //  \return True if the entry existed and was removed
        bool unset(const Key& key)
        {
            return this->parent_type::erase(key);
        }

        //- Attempts to extract entries from source parameter and insert them
        //- into \c this, does not overwrite existing entries.
        //- The source will contains any items that could not be merged.
        void merge(HashSet<Key, Hash>& source);

        //- Attempts to extract entries from source parameter and insert them
        //- into \c this, does not overwrite existing entries.
        //- The source will contains any items that could not be merged.
        void merge(HashSet<Key, Hash>&& source);


    // Convenience

        //- Insert keys from the input iterator range
        //  \return The number of new elements inserted
        template<class InputIter>
        inline label insert(InputIter first, InputIter last);

        //- Insert keys from a initializer list of Key
        //  \return The number of new elements inserted
        inline label insert(std::initializer_list<Key> list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        template<unsigned N>
        inline label insert(const FixedList<Key, N>& list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        inline label insert(const UList<Key>& list);

        //- Insert keys from the list of Key
        //  \return The number of new elements inserted
        template<class Addr>
        inline label insert(const IndirectListBase<Key, Addr>& list);

        //- Same as insert (no value to overwrite)
        template<class InputIter>
        inline label set(InputIter first, InputIter last)
        {
            return insert(first, last);
        }

        //- Same as insert (no value to overwrite)
        inline label set(std::initializer_list<Key> list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        template<unsigned N>
        inline label set(const FixedList<Key, N>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        inline label set(const UList<Key>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        template<class Addr>
        inline label set(const IndirectListBase<Key, Addr>& list)
        {
            return insert(list);
        }

        //- Same as insert (no value to overwrite)
        //  \note Method name compatibility with bitSet
        template<class InputIter>
        inline label setMany(InputIter first, InputIter last)
        {
            return insert(first, last);
        }

        //- Unset the keys listed in the input iterator range
        //  \return The number of items removed
        template<class InputIter>
        inline label unset(InputIter first, InputIter last);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        inline label unset(std::initializer_list<Key> list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        template<unsigned N>
        inline label unset(const FixedList<Key, N>& list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        inline label unset(const UList<Key>& list);

        //- Unset the listed keys - same as erase
        //  \return The number of items removed
        template<class Addr>
        inline label unset(const IndirectListBase<Key, Addr>& list);


    // STL iterators

        inline iterator begin();
        inline const_iterator begin() const;
        inline const_iterator cbegin() const;

        inline iterator end() noexcept;
        inline const_iterator end() const noexcept;
        inline constexpr const_iterator cend() const noexcept;


    // Writing

        //- Write unordered keys (list), with line-breaks
        //- when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const
        {
            return this->writeKeys(os, shortLen);
        }


    // Member Operators

        //- Return true if the entry exists, same as contains()
        //  \note this allows use of HashSet as a predicate test
        inline bool operator()(const Key& key) const noexcept;

        //- Return true if the entry exists, same as contains().
        inline bool operator[](const Key& key) const noexcept;

        //- Copy assign
        void operator=(const this_type& rhs)
        {
            parent_type::operator=(rhs);
        }

        //- Move assign
        void operator=(this_type&& rhs)
        {
            parent_type::operator=(std::move(rhs));
        }


    // Comparison

        //- Sets are equal if all keys are equal,
        //- independent of order or underlying storage size.
        bool operator==(const this_type& rhs) const;

        //- The opposite of the equality operation.
        bool operator!=(const this_type& rhs) const;


    // Assignment

        //- Assignment from a UList of keys
        void operator=(const UList<Key>& rhs);

        //- Assignment from a FixedList of keys
        template<unsigned N>
        void operator=(const FixedList<Key, N>& rhs);

        //- Assignment from an initializer list of keys
        void operator=(std::initializer_list<Key> rhs);


    // Logical and set operations

        //- Add entries to this HashSet
        this_type& operator|=(const this_type& rhs);

        //- Only retain entries contained in both HashSets
        inline this_type& operator&=(const this_type& rhs);

        //- Only retain unique entries (xor)
        this_type& operator^=(const this_type& rhs);

        //- Add entries to this HashSet. Same as the '|=' operator
        inline this_type& operator+=(const this_type& rhs);

        //- Remove entries from this HashSet. Uses erase()
        inline this_type& operator-=(const this_type& rhs);


    // Housekeeping

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        List<Key> tocValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        List<Key> tocEntries(const BinaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        label countValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        label countEntries(const BinaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class UnaryPredicate>
        label filterValues(const UnaryPredicate&, const bool) = delete;

        //- Not applicable for HashSet
        template<class BinaryPredicate>
        label filterEntries(const BinaryPredicate&, const bool) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Functions

//- Find the min value in labelHashSet, optionally limited by second argument.
//  For an empty set, returns the second argument (eg, labelMax).
label min(const labelHashSet& set, label minValue = labelMax);

//- Find the max value in labelHashSet, optionally limited by second argument.
//  For an empty set, returns the second argument (eg, labelMin).
label max(const labelHashSet& set, label maxValue = labelMin);

//- Find the min/max values of labelHashSet
MinMax<label> minMax(const labelHashSet& set);


// Global Operators

//- Write the list of HashSet keys
template<class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashSet<Key, Hash>& rhs);


//- Combine entries (OR) for two HashSets
//  See HashSet::operator|= for more details.
template<class Key, class Hash>
HashSet<Key, Hash> operator|
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
);

//- Subset (AND) intersection of two HashSet
//  See HashSet::operator&= for more details.
template<class Key, class Hash>
HashSet<Key, Hash> operator&
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
);

//- Create a HashSet that only contains unique entries (XOR)
//  See HashSet::operator^= for more details.
template<class Key, class Hash>
HashSet<Key, Hash> operator^
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
);

//- Subset difference of two HashSets
//  See HashSet::operator-= for more details.
template<class Key, class Hash>
HashSet<Key, Hash> operator-
(
    const HashSet<Key, Hash>& a,
    const HashSet<Key, Hash>& b
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashSet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
