/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::UPtrList

Description
    A list of pointers to objects of type \<T\>, without allocation/deallocation
    management of the pointers - this is to be done elsewhere.
    The operator[] returns a reference to the object (not the pointer).

    The iterators are similar to bitSet in that they skip nullptr entries,
    and also return a value (like the list operator[] does).

    When traversing lists, it possible to test the validity directly:
    \code
    forAll(interfaces, i)
    {
        if (interfaces.test(i))
        {
            // Interface is set, do something
            const auto& intf = interfaces[i];
            ...
        }
    }
    \endcode
    The lists can also be traversed with a for-range
    (in OpenFOAM-v2212 and earlier this would have failed on nullptr entries):
    \code
    for (const auto& intf : interfaces)
    {
        // Do something
        ...
    }
    \endcode

    It is also possible to traverse with non-null entries
    and use key/val access (like HashTable naming):
    \code
    forAllConstIters(interfaces, iter)
    {
        Info<< "entry " << iter.key() << " : " << iter.val() << nl;
    }
    \endcode

Note
    The class definition is such that it contains a list of pointers, but
    itself does not inherit from a list of pointers since this would
    wreak havoc later with inheritance resolution.

See Also
    Foam::PtrList
    Foam::PtrDynList

SourceFiles
    UPtrListI.H
    UPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPtrList_H
#define Foam_UPtrList_H

#include "PtrListDetail.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class PtrList;
template<class T> class UPtrList;
template<class T> Ostream& operator<<(Ostream& os, const UPtrList<T>& list);

/*---------------------------------------------------------------------------*\
                          Class UPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrList
{
protected:

    // Protected Member Data

        //- The list of pointers
        Detail::PtrListDetail<T> ptrs_;


    // Protected Member Functions

        //- Adjust addressable size
        inline void setAddressableSize(const label n) noexcept;

        //- The next non-null entry after the specified position
        inline label find_next(label pos) const;


    // Constructors

        //- Low-level move construct
        inline explicit UPtrList(Detail::PtrListDetail<T>&& ptrs) noexcept;


public:

    // STL type definitions

        //- Type of values the list contains
        typedef T value_type;

        //- A non-const reference to the value_type
        typedef T& reference;

        //- A const reference to the value_type
        typedef const T& const_reference;

        //- Forward iterator with non-const access
        class iterator;

        //- Forward iterator with const access
        class const_iterator;


    // Public Classes

        //- A wrapper for a binary comparison of values that interjects
        //- pointer dereferencing with null pointer guards.
        //  It will also sort any null pointers to the end
        //  (eg, rubbish that can be truncated)
        template<class Compare>
        struct value_compare
        {
            const Compare& comp;

            value_compare(const Compare& cmp)
            :
                comp(cmp)
            {}

            //- Compare dereferenced pointers
            bool operator()(const T* const a, const T* const b) const
            {
                return (a && b) ? comp(*a, *b) : !b;
            }
        };

        //- A UPtrList compare binary predicate for normal sort order.
        //- Null entries (if any) sort to the end.
        struct less
        {
            const UPtrList<T>& values;

            less(const UPtrList<T>& list)
            :
                values(list)
            {}

            //- Compare dereferenced pointer locations for normal sort.
            bool operator()(const label ai, const label bi) const
            {
                const T* const a = values.get(ai);
                const T* const b = values.get(bi);

                return (a && b) ? (*a < *b) : !b;
            }
        };

        //- A UPtrList compare binary predicate for reverse sort order.
        //  Null entries (if any) sort to the end.
        struct greater
        {
            const UPtrList<T>& values;

            greater(const UPtrList<T>& list)
            :
                values(list)
            {}

            //- Compare dereferenced pointer locations for reverse sort
            bool operator()(const label ai, const label bi) const
            {
                const T* const a = values.get(ai);
                const T* const b = values.get(bi);

                return (a && b) ? (*b < *a) : !a;
            }
        };


    // Constructors

        //- Default construct
        constexpr UPtrList() noexcept = default;

        //- Construct with specified size and set all entries to \c nullptr
        inline explicit UPtrList(const label len);

        //- Copy construct (shallow copies addresses)
        inline UPtrList(const UPtrList<T>& list);

        //- Move construct
        inline UPtrList(UPtrList<T>&& list) noexcept;

        //- Construct as shallow copy or re-use as specified
        inline UPtrList(UPtrList<T>& list, bool reuse);

        //- Shallow copy from PtrList.
        //  The argument is non-const to reflect that the UPtrList can change
        //  the values (not the addresses) within the original list.
        explicit UPtrList(PtrList<T>& list);

        //- Construct from UList of pointers (shallow copy)
        inline explicit UPtrList(const UList<T*>& list);

        //- Construct from UList, taking the address of each list element
        //  The argument is non-const to reflect that the UPtrList can change
        //  the values of the original list.
        inline explicit UPtrList(UList<T>& list);


    // Member Functions

    // Access

        //- True if the list is empty (ie, size() is zero)
        inline bool empty() const noexcept;

        //- The number of entries in the list
        inline label size() const noexcept;

        //- Size of the underlying storage.
        inline label capacity() const noexcept;

        //- The number of non-nullptr entries in the list
        inline label count_nonnull() const noexcept;

        //- Reference to the first element of the list
        inline T& front();

        //- Reference to first element of the list
        inline const T& front() const;

        //- Reference to the last element of the list
        inline T& back();

        //- Reference to the last element of the list
        inline const T& back() const;

        //- Return const pointer to element (can be nullptr),
        //- or nullptr for out-of-range access (ie, \em with bounds checking).
        //  The return value can be tested as a bool.
        inline const T* test(const label i) const;

        //- Return const pointer to element (can be nullptr),
        //- or nullptr for out-of-range access (ie, \em with bounds checking).
        //  The return value can be tested as a bool.
        inline const T* get(const label i) const;

        //- Return pointer to element (can be nullptr),
        //- or nullptr for out-of-range access (ie, \em with bounds checking).
        //  The return value can be tested as a bool.
        inline T* get(const label i);

        //- Return const pointer to element (can be nullptr),
        //- or nullptr for out-of-range access (ie, \em with bounds checking).
        //  The return value can be tested as a bool.
        const T* set(const label i) const { return this->get(i); }


    // Edit

        //- Set list size to zero.
        inline void clear();

        //- Nullify all entries. Does not change the list size.
        inline void free();

        //- Change the size of the list.
        //- Any new entries are \c nullptr.
        inline void resize(const label newLen);

        //- Set the list to the given size
        //- and set \em all entries to \c nullptr.
        inline void resize_null(const label newLen);

        //- Squeeze out nullptr entries in the list of pointers after which
        //- any null pointers will be at the end of the list
        //  \return the number of non-null entries
        label squeezeNull();

        //- Reduce addressable list size to ignore any trailing null pointers.
        //  The reduces the effective list length without reallocation
        void trimTrailingNull();

        //- Append an element to the end of the list
        inline void push_back(T* ptr);

        //- Move append another list to the end of this list.
        inline void push_back(UPtrList<T>&& other);

        //- Swap content
        inline void swap(UPtrList<T>& list) noexcept;

        //- Transfer contents into this list and annul the argument
        inline void transfer(UPtrList<T>& list);

        //- Set element to specified pointer and return the old list element,
        //- which can be a nullptr.
        //  No-op if the new pointer value is identical to the current content.
        inline T* set(const label i, T* ptr);

        //- Reorder elements.
        //- Reordering must be unique (ie, shuffle).
        //  Optionally check that all pointers have been set.
        void reorder(const labelUList& oldToNew, const bool check = false);

        //- Reorder elements according to new order mapping (newToOld).
        //- Reordering must be unique (ie, shuffle).
        //  Optionally check that all pointers have been set.
        void sortOrder(const labelUList& order, const bool check = false);


    // Checks

        //- Check and raise FatalError if any nullptr exists in the list
        inline void checkNonNull() const;

        //- Return const reference to the element at given position.
        //- FatalError for bounds problem or nullptr.
        inline const T& at(const label i) const;

        //- Return reference to the element at given position.
        //- FatalError for bounds problem or nullptr.
        inline T& at(const label i);


    // Member Operators

        //- Return const reference to the element at given position.
        //- FatalError for bounds problem or nullptr. Same as at().
        inline const T& operator[](const label i) const;

        //- Return reference to the element at given position.
        //- FatalError for bounds problem or nullptr. Same as at().
        inline T& operator[](const label i);

        //- Deprecated(2022-09) - same as get()
        //  \deprecated(2022-09) - use get(), set() or test() methods
        FOAM_DEPRECATED_FOR(2022-09, "get(), set() or test() methods")
        const T* operator()(const label i) const { return this->get(i); }

        //- Copy assignment (shallow copies addresses)
        inline void operator=(const UPtrList<T>& list);

        //- Move assignment
        inline void operator=(UPtrList<T>&& list);


    // IOstream Operators

        //- Write UPtrList to Ostream
        friend Ostream& operator<< <T>(Ostream& os, const UPtrList<T>& list);


protected:

    // Iterators and helpers

        //- Internally used base for iterator and const_iterator
        template<bool Const> class Iterator;

        //- Allow iterator access to internals
        friend class Iterator<true>;

        //- Allow iterator access to internals
        friend class Iterator<false>;


        //- The iterator base for UPtrList (internal use only).
        //  Iterates non-nullptr entries.

        template<bool Const>
        class Iterator
        {
        public:

        // Typedefs

            //- The list container type
            using list_type = typename std::conditional
            <
                Const,
                const UPtrList<T>,
                UPtrList<T>
            >::type;


        protected:

        // Protected Data

            //- The parent being iterated
            //  Uses pointer for default copy/assignment
            list_type* list_;

            //- The position within the list
            label pos_;

            // Friendship with UPtrList, for begin constructor
            friend class UPtrList;


        // Protected Constructors

            //- Default construct. Also the same as the end iterator
            inline constexpr Iterator() noexcept;

            //- Construct begin iterator
            inline explicit Iterator(list_type* list);


        // Protected Member Functions

            //- Increment to the next non-null position
            inline void increment();

            //- Permit explicit cast to the other (const/non-const) iterator
            template<bool Any>
            explicit operator const Iterator<Any>&() const
            {
                return *reinterpret_cast<const Iterator<Any>*>(this);
            }


        public:

        // Member Functions/Operators

            //- True if iterator points to a non-null entry
            bool good() const noexcept { return (list_ && pos_ >= 0); }

            //- The iterator position/index within the list
            label key() const noexcept { return pos_; }

            //- Compare hash-entry element pointers.
            //  Independent of const/non-const access
            template<bool Any>
            bool operator==(const Iterator<Any>& iter) const noexcept
            {
                return (pos_ == iter.pos_);
            }

            template<bool Any>
            bool operator!=(const Iterator<Any>& iter) const noexcept
            {
                return (pos_ != iter.pos_);
            }
        };


public:

    // Iterators

        //- Forward iterator with non-const access
        class iterator : public Iterator<false>
        {
        public:
            // using iterator_category = std::forward_iterator_tag;
            // using difference_type = label;
            using pointer = T*;
            using reference = T&;


        // Constructors

            //- Default construct - an end iterator
            constexpr iterator() noexcept = default;

            //- Copy construct from similar access type
            explicit iterator(const Iterator<false>& iter)
            :
                Iterator<false>(iter)
            {}


        // Member Functions/Operators

            //- Pointer to the referenced object (failsafe)
            inline pointer get() const;

            //- Reference to the object
            inline reference val() const;

            //- Pointer to the referenced object
            pointer operator->() const { return this->get(); }

            //- Reference to the object
            reference operator*() const { return this->val(); }

            //- Legacy call operator: reference to the object
            reference operator()() const { return this->val(); }

            //- Move to the next non-nullptr entry
            inline iterator& operator++();
            inline iterator operator++(int);
        };


        //- Forward iterator with const access
        class const_iterator : public Iterator<true>
        {
        public:
            // using iterator_category = std::forward_iterator_tag;
            // using difference_type = label;
            using pointer = const T*;
            using reference = const T&;


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

            //- Pointer to the referenced object (failsafe)
            inline pointer get() const;

            //- Reference to the object
            inline reference val() const;

            //- Pointer to the referenced object
            pointer operator->() const { return this->get(); }

            //- Reference to the object
            reference operator*() const { return this->val(); }

            //- Legacy call operator: reference to the object
            reference operator()() const { return this->val(); }

            //- Move to the next non-nullptr entry
            inline const_iterator& operator++();
            inline const_iterator operator++(int);
        };


        //- Iterator to begin of raw pointers traversal (use with caution)
        T** begin_ptr() noexcept { return ptrs_.begin(); }

        //- Iterator beyond end of raw pointers traversal (use with caution)
        T** end_ptr() noexcept { return ptrs_.end(); }


        //- Return iterator to begin traversal of non-nullptr entries.
        inline iterator begin();

        //- Return iterator beyond end of UPtrList traversal
        inline iterator end() noexcept;

        //- Return const_iterator to begin traversal of non-nullptr entries.
        inline const_iterator cbegin() const;

        //- Return const_iterator beyond end of UPtrList traversal
        inline const_iterator cend() const noexcept;

        //- Return const_iterator to begin traversal of non-nullptr entries.
        inline const_iterator begin() const;

        //- Return const_iterator beyond end of UPtrList traversal
        inline const_iterator end() const noexcept;


    // Housekeeping

        //- Alias for resize()
        void setSize(const label n) { this->resize(n); }

        //- Reference to the first element of the list
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        T& first() { return front(); }

        //- Return reference to first element of the list
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const T& first() const { return front(); }

        //- Return reference to the last element of the list
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        T& last() { return back(); }

        //- Return reference to the last element of the list
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const T& last() const{ return back(); }

        //- Append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T* ptr) { this->push_back(ptr); }

        //- Move append another list to the end of this list.
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(UPtrList<T>&& other) { this->push_back(std::move(other)); }

        //- The number of non-nullptr entries in the list
        FOAM_DEPRECATED_FOR(2024-01, "count_nonnull()")
        label count() const noexcept { return count_nonnull(); }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Inplace (stable) sorting of pointer list.
//  This sort function includes null pointer guards and will also sort
//  any null pointers to the end (eg, rubbish that can be truncated)
template<class T>
void sort(UPtrList<T>& list);

//- Inplace (stable) sorting of pointer list using given comparator,
//- which compares objects, not pointers.
//  This sort function includes null pointer guards and will also sort
//  any null pointers to the end (eg, rubbish that can be truncated)
template<class T, class Compare>
void sort(UPtrList<T>& list, const Compare& comp);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
