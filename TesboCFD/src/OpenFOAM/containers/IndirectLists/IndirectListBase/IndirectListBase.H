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

Class
    Foam::IndirectListBase

Description
    Base for lists with indirect addressing, templated on the list contents
    type and the addressing type. Storage for both values and addressing
    is held outside of the class.

SourceFiles
    IndirectListBase.C
    IndirectListBaseI.H
    IndirectListBaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IndirectListBase_H
#define Foam_IndirectListBase_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class IndirectListBase Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Addr>
class IndirectListBase
{
    // Private Data

        //- The list values
        UList<T>& values_;

        //- Reference to the addressing for the list values
        const Addr& addr_;


protected:

    // Protected Member Functions

        //- Deep copy values from the list
        template<class ListType>
        inline void copyList(const ListType& rhs);


public:

    // Type definitions (STL)

        //- Type of values the list contains.
        typedef T value_type;

        //- The pointer type for non-const access to value_type items
        typedef T* pointer;

        //- The pointer type for const access to value_type items
        typedef const T* const_pointer;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type to represent the size of a UList
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;

        //- Forward iterator with non-const access
        class iterator;

        //- Forward iterator with const access
        class const_iterator;

        //- The addressing type (non-stl definition)
        typedef Addr addressing_type;


    // Constructors

        //- No default construct
        IndirectListBase() = delete;

        //- Store references to the values list and the addressing array
        inline IndirectListBase(const UList<T>& values, const Addr& addr);


    // Member Functions

    // Access

        //- True if the list is empty (ie, size() is zero).
        bool empty() const noexcept { return addr_.empty(); }

        //- The number of elements in the list
        label size() const noexcept { return addr_.size(); }

        //- The list of values (without addressing)
        const UList<T>& values() const noexcept { return values_; }

        //- The list of values (without addressing)
        UList<T>& values() noexcept { return values_; }

        //- The addressing used for the list
        const Addr& addressing() const noexcept { return addr_; }


        //- True if all entries have identical values, and list is non-empty
        inline bool uniform() const;

        //- The first element of the list.
        inline const T& front() const;

        //- The first element of the list.
        inline T& front();

        //- The last element of the list.
        inline const T& back() const;

        //- The last element of the list.
        inline T& back();

        //- The forward circular index. The next index in the list
        //- which returns to the first at the end of the list
        inline label fcIndex(const label i) const;

        //- The reverse circular index. The previous index in the list
        //- which returns to the last at the beginning of the list
        inline label rcIndex(const label i) const;

        //- Return forward circular value (ie, next value in the list)
        inline const T& fcValue(const label i) const;

        //- Return forward circular value (ie, next value in the list)
        inline T& fcValue(const label i);

        //- Return reverse circular value (ie, previous value in the list)
        inline const T& rcValue(const label i) const;

        //- Return reverse circular value (ie, previous value in the list)
        inline T& rcValue(const label i);

        //- Return the addressed elements as a List
        inline List<T> list() const;


    // Search

        //- Is the value contained in the list?
        //  \param val The value to search for
        //  \param pos The first position to examine (default: 0, no-op if -ve)
        //  \param len The length of the search region (-ve until the end)
        //  \return true if found.
        inline bool contains(const T& val, label pos = 0, label len = -1) const;

        //- Find index of the first occurrence of the value.
        //  \param val The value to search for
        //  \param pos The first position to examine (default: 0, no-op if -ve)
        //  \param len The length of the search region (-ve until the end)
        //  \return position in list or -1 if not found.
        label find(const T& val, label pos = 0, label len = -1) const;

        //- Find index of the last occurrence of the value.
        //  Any occurrences after the end pos are ignored.
        //  Linear search.
        //  \return -1 if not found.
        label rfind(const T& val, label pos = -1) const;


    // Member Operators

        //- Return the addressed elements as a List
        List<T> operator()() const { return this->list(); }

        //- Non-const access to an element in the list
        inline T& operator[](const label i);

        //- Const access to an element in the list
        inline const T& operator[](const label i) const;

        //- Assign all addressed elements to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const Foam::zero);

        //- Deep copy values from a list of the addressed elements
        //  Fatal if list sizes are not identical
        inline void operator=(const UList<T>& rhs);

        //- Deep copy values from a list of the addressed elements
        //  Fatal if list sizes are not identical
        inline void operator=(const IndirectListBase<T, Addr>& rhs);

        //- Deep copy values from a list of the addressed elements
        //  Fatal if list sizes are not identical
        template<class AnyAddr>
        inline void operator=(const IndirectListBase<T, AnyAddr>& rhs);


    // Iterators

        //- A non-const iterator for an indirect list
        //  Only supports forward prefix increment, since the addressing
        //  may/may not support postfix or decrement.
        class iterator
        {
            T* begin_;
            typename addressing_type::const_iterator iter_;

        public:

            using difference_type = label;
            using value_type = T;
            using pointer = T*;
            using reference = T&;
            using iterator_category = std::forward_iterator_tag;

            iterator
            (
                UList<T>& list,
                typename addressing_type::const_iterator addrIter
            )
            :
                begin_(list.data()),
                iter_(addrIter)
            {}

            reference operator*() const { return *(begin_ + *iter_); }

            iterator& operator++()
            {
                ++iter_;
                return *this;
            }

            bool operator==(const iterator& rhs) const
            {
                return (iter_ == rhs.iter_);
            }

            bool operator!=(const iterator& rhs) const
            {
                return (iter_ != rhs.iter_);
            }
        };


        //- A const iterator for an indirect list
        //  Only supports forward prefix increment, since the addressing
        //  may/may not support postfix or decrement.
        class const_iterator
        {
            const T* begin_;
            typename addressing_type::const_iterator iter_;

        public:

            using difference_type = label;
            using value_type = const T;
            using pointer = const T*;
            using reference = const T&;
            using iterator_category = std::forward_iterator_tag;

            const_iterator
            (
                const UList<T>& list,
                typename addressing_type::const_iterator addrIter
            )
            :
                begin_(list.cdata()),
                iter_(addrIter)
            {}

            reference operator*() const { return *(begin_ + *iter_); }

            const_iterator& operator++()
            {
                ++iter_;
                return *this;
            }

            bool operator==(const const_iterator& rhs) const
            {
                return (iter_ == rhs.iter_);
            }

            bool operator!=(const const_iterator& rhs) const
            {
                return (iter_ != rhs.iter_);
            }
        };


    // Iterator (non-const)

        //- Return an iterator at begin of list
        inline iterator begin()
        {
            return iterator(values_, addr_.cbegin());
        }

        //- Return an iterator at end of list
        inline iterator end()
        {
            return iterator(values_, addr_.cend());
        }


    // Iterator (const)

        //- Return a const_iterator at begin of list
        inline const_iterator cbegin() const
        {
            return const_iterator(values_, addr_.cbegin());
        }

        //- Return a const_iterator at end of list
        inline const_iterator cend() const
        {
            return const_iterator(values_, addr_.cend());
        }

        //- Return a const_iterator at begin of list
        const_iterator begin() const { return cbegin(); }

        //- Return a const_iterator at end of list
        const_iterator end() const { return cend(); }


    // Writing

        //- Write List, with line-breaks in ASCII when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;


    // Housekeeping

        //- Access first element of the list, position [0]
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        T& first() { return front(); }

        //- Access first element of the list
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const T& first() const { return front(); };

        //- Access last element of the list, position [size()-1]
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        T& last() { return back(); }

        //- Access last element of the list, position [size()-1]
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const T& last() const { return back(); };

        //- Same as contains()
        bool found(const T& val, label pos = 0) const
        {
            return this->contains(val, pos);
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Write List to Ostream, as per UList::writeList() with default length.
//  The default short-length is given by Detail::ListPolicy::short_length
template<class T, class Addr>
Ostream& operator<<(Ostream& os, const IndirectListBase<T, Addr>& list)
{
    return list.writeList(os, Detail::ListPolicy::short_length<T>::value);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListBaseI.H"

#ifdef NoRepository
    #include "IndirectListBase.C"
    #include "IndirectListBaseIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
