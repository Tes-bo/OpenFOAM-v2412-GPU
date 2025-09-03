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
    Foam::List

Description
    A 1D array of objects of type \<T\>, where the size of the vector
    is known and used for subscript bounds checking, etc.

    Storage is allocated on free-store during construction.

SourceFiles
    List.C
    ListI.H
    ListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_List_H
#define Foam_List_H

#include "autoPtr.H"
#include "UList.H"
#include "SLListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class List;
template<class T, unsigned N> class FixedList;
template<class T, int SizeMin> class DynamicList;

template<class T> class PtrList;

template<class T> Istream& operator>>(Istream& is, List<T>& list);

// Common list types
typedef List<bool> boolList;        //!< A List of bools
typedef List<char> charList;        //!< A List of chars
typedef List<label> labelList;      //!< A List of labels


/*---------------------------------------------------------------------------*\
                           Class List Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class List
:
    public UList<T>
{
    // Private Member Functions

        //- Allocate list storage
        inline void doAlloc();

        //- Reallocate list storage to the given size
        //  Discards old storage (if any). Does not copy old contents
        inline void reAlloc(const label len);

        //- Copy all list contents. Uses operator[] on the input list
        template<class ListType>
        inline void copyList(const ListType& list);

        //- Copy list contents via indirect indices.
        //  Uses operator[] on the values and indices lists
        template<class ListType, class ListIndices>
        inline void copyList(const ListType& list, const ListIndices& indices);

        //- Change allocation size of List, retaining old contents.
        //  Backend for resize
        void doResize(const label len);

        //- Construct given begin/end iterators and number of elements
        //  Since the size is provided, the end iterator is actually ignored.
        template<class InputIterator>
        inline List
        (
            InputIterator firstIter,
            InputIterator lastIter,  // (unused)
            const label len
        );

        //- Read List from Istream between '(' and ')' delimiters.
        //- The size is not known a priori.
        bool readBracketList(Istream& is);


    // Methods as per DynamicList to simplify code maintenance

        //- Stub method for internal naming as per DynamicList
        void setCapacity_nocopy(const label len) { resize_nocopy(len); }


public:

    // Related types

        //- Declare type of subList
        typedef SubList<T> subList;


    // Static Member Functions

        //- Return a null List (reference to a nullObject).
        //- Behaves like an empty List.
        static const List<T>& null() noexcept
        {
            return NullObjectRef<List<T>>();
        }


    // Constructors

        //- Default construct
        inline constexpr List() noexcept;

        //- Construct with given size
        explicit List(const label len);

        //- Construct with given size and value for all elements
        List(const label len, const T& val);

        //- Construct with given size initializing all elements to zero
        List(const label len, const Foam::zero);

        //- Construct with length=1, copying the value as the only content
        List(const Foam::one, const T& val);

        //- Construct with length=1, moving the value as the only content
        List(const Foam::one, T&& val);

        //- Construct with length=1, initializing content to zero
        List(const Foam::one, const Foam::zero);

        //- Copy construct from list
        List(const List<T>& list);

        //- Copy construct contents from list
        explicit List(const UList<T>& list);

        //- Construct as copy or re-use as specified
        List(List<T>& list, bool reuse);

        //- Copy construct subset of list
        List(const UList<T>& list, const labelUList& indices);

        //- Copy construct subset of list
        template<unsigned N>
        List(const UList<T>& list, const FixedList<label, N>& indices);

        //- Construct as copy of FixedList\<T, N\>
        template<unsigned N>
        explicit List(const FixedList<T, N>& list);

        //- Construct as copy of PtrList<T>
        explicit List(const PtrList<T>& list);

        //- Construct as copy of IndirectList contents
        template<class Addr>
        explicit List(const IndirectListBase<T, Addr>& list);

        //- Construct from an initializer list
        List(std::initializer_list<T> list);

        //- Move construct from List
        List(List<T>&& list) noexcept;

        //- Move construct from DynamicList
        template<int SizeMin>
        List(DynamicList<T, SizeMin>&& list);

        //- Construct from Istream
        List(Istream& is);

        //- Clone
        inline autoPtr<List<T>> clone() const;


    //- Destructor
    ~List();


    // Member Functions

    // Sizing

        //- Clear the list, i.e. set size to zero
        inline void clear();

        //- Adjust allocated size of list.
        //  The boolList version fills new memory with false.
        inline void resize(const label len);

        //- Adjust allocated size of list and set val for \em new elements
        void resize(const label len, const T& val);

        //- Adjust allocated size of list and set val for \em all elements
        inline void resize_fill(const label len, const T& val);

        //- Adjust allocated size of list \b without necessarily
        //  retaining old content.
        //  If no reallocation is required, the contents remain untouched.
        //  Otherwise the contents will be uninitialized.
        inline void resize_nocopy(const label len);

        //- Change the addressed list size directly without affecting
        //- any memory management (advanced usage).
        //
        //  It is left to the caller to avoid \em unsafe lengthening beyond
        //  the allocated memory region.
        inline void resize_unsafe(const label len) noexcept;

        //- Alias for resize()
        void setSize(const label n) { this->resize(n); }

        //- Alias for resize()
        void setSize(const label n, const T& val) { this->resize(n, val); }


    // Edit

        //- Transfer the contents of the argument List into this list
        //- and annul the argument list
        void transfer(List<T>& list);

        //- Transfer the contents of the argument List into this list
        //- and annul the argument list
        template<int SizeMin>
        void transfer(DynamicList<T, SizeMin>& list);

        //- Return subscript-checked element of UList and resizing the list
        //- if required.
        inline T& newElmt(const label i);


    // Edit

        //- Construct an element at the end of the list,
        //- return reference to the new list element.
        //  If this is frequently required, consider a DynamicList instead.
        template<class... Args>
        inline T& emplace_back(Args&&... args);

        //- Append an element at the end of the list
        //  If this is frequently required, consider a DynamicList
        inline void push_back(const T& val);

        //- Move append an element at the end of the list
        //  If this is frequently required, consider a DynamicList
        inline void push_back(T&& val);

        //- Append a List to the end of this list
        //  If this is frequently required, consider a DynamicList
        inline void push_back(const UList<T>& list);

        //- Append IndirectList contents at the end of this list
        //  If this is frequently required, consider a DynamicList
        template<class Addr>
        inline void push_back(const IndirectListBase<T, Addr>& list);

        //- Append an element if not already in the list.
        //  \return the change in list length
        inline label push_uniq(const T& val);

        //- Reduce size by 1 or more elements. Can be called on an empty list.
        inline void pop_back(label n = 1);


    // Member Operators

        //- Assignment to UList operator. Takes linear time
        void operator=(const UList<T>& list);

        //- Assignment operator. Takes linear time
        void operator=(const List<T>& list);

        //- Assignment from IndirectList. Takes linear time
        template<class Addr>
        void operator=(const IndirectListBase<T, Addr>& list);

        //- Copy assignment from FixedList
        template<unsigned N>
        void operator=(const FixedList<T, N>& list);

        //- Assignment to an initializer list
        void operator=(std::initializer_list<T> list);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const Foam::zero);

        //- Move assignment. Takes constant time
        void operator=(List<T>&& list);

        //- Move assignment. Takes constant time.
        template<int SizeMin>
        void operator=(DynamicList<T, SizeMin>&& list);


    // Reading/writing

        //- Read List from Istream, discarding contents of existing List
        Istream& readList(Istream& is);


    // IOstream Operators

        //- Use the readList() method to read contents from Istream.
        friend Istream& operator>> <T>
        (
            Istream& is,
            List<T>& list
        );


    // Housekeeping

        //- No shallowCopy permitted
        void shallowCopy(const UList<T>&) = delete;


    // Special Methods

        //- A bitSet::set() method for a list of bool
        //  Increases size when setting an out-of-bounds value.
        //
        //  \return True if value changed.
        template<class TypeT = T>
        typename std::enable_if<std::is_same<bool, TypeT>::value, bool>::type
        inline set(const label i, bool val = true)
        {
            if (i < 0)
            {
                return false;               // Out-of-bounds: ignore
            }
            else if (i >= this->size())
            {
                if (!val)                   // Unset out-of-bounds: ignore
                {
                    return false;
                }
                this->resize(i+1, false);   // Adjust size for assign, fill 0
            }

            (*this)[i] = val;
            return true;
        }


    // Housekeeping

        //- Append an element at the end of the list
        //  If this is frequently required, consider a DynamicList
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const T& val) { this->push_back(val); }

        //- Move append an element at the end of the list
        //  If this is frequently required, consider a DynamicList
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T&& val) { this->push_back(std::move(val)); }

        //- Append a List to the end of this list
        //  If this is frequently required, consider a DynamicList
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const UList<T>& list) { this->push_back(list); }

        //- Append IndirectList contents at the end of this list
        //  If this is frequently required, consider a DynamicList
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        template<class Addr>
        void append(const IndirectListBase<T, Addr>& list)
        {
            this->push_back(list);
        }

        //- Same as push_uniq()
        FOAM_DEPRECATED_STRICT(2022-10, "push_uniq()")
        label appendUniq(const T& val) { return this->push_uniq(val); }


        //- Copy construct from SLList
        explicit List(const SLList<T>& list);

        //- Copy assign from SLList in linear time
        void operator=(const SLList<T>& list);
};


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

//- Hashing for List data
template<class T>
struct Hash<List<T>> : List<T>::hasher {};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read List contents from Istream
template<class T>
Istream& operator>>(Istream& is, List<T>& list)
{
    return list.readList(is);
}


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return an identity map of the given length with (map[i] == i),
//- works like std::iota() but returning a list of label values.
//  Optionally with an alternative start index, so that (map[i] == i+start)
labelList identity(const label len, label start=0);

//- Return the (stable) sort order for the list
template<class T>
labelList sortedOrder(const UList<T>& input);

//- Generate the (stable) sort order for the list
template<class T>
void sortedOrder(const UList<T>& input, labelList& order);

//- Sort using specified list compare predicate
template<class T, class ListComparePredicate>
void sortedOrder
(
    const UList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "List.C"
    #include "ListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
