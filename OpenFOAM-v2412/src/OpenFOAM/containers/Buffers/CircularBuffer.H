/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::CircularBuffer

Description
    A simple list of objects of type \<T\> that is intended to be used
    as a circular buffer (eg, a FIFO) when the alloc/free overhead
    associated with a linked-list approach is to be avoided.

    The internal storage is addressed by independent begin/end markers.
    - The %begin marker points to the \em front.
    - The %end marker is a one-past the \em back.
    .
    This results in a variety ofr different possible buffer states:
    -# \em empty (\em begin == \em end)

    -# \em simple/linear (\em begin \< \em end) has no wrapping:
       \verbatim
       |.|.|.|a|b|c|d|.|.|.|
       beg ___^
       end ___________^
       \endverbatim

    -# \em split (\em begin \> \em end):
       \verbatim
       |f|g|h|i|.|.|.|a|b|c|d|e|
       end _____^
       beg ___________^
       \endverbatim
    .

    The methods range_one(), range_two() return the internal indexing and
    the methods array_one(), array_two() provide direct access to the
    internal contents.

    When filling the buffer, the internal storage will be resized
    (doubling strategy) as required. When this occurs, the new list
    will be linearized with \em begin = 0.

    Simultaneously when filling, the storage buffer will be over-allocated
    to avoid ambiguity when (\em begin == \em end), which represents an
    \em %empty buffer and not a \em %full buffer. Eg,
    \verbatim
        |c|d|.|a|b|
        end _^
        beg ___^
    \endverbatim
    after appending one more, it would be incorrect to simply fill
    the available space:
    \verbatim
        |c|d|e|a|b|
        end ___^        WRONG : would represent empty!
        beg ___^
    \endverbatim
    the storage is instead increased (doubled) and rebalanced before
    the append occurs (old capacity 5, new capacity 10):
    \verbatim
        |a|b|c|d|e|.|.|.|.|.|
        _^_ beg
        end _______^
    \endverbatim

SourceFiles
    CircularBuffer.C
    CircularBufferI.H
    CircularBufferIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_CircularBuffer_H
#define Foam_CircularBuffer_H

#include "labelRange.H"
#include "List.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class CircularBuffer;

template<class T>
Istream& operator>>(Istream& is, CircularBuffer<T>& rhs);

template<class T>
Ostream& operator<<(Ostream& os, const CircularBuffer<T>& rhs);


/*---------------------------------------------------------------------------*\
                       Class CircularBuffer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CircularBuffer
{
    // Private Data

        //- The allocated buffer storage
        List<T> storage_;

        //- The first addressable element
        label begin_;

        //- One past last addressable element
        label end_;


    // Private Member Functions

        //- Map the logical location to the buffer location
        inline label toGlobal(const label i) const;

        //- Length of array one
        inline label size_one() const noexcept;

        //- Length of array two
        inline label size_two() const noexcept;

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        //  The 'nocopy' option will not attempt to recover old content
        void doReserve(const bool nocopy, const label len);

        //- Copy all list contents
        template<class OtherListType>
        inline void copyList(const OtherListType& rhs);


public:

    // STL type definitions

        //- The value type the list contains
        typedef T value_type;

        //- The pointer type for non-const access to value_type items
        typedef T* pointer;

        //- The pointer type for const access to value_type items
        typedef const T* const_pointer;

        //- The type used for storing into value_type objects
        typedef T& reference;

        //- The type used for reading from constant value_type objects
        typedef const T& const_reference;

        //- The type to represent the size of a buffer
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;

        //- Forward iterator with const access
        class const_iterator;


    // Constructors

        //- Default construct, empty buffer without allocation
        inline constexpr CircularBuffer() noexcept;

        //- Construct an empty buffer with given reserve size
        inline explicit CircularBuffer(const label len);

        //- Copy construct
        inline CircularBuffer(const CircularBuffer<T>& list);

        //- Move construct
        inline CircularBuffer(CircularBuffer<T>&& list);

        //- Construct from Istream - uses readList
        explicit CircularBuffer(Istream& is);


    // Member Functions

    // Characteristics

        //- Lower capacity limit
        static constexpr label min_size() noexcept { return 16; }

        //- Size of the underlying storage.
        inline label capacity() const noexcept;

        //- Empty or exhausted buffer
        inline bool empty() const noexcept;

        //- The current number of buffer items
        inline label size() const noexcept;


    // Internal Access

        //- The nominal space available to fill.
        //- Subtract 1 for the number to append before re-balancing is needed.
        inline label space() const noexcept;

        //- The addressing range covered by array_one()
        inline labelRange range_one() const noexcept;

        //- The addressing range covered by array_two()
        inline labelRange range_two() const noexcept;

        //- The contents of the first internal array
        SubList<T> array_one();

        //- The contents of the second internal array
        SubList<T> array_two();

        //- The contents of the first internal array
        const SubList<T> array_one() const;

        //- The contents of the second internal array
        const SubList<T> array_two() const;


    // Access

        //- Access the first element (front). Requires !empty().
        T& front();

        //- Access the last element (back). Requires !empty().
        T& back();

        //- Const access to the first element (front). Requires !empty().
        const T& front() const;

        //- Const access to the last element (back). Requires !empty().
        const T& back() const;


    // Sizing

        //- Reserve allocation space for at least this size, allocating new
        //- space if required and \em retaining old content.
        //  Never shrinks.
        inline void reserve(const label len);

        //- Reserve allocation space for at least this size, allocating new
        //- space if required \em without retaining old content.
        //  Never shrinks.
        inline void reserve_nocopy(const label len);

        //- Clear the addressed buffer, does not change allocation
        inline void clear() noexcept;

        //- Clear the buffer and delete storage.
        inline void clearStorage();

        //- Swap content, independent of sizing parameter
        inline void swap(CircularBuffer<T>& other);


    // Search

        //- True if the value is contained in the list.
        inline bool contains(const T& val) const;

        //- Is the value contained in the list?
        //  \param val The value to search for
        //  \param pos The first position to examine (no-op if -ve)
        //  \return true if found.
        inline bool contains(const T& val, label pos) const;

        //- Find index of the first occurrence of the value.
        //  Any occurrences before the start pos are ignored.
        //  Linear search.
        //  \return position in list or -1 if not found.
        label find(const T& val, label pos = 0) const;


    // Stack-like Operations

        //- Copy prepend an element to the front of the buffer
        inline void push_front(const T& val);

        //- Move prepend an element to the front of the buffer
        inline void push_front(T&& val);

        //- Construct an element at the front of the buffer,
        //- return reference to the new element
        template<class... Args>
        inline T& emplace_front(Args&&... args);

        //- Copy append an element to the end of the buffer
        inline void push_back(const T& val);

        //- Move append an element to the end of the buffer
        inline void push_back(T&& val);

        //- Construct an element at the end of the buffer,
        //- return reference to the new element
        template<class... Args>
        inline T& emplace_back(Args&&... args);

        //- Shrink by moving the front of the buffer 1 or more times
        inline void pop_front(label n = 1);

        //- Shrink by moving the end of the buffer 1 or more times
        inline void pop_back(label n = 1);

        //- Copy append multiple elements the end of the buffer
        inline void push_back(const UList<T>& list);

        //- Copy append IndirectList elements the end of the buffer
        template<class Addr>
        inline void push_back(const IndirectListBase<T, Addr>& list);

        //- Append an element if not already in the buffer.
        //  \return the change in the buffer length
        inline label push_uniq(const T& val);


    // Other Operations

        //- Return a copy of the buffer flattened into a single List.
        //- Use sparingly!
        List<T> list() const;

        //- Reverse the buffer order, swapping elements
        void reverse();


    // Member Operators

        //- Non-const access to an element in the list.
        //  The index is allowed to wrap in both directions
        inline T& operator[](const label i);

        //- Const access to an element in the list
        //  The index is allowed to wrap in both directions
        inline const T& operator[](const label i) const;

        //- Copy construct
        inline void operator=(const CircularBuffer<T>& list);

        //- Move construct
        inline void operator=(CircularBuffer<T>&& list);

        //- Assign all addressed elements to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const Foam::zero);

        //- Deep copy values from a list of the addressed elements
        inline void operator=(const UList<T>& rhs);

        //- Deep copy values from a list of the addressed elements
        template<class AnyAddr>
        inline void operator=(const IndirectListBase<T, AnyAddr>& rhs);


    // IOstream Operators

        //- Print information
        Ostream& info(Ostream& os) const;

        //- Read buffer contents from Istream.
        Istream& readList(Istream& is);

        //- Write buffer contents with line-breaks in ASCII
        //- when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;

        //- Use the readList() method to read contents from Istream.
        friend Istream& operator>> <T>
        (
            Istream& is,
            CircularBuffer<T>& list
        );

        //- Write to Ostream
        friend Ostream& operator<< <T>
        (
            Ostream& os,
            const CircularBuffer<T>& list
        );


    // Iterators

        //- A simple forward const iterator for a circular buffer
        class const_iterator
        {
            const CircularBuffer<T>* container_;
            label iter_;

        public:

            using difference_type = label;
            using value_type = const T;
            using pointer = const T*;
            using reference = const T&;
            using iterator_category = std::forward_iterator_tag;

            const_iterator(const const_iterator&) = default;
            const_iterator& operator=(const const_iterator&) = default;

            const_iterator
            (
                const CircularBuffer<T>* buffer,
                label i
            )
            :
                container_(buffer),
                iter_(i)
            {}

            reference operator*() const
            {
                return (*container_)[iter_];
            }

            const_iterator& operator++()
            {
                ++iter_;
                return *this;
            }

            const_iterator operator++(int)
            {
                auto old(*this);
                ++iter_;
                return old;
            }

            bool operator==(const const_iterator& rhs) const
            {
                return iter_ == rhs.iter_;
            }

            bool operator!=(const const_iterator& rhs) const
            {
                return iter_ != rhs.iter_;
            }
        };


    // Iterator (const)

        //- Return a const_iterator at begin of buffer
        inline const_iterator cbegin() const
        {
            return const_iterator(this, 0);
        }

        //- Return a const_iterator at end of buffer
        inline const_iterator cend() const
        {
            return const_iterator(this, this->size());
        }

        //- Return a const_iterator at begin of buffer
        inline const_iterator begin() const { return cbegin(); }

        //- Return a const_iterator at end of buffer
        inline const_iterator end() const { return cend(); }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read buffer contents from Istream
template<class T>
Istream& operator>>(Istream& is, CircularBuffer<T>& rhs)
{
    return rhs.readList(is);
}


//- Write buffer contents to Ostream,
//- as per CircularBuffer::writeList() with default length.
template<class T>
Ostream& operator<<(Ostream& os, const CircularBuffer<T>& rhs)
{
    return rhs.writeList(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CircularBufferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CircularBuffer.C"
    #include "CircularBufferIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
