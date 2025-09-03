/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::sliceRange

Description
    A set of labels defined by a start, a length and a stride.

SourceFiles
    sliceRange.C
    sliceRangeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_sliceRange_H
#define Foam_sliceRange_H

#include "labelFwd.H"
#include <iterator>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
template<class T> class List;
template<class T, unsigned N> class FixedList;

/*---------------------------------------------------------------------------*\
                         Class sliceRange Declaration
\*---------------------------------------------------------------------------*/

class sliceRange
{
    // Private Data

        //- The start of the interval
        label start_;

        //- The length of the interval
        label size_;

        //- The stride for the interval
        label stride_;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef label value_type;

        //- The type that can represent the size of the range
        typedef label size_type;

        //- Input iterator with const access
        class const_iterator;

        //- Reverse input iterator with const access
        class const_reverse_iterator;


    // Generated Methods: copy/move construct, copy/move assignment


    // Constructors

        //- Default construct an empty slice (0,0,0)
        inline constexpr sliceRange() noexcept;

        //- Construct slice from start/size/stride, no checks
        inline constexpr sliceRange
        (
            const label beg,
            const label len,
            const label stride
        ) noexcept;

        //- Construct slice from start/size/stride coefficients,
        //- enforce non-negative size and stride.
        explicit sliceRange(const FixedList<label,3>& coeffs);


    // Member Functions

    // Access

        //- True if range is empty (zero-sized)
        bool empty() const noexcept { return !size_; }

        //- True if range has size greater than zero
        bool good() const noexcept { return (size_ > 0); }

        //- The size of the range
        label size() const noexcept { return size_; }

        //- Non-const access to size of the range
        label& size() noexcept { return size_; }

        //- The (inclusive) lower value of the range
        label start() const noexcept { return start_; }

        //- Non-const access to start of the range
        label& start() noexcept { return start_; }

        //- The stride for the range
        label stride() const noexcept { return stride_; }

        //- Non-const access to the stride
        label& stride() noexcept { return stride_; }

        //- The (inclusive) lower value of the range,
        //- same as start()
        label min() const noexcept { return start_; }

        //- The (inclusive) upper value of the range.
        //- Ill-defined for an empty range
        label max() const noexcept { return start_ + (size_-1) * stride_; }

        //- Return list of labels corresponding to the slice
        List<label> labels() const;


    // Member Operators

        //- Return element in the range, without bounds checking
        label operator[](const label i) const noexcept
        {
            return start_ + i * stride_;
        }

        //- True if range has size greater than zero. Same as good()
        explicit operator bool() const noexcept { return (size_ > 0); }


    // Iteration / Generation

        //- A value indexer, for iteration or generation
        class indexer;

        //- Return a forward values generator
        inline indexer generator() const;

        //- Return const_iterator to a position within the range,
        //- with bounds checking.
        //  \return iterator at the requested position, or end() for
        //      out-of-bounds
        inline const_iterator at(const label i) const;

        //- A const_iterator set to the beginning of the range
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to the beginning of the range
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator cend() const noexcept;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator end() const noexcept;

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator rbegin() const noexcept;

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator crbegin() const noexcept;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator rend() const noexcept;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator crend() const noexcept;


    // Iterators

        //- A value indexer, for iteration or generation
        class indexer
        {
            //- The stride when indexing
            label stride_;

            //- The global value
            label value_;

        public:

            // STL definitions (as per std::iterator)
            typedef label value_type;
            typedef label difference_type;
            typedef const label* pointer;
            typedef label reference;


        // Constructors

            //- Default construct with zero value and stride = 1
            inline constexpr indexer() noexcept;

            //- Construct with specified value and stride
            inline constexpr indexer
            (
                const label val,
                const label stride
            ) noexcept;


        // Member Functions

            //- The current value
            constexpr label value() const noexcept { return value_; }

            //- The stride
            constexpr label stride() const noexcept { return stride_; }

            //- Value with offset
            constexpr label value(const label n) const noexcept
            {
                return value_ + (n * stride_);
            }

            //- Decrement value
            void prev() noexcept { value_ -= stride_; }

            //- Decrease value
            void prev(const label n) noexcept { value_ -= (n * stride_); }

            //- Increment value
            void next() noexcept { value_ += stride_; }

            //- Increase value
            void next(const label n) noexcept { value_ += (n * stride_); }


        // Member Operators

            //- Return the value
            constexpr label operator*() const noexcept { return value_; }

            //- Apply a postfix increment and return the current value.
            //  This operator definition is required for a generator -
            //  see std::generate()
            inline label operator()() noexcept;
        };


        //- Bidirectional input iterator with const access
        class const_iterator
        :
            public indexer
        {
        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;


        // Constructors

            //- Inherit constructors from indexer
            using indexer::indexer;


        // Member Operators

            //- Offset dereference operator
            inline constexpr label operator[](const label n) const noexcept;

            //- Prefix increment
            inline const_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_iterator operator++(int) noexcept;

            //- Prefix decrement
            inline const_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_iterator operator--(int) noexcept;

            //- Arbitrary increment
            inline const_iterator& operator+=(const label n) noexcept;

            //- Arbitrary decrement
            inline const_iterator& operator-=(const label n) noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator+
            (
                const label n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator-
            (
                const label n
            ) const noexcept;

            //- Difference operator
            inline constexpr label operator-
            (
                const const_iterator& iter
            ) const noexcept;


        // Comparison

            //- Test for equality of values (ignore stride)
            inline constexpr bool operator==(const const_iterator& iter)
            const noexcept;

            //- Compare less-than values (ignore stride)
            inline constexpr bool operator<(const const_iterator& iter)
            const noexcept;


        // Derived comparisons

            constexpr bool operator!=(const const_iterator& iter)
            const noexcept
            {
                return !(*this == iter);
            }

            constexpr bool operator<=(const const_iterator& iter)
            const noexcept
            {
                return !(iter < *this);
            }

            constexpr bool operator>(const const_iterator& iter)
            const noexcept
            {
                return (iter < *this);
            }

            constexpr bool operator>=(const const_iterator& iter)
            const noexcept
            {
                return !(*this < iter);
            }
        };


        //- Bidirectional reverse input iterator with const access
        class const_reverse_iterator
        :
            public indexer
        {
        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;


        // Constructors

            //- Inherit constructors from indexer
            using indexer::indexer;


        // Member Operators

            //- Offset dereference operator
            inline constexpr label operator[](const label n) const noexcept;

            //- Prefix increment
            inline const_reverse_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_reverse_iterator operator++(int) noexcept;

            //- Prefix decrement
            inline const_reverse_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_reverse_iterator operator--(int) noexcept;

            //- Arbitrary increment
            inline const_reverse_iterator& operator+=(const label n) noexcept;

            //- Arbitrary decrement
            inline const_reverse_iterator& operator-=(const label n) noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator+
            (
                const label n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator-
            (
                const label n
            ) const noexcept;

            //- Difference operator
            inline constexpr label operator-
            (
                const const_reverse_iterator& iter
            ) const noexcept;


        // Comparison

            //- Test for equality of values (ignore stride)
            inline constexpr bool operator==(const const_reverse_iterator& iter)
            const noexcept;

            //- Reverse compare less-than values (ignore stride)
            inline constexpr bool operator<(const const_reverse_iterator& iter)
            const noexcept;


        // Derived comparisons

            constexpr bool operator!=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(*this == iter);
            }

            constexpr bool operator<=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(iter < *this);
            }

            constexpr bool operator>(const const_reverse_iterator& iter)
            const noexcept
            {
                return (iter < *this);
            }

            constexpr bool operator>=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(*this < iter);
            }
        };
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read sliceRange from Istream as (start size stride) tuple, no checks
Istream& operator>>(Istream& is, sliceRange& range);

//- Write sliceRange to Ostream as (start size stride) tuple
Ostream& operator<<(Ostream& os, const sliceRange& range);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Test for equality of begin/size/slice values
inline bool operator==
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return
    (
        a.start() == b.start()
     && a.size() == b.size()
     && a.stride() == b.stride()
    );
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, compare the sizes.
//  If the sizes are also equal, compare the strides.
inline bool operator<
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return
    (
        a.start() < b.start()
     ||
        (
            !(b.start() < a.start())
         &&
            (
                a.size() < b.size()
             ||
                (
                    !(b.size() < a.size())
                 && a.stride() < b.stride()
                )
            )
        )
    );
}


// Derived comparisons

inline bool operator!=
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return !(a == b);
}


inline bool operator<=
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return !(b < a);
}

inline bool operator>
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return (b < a);
}

inline bool operator>=
(
    const sliceRange& a,
    const sliceRange& b
) noexcept
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sliceRangeI.H"

#endif

// ************************************************************************* //
