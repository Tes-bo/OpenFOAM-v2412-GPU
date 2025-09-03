/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::IntRange

Description
    An interval of (signed) integers defined by a start and a size.

Note
    Only a minimum of IO operators are defined, to avoid incurring
    too many dependencies or cyclic dependencies.

SourceFiles
    IntRangeI.H
    IntRanges.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IntRange_H
#define Foam_IntRange_H

#include "labelFwd.H"
#include <iterator>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
template<class T> class List;

/*---------------------------------------------------------------------------*\
                          Class IntRange Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class IntRange
{
    static_assert(std::is_integral<IntType>::value, "Integral required");

    // Private Data

        //- The start of the interval
        IntType start_;

        //- The length of the interval
        IntType size_;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef IntType value_type;

        //- The type that can represent the size of the range
        typedef IntType size_type;

        //- Input iterator with const access
        class const_iterator;

        //- Reverse input iterator with const access
        class const_reverse_iterator;


    // Public Classes

        //- Unary predicate for greater than 0 (int values)
        struct gt0
        {
            template<class Int>
            constexpr bool operator()(const Int val) const noexcept
            {
                return (val > 0);
            }
        };

        //- Unary predicate for less than 0 (int values)
        struct lt0
        {
            template<class Int>
            constexpr bool operator()(const Int val) const noexcept
            {
                return (val < 0);
            }
        };

        //- Unary predicate for greater-equal 0 (int values)
        struct ge0
        {
            template<class Int>
            constexpr bool operator()(const Int val) const noexcept
            {
                return (val >= 0);
            }
        };

        //- Unary predicate for less-equal 0 (int values)
        struct le0
        {
            template<class Int>
            constexpr bool operator()(const Int val) const noexcept
            {
                return (val <= 0);
            }
        };


    // Generated Methods: copy/move construct, copy/move assignment


    // Constructors

        //- Default construct an empty range (0,0)
        inline constexpr IntRange() noexcept;

        //- Construct a range with specified length, starting at zero (0,len)
        inline explicit constexpr IntRange(const IntType len) noexcept;

        //- Construct a range from start/length, no checks
        inline constexpr IntRange
        (
            const IntType beg,
            const IntType len
        ) noexcept;


    // Member Functions

    // Access

        //- True if range is empty (zero-sized)
        bool empty() const noexcept { return !size_; }

        //- True if range has size greater than zero
        bool good() const noexcept { return (size_ > 0); }

        //- The size of the range
        IntType size() const noexcept { return size_; }

        //- Non-const access to size of the range
        IntType& size() noexcept { return size_; }

        //- The (inclusive) lower value of the range
        IntType start() const noexcept { return start_; }

        //- Non-const access to start of the range
        IntType& start() noexcept { return start_; }

        //- The (inclusive) lower value of the range,
        //- same as start(), begin_value()
        IntType min() const noexcept { return start_; }

        //- The (inclusive) upper value of the range,
        //- same as rbegin_value(). Ill-defined for an empty range
        IntType max() const noexcept { return (start_ + (size_ - 1)); }


    // Edit

        //- Reset to zero start and zero size
        inline void reset() noexcept;

        //- Reset start and length, no checks
        inline void reset(const IntType beg, const IntType len) noexcept;

        //- Same as reset() - reset to zero start and zero size
        void clear() noexcept { reset(); }

        //- Set the start position, no checks
        inline void setStart(const IntType i) noexcept;

        //- Change the size, no checks. Identical to resize()
        inline void setSize(const IntType n) noexcept;

        //- Change the size, no checks. Identical to setSize()
        inline void resize(const IntType n) noexcept;

        //- Enforce non-negative size
        inline void clampSize() noexcept;


    // Search

        //- True if the (global) value is within the range
        inline bool contains(const IntType value) const noexcept;

        //- True if the (global) value is within the range
        bool found(const IntType val) const noexcept { return contains(val); }


    // Member Operators

        //- Return const_iterator to a position within the range,
        //- with bounds checking.
        //  \return iterator at the requested position, or end() for
        //      out-of-bounds
        inline const_iterator at(const IntType i) const;

        //- Offset dereference, without bounds checking
        inline constexpr IntType operator[](const IntType i) const noexcept;

        //- True if the global value is located within the range.
        //  Behaviour identical to contains() - usable as a predicate
        bool operator()(const IntType i) const noexcept { return contains(i); }

        //- Increase the size by 1.
        inline IntType operator++() noexcept;
        inline IntType operator++(int) noexcept;

        //- Increase the size by n.
        inline IntType operator+=(const IntType n) noexcept;

        //- Decrease the size by 1, but never below 0.
        inline IntType operator--() noexcept;
        inline IntType operator--(int) noexcept;

        //- Decrease the size by n, but never below 0.
        inline IntType operator-=(const IntType n) noexcept;

        //- True if range has size greater than zero. Same as good()
        explicit operator bool() const noexcept { return (size_ > 0); }


    // Iterator ranges

        //- The value at the beginning of the range - same as min(), start()
        inline IntType begin_value() const noexcept;

        //- The value 1 beyond the end of the range.
        inline IntType end_value() const noexcept;

        //- The max value of the range.
        inline IntType rbegin_value() const noexcept;

        //- The value 1 before the begin of range
        inline IntType rend_value() const noexcept;


    // Bidirectional input iterators (const)

        //- Random-access input iterator with const access
        //  \note No operator+(IntType, iterator) since this provokes
        //      misleading resolution errors
        class const_iterator
        {
            //- The global value
            IntType value_;

        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;
            typedef IntType value_type;
            typedef IntType difference_type;
            typedef const IntType* pointer;
            typedef IntType reference;


        // Constructors

            //- Construct with specified value, or default construct
            explicit constexpr const_iterator(IntType val = 0) noexcept
            :
                value_(val)
            {}


        // Member Operators

            //- Return the value
            constexpr IntType operator*() const noexcept { return value_; }

            //- Offset dereference operator
            inline constexpr IntType operator[](const IntType n) const noexcept;

            //- Prefix increment
            inline const_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_iterator operator++(int) noexcept;

            //- Prefix decrement
            inline const_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_iterator operator--(int) noexcept;

            //- Arbitrary increment
            inline const_iterator& operator+=(const IntType n) noexcept;

            //- Arbitrary decrement
            inline const_iterator& operator-=(const IntType n) noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator+
            (
                const IntType n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator-
            (
                const IntType n
            ) const noexcept;

            //- Difference operator
            inline constexpr IntType operator-
            (
                const const_iterator& iter
            ) const noexcept;


        // Comparison

            //- Test for equality of values
            inline constexpr bool operator==(const const_iterator& iter)
            const noexcept;

            //- Compare less-than
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


        //- Random-access reverse input iterator with const access
        //  \note No operator+(IntType, iterator) since this provokes
        //      misleading resolution errors
        class const_reverse_iterator
        {
            //- The global value
            IntType value_;

        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;
            typedef IntType value_type;
            typedef IntType difference_type;
            typedef const IntType* pointer;
            typedef IntType reference;


        // Constructors

            //- Construct with specified value, or default construct
            explicit constexpr const_reverse_iterator(IntType val = 0) noexcept
            :
                value_(val)
            {}


        // Member Operators

            //- Return the value
            constexpr IntType operator*() const noexcept { return value_; }

            //- Offset dereference operator
            inline constexpr IntType operator[](const IntType n) const noexcept;

            //- Prefix increment
            inline const_reverse_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_reverse_iterator operator++(int) noexcept;

            //- Prefix decrement
            inline const_reverse_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_reverse_iterator operator--(int) noexcept;

            //- Arbitrary increment
            inline const_reverse_iterator& operator+=(const IntType n) noexcept;

            //- Arbitrary decrement
            inline const_reverse_iterator& operator-=(const IntType n) noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator+
            (
                const IntType n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator-
            (
                const IntType n
            ) const noexcept;

            //- Difference operator
            inline constexpr IntType operator-
            (
                const const_reverse_iterator& iter
            )const noexcept;


        // Comparison

            //- Test for equality of values
            inline constexpr bool operator==(const const_reverse_iterator& iter)
            const noexcept;

            //- Reverse compare less-than
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


    // Bidirectional input iterators (const)

        //- A const_iterator set to the beginning of the range
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to the beginning of the range
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator cend() const noexcept;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator end() const noexcept;


    // Bidirectional reverse input iterators (const)

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator rbegin() const noexcept;

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator crbegin() const noexcept;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator rend() const noexcept;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator crend() const noexcept;
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Identity function for common integer types

//- Identity map from an int32_t IntRange
List<label> identity(const IntRange<int32_t>& range);

#if defined(WM_LABEL_SIZE) && (WM_LABEL_SIZE >= 64)
//- Identity map from an int64_t IntRange
List<label> identity(const IntRange<int64_t>& range);
#endif


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Input operators for common integer types

//- Read IntRange from Istream as bracketed (start size) tuple, no checks
Istream& operator>>(Istream& os, IntRange<int32_t>& range);

//- Read IntRange from Istream as bracketed (start size) tuple, no checks
Istream& operator>>(Istream& os, IntRange<int64_t>& range);


// Output operators for common integer types

//- Write IntRange to Ostream as bracketed (start size) tuple
Ostream& operator<<(Ostream& os, const IntRange<int32_t>& range);

//- Write IntRange to Ostream as bracketed (start size) tuple
Ostream& operator<<(Ostream& os, const IntRange<int64_t>& range);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Test for equality of begin/size values
template<class IntType>
inline constexpr bool operator==
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return (a.start() == b.start() && a.size() == b.size());
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size.
template<class IntType>
inline constexpr bool operator<
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return
    (
        a.start() < b.start()
     ||
        (
            !(b.start() < a.start())
         && a.size() < b.size()
        )
    );
}


// Derived comparisons

template<class IntType>
inline constexpr bool operator!=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return !(a == b);
}

template<class IntType>
inline constexpr bool operator<=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return !(b < a);
}

template<class IntType>
inline constexpr bool operator>
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return (b < a);
}

template<class IntType>
inline constexpr bool operator>=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept

{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IntRangeI.H"

#endif

// ************************************************************************* //
