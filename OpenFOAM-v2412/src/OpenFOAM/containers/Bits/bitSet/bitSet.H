/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::bitSet

Description
    A bitSet stores bits (elements with only two states) in packed internal
    format and supports a variety of bit-set operations.
    Its behaviour is largely list-like, with some HashSet features.

SourceFiles
    bitSetI.H
    bitSet.C
    bitSetIO.C
    bitSetTemplates.C

See also
    Foam::BitOps
    Foam::PackedList

\*---------------------------------------------------------------------------*/

#ifndef Foam_bitSet_H
#define Foam_bitSet_H

#include "className.H"
#include "PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                           Class bitSet Declaration
\*---------------------------------------------------------------------------*/

class bitSet
:
    public PackedList<1>
{
protected:

    // Logic/Set Operations

        //- The set difference
        //  \code
        //      A = (A - B)
        //      A = (A & !B)
        //      A = (A & ~B)
        //  \endcode
        //  A and B can have different sizes.
        //  Never changes the original set size.
        //  Non-overlapping parts are considered \em off.
        bitSet& minusEq(const bitSet& other);

        //- The set logical AND
        //  \code
        //      A = (A & B)
        //
        //  \endcode
        //  A and B can have different sizes.
        //  Never changes the original set size.
        //  Non-overlapping parts are considered \em off.
        bitSet& andEq(const bitSet& other);

        //- The set logical OR
        //  \code
        //      A = (A | B)
        //  \endcode
        //  A and B can have different sizes
        //
        //  The size grows to accommodate new \em on bits.
        bitSet& orEq(const bitSet& other);

        //- The set logical XOR
        //  \code
        //      A = (A ^ B)
        //  \endcode
        //  A and B can have different sizes.
        //
        //  The size grows to accommodate new \em on bits.
        bitSet& xorEq(const bitSet& other);


public:

    // Forward declaration of access classes

        class reference;
        class const_iterator;
        typedef unsigned int const_reference;


    // Static Member Functions

        //- Return a null bitSet (reference to a nullObject).
        static const bitSet& null() noexcept
        {
            return NullObjectRef<bitSet>();
        }


    //- Declare type-name (with debug switch)
    ClassName("bitSet");


    // Constructors

        //- Default construct an empty, zero-sized bitSet
        inline constexpr bitSet() noexcept;

        //- Construct from Istream
        explicit bitSet(Istream& is);

        //- Construct with given size, with all bits set to 0
        inline explicit bitSet(const label n);

        //- Construct with given size and value for all elements
        inline bitSet(const label n, const bool val);

        //- Copy construct
        inline bitSet(const bitSet& bitset);

        //- Move construct
        inline bitSet(bitSet&& bitset);

        //- Construct a new bitSet by extracting the specified (unique)
        //- locations of an existing bitSet.
        bitSet(const bitSet& bitset, const labelUList& addr);

        //- Construct a new set by extracting the specified (unique)
        //- locations of an existing bitSet.
        template<class Addr>
        bitSet
        (
            const bitSet& bitset,
            const IndirectListBase<label, Addr>& addr
        );

        //- Construct a new bitSet by extracting the specified range
        //- locations of an existing bitSet.
        bitSet(const bitSet& bitset, const labelRange& range);

        //- Construct from a list of bools
        inline explicit bitSet(const UList<bool>& bools);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        explicit bitSet(const label n, const labelRange& range);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        inline bitSet(const label n, const labelUList& locations);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        template<class Addr>
        bitSet
        (
            const label n,
            const IndirectListBase<label, Addr>& locations
        );

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        template<unsigned N>
        bitSet(const label n, const FixedList<label, N>& locations);

        //- Construct with given pre-size (filled with 0),
        //- subsequently add specified locations as 1,
        //- auto-vivify entries if needed.
        inline bitSet(const label n, std::initializer_list<label> locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        explicit bitSet(const labelRange& range);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        inline explicit bitSet(const labelUList& locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        template<class Addr>
        inline explicit bitSet(const IndirectListBase<label, Addr>& locations);

        //- Construct with automatic sizing (filled with 0),
        //- and set the specified locations as 1.
        template<unsigned N>
        explicit bitSet(const FixedList<label, N>& locations);

        //- Clone
        inline autoPtr<bitSet> clone() const;


    // Member Functions

    // Query

        //- True if all bits in this bitset are set or if the set is \b empty.
        //  Returning true for an empty set may not seem intuitive, but
        //  conforms with boost definitions and std::all_of behaviour.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool all() const;

        //- True if any bits in this bitset are set.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool any() const;

        //- True if no bits in this bitset are set.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline bool none() const;

        //- Count number of bits set.
        //  \param on can be set to false to count the number of unset bits
        //     instead.
        //  \note Method name compatibility with boost::dynamic_bitset
        inline unsigned int count(const bool on=true) const;

        //- True if any bits in the other bitset intersect (are the same).
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        bool intersects(const bitSet& other) const;

        //- Test for \em True value at specified position,
        //- never auto-vivify entries.
        //
        //  \note Method name compatibility with std::bitset
        bool test(const label pos) const { return this->get(pos); }

        //- Test for \em True value at specified position,
        //- never auto-vivify entries.
        //
        //  \note Method name compatibility with HashSet
        bool contains(const label pos) const { return this->get(pos); }

        //- Locate the first bit that is set.
        //  \return the location or -1 if there are no bits set.
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        inline label find_first() const;

        //- Locate the first bit that is unset.
        //  \return the location or -1 if the set is empty or all bits are on.
        //
        //  \note Provided for symmetry with find_first()
        inline label find_first_not() const;

        //- Locate the last bit set.
        //  \return the location or -1 if there are no bits set.
        //
        //  \note Provided for symmetry with find_first()
        inline label find_last() const;

        //- Locate the next bit set, starting one beyond the specified position
        //  \return the location or -1 if there are no further bits set.
        //
        //  \note Method name compatibility with boost::dynamic_bitset
        inline label find_next(label pos) const;

        //- The indices of the \a on bits as a sorted labelList.
        //
        //  \note Method name compatibility with HashSet
        labelList toc() const;

        //- The indices of the \a on bits as a sorted labelList.
        //  This is identical to toc(), which is always sorted.
        //
        //  \note Method name compatibility with HashSet
        inline labelList sortedToc() const;

        //- Return the bitset values as a boolList.
        //  When the output is a bool, this is more efficient than unpack()
        List<bool> values() const;


    // Assignment

        //- Assign all entries to the given value.
        inline void fill(const bool val);

        //- Copy assign all entries from a list of bools.
        void assign(const UList<bool>& bools);


    // Setting single or multiple values

        //- Assign a single index/value
        using PackedList<1>::set;

        //- Set specified bits from another bitset.
        //  The current set size may grow to accommodate any new bits
        //  (auto-vivifies).
        inline void set(const bitSet& bitset);

        //- Set the specified range of bits
        //  The current set size may grow to accommodate any new bits
        //  (auto-vivifies).
        //  \note this operation is generally more efficient than calling
        //      set(pos) on individual bits.
        void set(const labelRange& range);


    // Unsetting single or multiple values

        //- Unset a single index
        using PackedList<1>::unset;

        //- Unset (subtract) the bits specified in the other bitset, which is
        //- a set difference corresponds to the logical operation
        //  \code
        //      A = (A & !B)
        //  \endcode
        //  The result is comparable to 'operator-='
        //  \endcode
        //
        //  A and B can have different sizes.
        //  Does not change the original set size.
        inline bitSet& unset(const bitSet& other);

        //- Unset the specified range of bits specified, never auto-vivifies.
        //  \note this operation can be more efficient than calling
        //      unset(pos) on individual bits.
        void unset(const labelRange& range);


    // Edit

        //- Invert all bits in the addressable region
        inline void flip();

        //- Invert bit at the specified position.
        //  A no-op if the position is out-of-range
        inline void flip(const label pos);

        //- Resize to include the last \em on bit only.
        //  Functionally identical to resize(find_last()+1)
        inline void resize_last();

        //- Swap contents
        inline void swap(bitSet& bitset);

        //- Transfer the contents of the argument list into this list
        //- and annul the argument list.
        inline void transfer(bitSet& bitset);


    // Convenience methods

        //- Ensure the addressable range does not exceed maxSize.
        //  Either decreases the size of the bitSet or is a no-op.
        //
        //  \code
        //      pointField& pts = mesh.points();
        //      bitset.bound(pts.size());
        //
        //      for (const label pointi : bitset)
        //      {
        //          pts[pointi]  ...
        //      }
        //  \endcode
        inline bitSet& bound(const label maxSize);

        //- Ensure the addressable range does not exceed that of other.
        //  Either decreases the size of the bitSet or is a no-op.
        inline bitSet& bound(const bitSet& other);

        //- Ensure that minSize is covered by the bitSet.
        //  Either increases the size of the bitSet or is a no-op.
        inline bitSet& extend(const label minSize);

        //- Ensure the bitset is addressable throughout the range of other.
        //  Either increases the size of the bitSet or is a no-op.
        inline bitSet& extend(const bitSet& other);

        //- Set the locations listed by the iterator range,
        //- auto-vivify entries if needed.
        //
        //  \return number of locations changed
        template<class InputIter>
        label setMany(InputIter first, InputIter last);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        inline label set(const labelUList& locations);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        template<class Addr>
        inline label set(const IndirectListBase<label, Addr>& locations);

        //- Set the listed locations to 1.
        //  Does auto-vivify for non-existent entries.
        //
        //  \return number of locations changed
        template<unsigned N>
        label set(const FixedList<label, N>& locations);

        //- Unset the locations listed by the iterator range,
        //- never auto-vivify entries.
        //
        //  \return number of locations changed
        template<class InputIter>
        label unset(InputIter first, InputIter last);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        inline label unset(const labelUList& locations);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        template<class Addr>
        inline label unset(const IndirectListBase<label, Addr>& locations);

        //- Unset the listed locations, never auto-vivifies.
        //
        //  \return number of locations changed
        template<unsigned N>
        label unset(const FixedList<label, N>& locations);


    // Access helpers

        //- A reference supporting read/write access to an entry
        class reference
        :
            public PackedList<1>::reference
        {
        protected:

            friend class bitSet;        // Access for parent
            void operator&() = delete;  // Refuse to provide its address

            //- Construct by taking reference of block from within
            //- the list and the specified index.
            inline reference(bitSet* parent, const label index);

        public:

            //- Copy construct
            reference(const reference&) = default;

            //- Move construct
            reference(reference&&) = default;

            //- Flip the bit at the position, no range-checking
            inline void flip();

            //- Value assignment
            inline void operator=(const reference& other);

            //- Value assignment
            inline void operator=(const unsigned int val);

            //- Conversion operator
            inline operator unsigned int () const;
        };


    // Iteration

        //- A const_iterator for iterating across \a on values
        class const_iterator
        {
            friend class bitSet;

            //- The parent being iterated
            const bitSet* set_;

            //- Global position of the current \a on bit
            label pos_;

            //- Default construct - an end iterator
            inline const_iterator() noexcept;

            //- Construct begin iterator
            inline const_iterator(const bitSet* bitset);

            //- Construct iterator, starting at or beyond the given position
            inline const_iterator(const bitSet* bitset, label pos);

        public:

            //- Return the current \a on position
            inline label operator*() const noexcept;

            //- Move to the next \a on position
            inline const_iterator& operator++();

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- Iterator set to the position of the first \a on bit
        inline const_iterator begin() const;

        //- Iterator set to the position of the first \a on bit
        inline const_iterator cbegin() const;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator begin(label pos) const;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator cbegin(label pos) const;

        //- Iterator beyond the end of the bitSet
        inline const_iterator end() const noexcept;

        //- Iterator beyond the end of the bitSet
        inline const_iterator cend() const noexcept;


    // Member Operators

        //- Test value at specified position, same as test()
        //  Enables use as a predicate
        inline bool operator()(const label pos) const;

        //- Identical to get() - get value at index.
        //  Never auto-vivify entries.
        inline unsigned int operator[](const label i) const;

        //- Non-const access to value at index.
        //  Fatal for out-of-range indices
        inline reference operator[](const label i);

        //- Copy assignment
        inline bitSet& operator=(const bitSet& bitset);

        //- Move assignment
        inline bitSet& operator=(bitSet&& bitset);

        //- Assign all entries to the given value. fill()
        inline bitSet& operator=(const bool val);

        //- Bitwise-AND all the bits in other with the bits in this bitset.
        //  The operands may have dissimilar sizes,
        //  never changes the original set size.
        //  Non-overlapping elements are considered \em off.
        inline bitSet& operator&=(const bitSet& other);

        //- Bitwise-OR operator - similar to the set() method.
        //  The operands may have dissimilar sizes,
        //  the set grows to accommodate new \em on bits.
        inline bitSet& operator|=(const bitSet& other);

        //- Bitwise-XOR operator - retains unique entries.
        //  The operands may have dissimilar sizes,
        //  the set grows to accommodate new \em on bits.
        inline bitSet& operator^=(const bitSet& other);

        //- Remove entries from this list - identical to the unset() method.
        //  The operands may have dissimilar sizes,
        //  never changes the original set size.
        inline bitSet& operator-=(const bitSet& other);


    // IOstream Operators

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<bitSet> info() const noexcept
        {
            return *this;
        }


    // Housekeeping

        //- Same as contains()
        bool found(const label pos) const { return this->contains(pos); }

        //- Deprecated(2020-11) use fill()
        //  \deprecated(2020-11) use fill()
        void assign(const unsigned int val) { this->fill(val); }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Hashing for bitSet data
template<> struct Hash<bitSet> : bitSet::hasher {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Write bitset to Ostream with 40 items per line.
Ostream& operator<<(Ostream& os, const bitSet& bitset);

//- Output bitset information
Ostream& operator<<(Ostream& os, const InfoProxy<bitSet>& iproxy);


//- Bitset complement, returns a copy of the bitset with all its bits flipped
inline bitSet operator~(const bitSet& bitset);

//- Bitwise-AND of two bitsets.
//  See bitSet::operator&= for more details.
inline bitSet operator&(const bitSet& a, const bitSet& b);

//- Bitwise-OR of two bitsets
//  See bitSet::operator|= for more details.
inline bitSet operator|(const bitSet& a, const bitSet& b);

//- Bitwise-XOR of two bitsets to form a unique bit-set
//  See bitSet::operator^= for more details.
inline bitSet operator^(const bitSet& a, const bitSet& b);

//- Bitwise difference (subset) of two bitsets to form a unique bit-set
//  See bitSet::operator-= for more details.
inline bitSet operator-(const bitSet& a, const bitSet& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bitSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bitSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
