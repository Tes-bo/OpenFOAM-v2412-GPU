/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::CompactListList

Description
    A packed storage of objects of type \<T\>
    using an offset table for access.

    The offset table is the size of the number of rows+1
    whose elements are the
    accumulated sizes of the rows, i.e.
      - offset[i] gives the index of first element of row i
      - offset[i+1] - offset[i] is the number of elements in row i

    Note that an empty CompactListList should have empty offsets
    (not size 1).

SourceFiles
    CompactListList.C
    CompactListListI.H
    CompactListListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_CompactListList_H
#define Foam_CompactListList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class CompactListList;

template<class T> Istream& operator>>(Istream&, CompactListList<T>&);
template<class T> Ostream& operator<<(Ostream&, const CompactListList<T>&);

// Common list types

//! A CompactListList of labels
typedef CompactListList<label> labelCompactListList;


/*---------------------------------------------------------------------------*\
                       Class CompactListList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CompactListList
{
    // Private Data

        //- Offset (addressing) table
        labelList offsets_;

        //- Packed (matrix) of values
        List<T> values_;


    // Private Member Functions

        //- Report overflow at specified (non-negative) index
        static void reportOverflowAndExit
        (
            const label idx,
            const label prevOffset = -1,  // The last valid offset value
            const label count = 0         // The count to add to prevOffset
        );

        //- Construct by packing together the list of lists
        template<class ListListType>
        static CompactListList<T> pack_impl
        (
            const ListListType& lists,
            const bool checkOverflow = false
        );

        //- Avoid poorly sized offsets/values
        inline void enforceSizeSanity();

        //- The max of localSizes, excluding the specified row
        label maxNonLocalSize(const label rowi) const;

        //- Change the offset table based on the sizes
        //- and return the total number of values to be stored but
        //- does not touch the values array (unless the list is empty)
        label resize_offsets
        (
            const labelUList& listSizes,
            const bool checkOverflow = false
        );

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

        //- The type used for reading from constant value_type objects.
        typedef const T& const_reference;

        //- The type to represent the size of a CompactListList
        typedef label size_type;


    // Static Member Functions

        //- Return a null CompactListList (reference to a nullObject).
        static const CompactListList<T>& null() noexcept
        {
            return NullObjectRef<CompactListList<T>>();
        }

        //- Construct by packing together the list of lists
        template<class SubListType = List<T>>
        static CompactListList<T> pack
        (
            const UList<SubListType>& lists,
            const bool checkOverflow = false
        );

        //- Construct by packing together an indirect list of lists
        template<class SubListType, class Addr>
        static CompactListList<T> pack
        (
            const IndirectListBase<SubListType, Addr>& lists,
            const bool checkOverflow = false
        );


    // Constructors

        //- Default construct
        CompactListList() noexcept = default;

        //- Copy construct
        inline CompactListList(const CompactListList<T>& list);

        //- Move construct
        inline CompactListList(CompactListList<T>&& list);

        //- Copy/move construct as specified.
        inline CompactListList(CompactListList<T>& list, bool reuse);

        //- Construct from number of rows and number of values.
        inline CompactListList(const label mRows, const label nVals);

        //- Construct from number of rows, number of values
        //- initializing all elements to zero
        inline CompactListList
        (
            const label mRows, const label nVals, const Foam::zero
        );

        //- Construct from number of rows, number of values
        //- and uniform value for all elements.
        inline CompactListList(const label mRows, const label nVals, const T&);

        //- Construct given list of row-sizes.
        explicit CompactListList(const labelUList& listSizes);

        //- Construct given list of row-sizes and a uniform value
        CompactListList(const labelUList& listSizes, const T& val);

        //- Construct from Istream.
        explicit CompactListList(Istream& is);

        //- Clone
        inline autoPtr<CompactListList<T>> clone() const;


    // Member Functions

    // Access

        //- True if the number of rows/sublists is zero
        inline bool empty() const noexcept;

        //- True if content is a single row/sublist only.
        //- Such content could be flattened out into a straight list
        //- (for example).
        inline bool single() const noexcept;

        //- The primary size (the number of rows/sublists)
        inline label length() const noexcept;

        //- The primary size (the number of rows/sublists)
        inline label size() const noexcept;

        //- The total addressed size, which corresponds to the
        //- end (back) offset and also the sum of all localSizes.
        inline label totalSize() const noexcept;

        //- Return the offset table (= size()+1)
        const labelList& offsets() const noexcept { return offsets_; }

        //- Return non-const access to the offset table
        labelList& offsets() noexcept { return offsets_; }

        //- Return the packed values
        const List<T>& values() const noexcept { return values_; }

        //- Return non-const access to the packed values
        List<T>& values() noexcept { return values_; }


        //- Return const pointer to the first data in values()
        inline const T* cdata() const noexcept;

        //- Return pointer to the first data in values()
        inline T* data() noexcept;

        //- Return const pointer to underlying values storage,
        //- reinterpreted as byte data.
        //  \note Only meaningful for contiguous data
        inline const char* cdata_bytes() const noexcept;

        //- Return pointer to underlying values storage,
        //- reinterpreted as byte data.
        //  \note Only meaningful for contiguous data
        inline char* data_bytes() noexcept;

        //- Number of contiguous bytes for the values data,
        //- no runtime check that the type is actually contiguous
        //  \note Only meaningful for contiguous data
        inline std::streamsize size_bytes() const noexcept;

        //- Number of contiguous bytes for the values data,
        //- runtime FatalError if type is not contiguous
        std::streamsize byteSize() const;


    // Queries

        //- The local row sizes. Same as localSizes
        inline labelList sizes() const;

        //- The local row starts
        inline const labelUList localStarts() const;

        //- The local row sizes
        labelList localSizes() const;

        //- Starting offset for given row
        inline label localStart(const label i) const;

        //- End offset (exclusive) for given row
        inline label localEnd(const label i) const;

        //- Size of given row
        inline label localSize(const label i) const;

        //- Return start/size ranges for all sub-lists
        List<labelRange> ranges() const;

        //- Return start/size range for given sub-list
        labelRange range(const label i) const;

        //- The max row length used
        inline label maxSize() const;


    // Edit

        //- Clear addressing and contents
        inline void clear();

        //- Reset size of CompactListList.
        //  \note this form only allows truncation of the CompactListList.
        inline void resize(const label mRows);

        //- Redimension CompactListList
        inline void resize(const label mRows, const label nVals);

        //- Redimension \em without preserving existing content
        inline void resize_nocopy(const label mRows, const label nVals);

        //- Redimension CompactListList and fill new elements with value.
        inline void resize(const label mRows, const label nVals, const T&);

        //- Reset dimensions of CompactListList
        void resize(const labelUList& listSizes);

        //- Reset dimensions of CompactListList
        //- \em without preserving existing content
        void resize_nocopy(const labelUList& listSizes);

        //- Alter local addressing size for given row, does not change content
        void setLocalSize(const label rowi, const label len);


        //- Swap contents
        void swap(CompactListList<T>& other);

        //- Transfer contents into this and annul the argument
        void transfer(CompactListList<T>& list);


    // Global addressing queries

        //- From local index on rowi to global (flat) indexing
        //- into packed values
        inline label toGlobal(const label rowi, const label i) const;

        //- From global to local index on rowi
        inline label toLocal(const label rowi, const label i) const;

        //- Find row where global index comes from. Binary search.
        //  \return -1 if out-of-bounds
        inline label findRow(const label i) const;

        //- Which row does global index come from? Binary search.
        //  FatalError if out-of-bounds
        inline label whichRow(const label i) const;


    // Pack / Unpack

        //- Unpack sub-list copies in the range defined by \p pos and \p len
        //- with bounding behaviour like List::slice() by copy constructing
        //- begin at the destination iterator \p d_iter.
        //
        // \returns Output iterator to the element in the destination range,
        //  one past the last element copied.
        template<class SubListType, class OutputIter>
        OutputIter copy_unpack
        (
            //! [out] The output destination
            OutputIter d_iter,
            //! The start of sub-region to copy (no-op if -ve or out-of-range)
            const label pos = 0,
            //! The length of sub-region to copy (-ve = until the end)
            label len = -1
        ) const;

        //- Unpack sub-list copies in the specified range.
        //
        // \returns Output iterator to the element in the destination range,
        //  one past the last element copied.
        template<class SubListType, class OutputIter>
        OutputIter copy_unpack
        (
            //! [out] The output destination
            OutputIter d_iter,
            //! The sub-region to copy
            const labelRange& range
        ) const;

        //- Unpack sub-list copies for the specified indices
        //
        // \returns Output iterator to the element in the destination range,
        //  one past the last element copied.
        template<class SubListType, class OutputIter>
        OutputIter copy_unpack
        (
            //! [out] The output destination
            OutputIter d_iter,
            //! The sub-regions to copy
            const labelUList& indices
        ) const;

        //- Return non-compact list of lists
        template<class SubListType = List<T>>
        List<SubListType> unpack() const;

        //- Return non-compact list of lists for the range of sub-lists
        template<class SubListType = List<T>>
        List<SubListType> unpack(const labelRange& range) const;

        //- Return non-compact list of lists for specified indices
        template<class SubListType = List<T>>
        List<SubListType> unpack(const labelUList& indices) const;


    // Assignment

        //- Copy assignment
        inline void operator=(const CompactListList<T>& list);

        //- Move assignment
        inline void operator=(CompactListList<T>&& list);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const Foam::zero);


    // Row-based access

        //- Return const access to sub-list (no subscript checking)
        inline const SubList<T> localList(const label i) const;

        //- Return non-const access to sub-list (no subscript checking)
        inline SubList<T> localList(const label i);

        //- Return const access to sub-list (no subscript checking)
        inline const SubList<T> operator[](const label i) const;

        //- Return non-const access to sub-list (no subscript checking)
        inline SubList<T> operator[](const label i);


    // Element access

        //- Return subscript-checked element.
        inline T& operator()(const label i, const label j);

        //- Return const subscript-checked element.
        inline const T& operator()(const label i, const label j) const;


    // Iteration (FUTURE)


    // Reading/writing

        //- Read CompactListList as offsets/values pair from Istream,
        //- discards current list contents
        Istream& readList(Istream& is);

        //- Write CompactListList as offsets/values pair
        Ostream& writeList(Ostream& os, const label shortLen=0) const;

        //- Write CompactListList as a formatted matrix of values (ASCII)
        Ostream& writeMatrix(Ostream& os, const label shortLen=0) const;


    // IO Operators

        //- Read CompactListList offsets/values pair from Istream,
        //- discarding existing contents
        friend Istream& operator>> <T>
        (
            Istream&,
            CompactListList<T>&
        );

        //- Write CompactListList as offsets/values pair
        friend Ostream& operator<< <T>
        (
            Ostream&,
            const CompactListList<T>&
        );


    // Housekeeping

        //- Const access to the packed values
        FOAM_DEPRECATED_STRICT(2022-03, "values()")
        const List<T>& m() const noexcept { return values_; }

        //- Non-const access to the packed values
        FOAM_DEPRECATED_STRICT(2022-03, "values()")
        List<T>& m() noexcept { return values_; }

        //- Return flat index into packed values
        label index(const label rowi, const label colj) const
        {
            return this->toGlobal(rowi, colj);
        }

        //- Get column within specified row that corresponds to global index
        label whichColumn(const label rowi, const label i) const
        {
            return this->toLocal(rowi, i);
        }

        //- Redimension - same as resize()
        void setSize(const label mRows)
        {
            this->resize(mRows);
        }

        //- Redimension - same as resize()
        void setSize(const label mRows, const label nVals)
        {
            this->resize(mRows, nVals);
        }

        //- Redimension - same as resize()
        void setSize(const label mRows, const label nVals, const T& val)
        {
            this->resize(mRows, nVals, val);
        }

        //- Reset sizes - same as resize()
        void setSize(const labelUList& listSizes)
        {
            this->resize(listSizes);
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Istream& operator>>(Istream& is, CompactListList<T>& list)
{
    return list.readList(is);
}

template<class T>
Ostream& operator<<(Ostream& os, const CompactListList<T>& list)
{
    return list.writeList(os, Detail::ListPolicy::short_length<T>::value);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Does not need std::swap or Foam::Swap() specialization
// since CompactListList is MoveConstructible and MoveAssignable


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CompactListListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactListList.C"
    #include "CompactListListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
