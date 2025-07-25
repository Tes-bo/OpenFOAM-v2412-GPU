/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::labelRanges

Description
    A list of labelRange with constrained list capabilities.

SourceFiles
    labelRangesI.H
    labelRanges.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_labelRanges_H
#define Foam_labelRanges_H

#include "labelRange.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class labelRanges;
Istream& operator>>(Istream& is, labelRanges& ranges);
Ostream& operator<<(Ostream& is, const labelRanges& ranges);

/*---------------------------------------------------------------------------*\
                         Class labelRanges Declaration
\*---------------------------------------------------------------------------*/

class labelRanges
{
    // Private Data

        //- The list of ranges
        DynamicList<labelRange> ranges_;


    // Private Member Functions

        //- Insert range before specified insertion index, by copying up
        void insertBefore(const label insert, const labelRange& range);

        //- Purge empty ranges, by copying down
        void purgeEmpty();

public:

    // STL type definitions

        //- The value type the list contains
        typedef labelRange value_type;

        //- Input iterator with const access
        class const_iterator;


    // Generated Methods

        //- Default construct
        labelRanges() = default;

        //- Default copy construct
        labelRanges(const labelRanges&) = default;

        //- Default move construct
        labelRanges(labelRanges&&) = default;

        //- Default copy assignment
        labelRanges& operator=(const labelRanges&) = default;

        //- Default move assignment
        labelRanges& operator=(labelRanges&&) = default;


    // Constructors

        //- Construct an empty list with given initial capacity
        inline explicit labelRanges(const label initialCapacity);

        //- Copy construct from list of ranges
        inline explicit labelRanges(const UList<labelRange>& list);

        //- Move construct from list of ranges
        inline labelRanges(List<labelRange>&& list);

        //- Move construct from list of ranges
        template<int AnySizeMin>
        inline labelRanges(DynamicList<labelRange, AnySizeMin>&& list);

        //- Construct from Istream.
        explicit labelRanges(Istream& is);


    // Member Functions

        //- The list of ranges
        const UList<labelRange>& ranges() const noexcept { return ranges_; }

        //- Clear the addressable list of ranges
        void clear() noexcept { ranges_.clear(); }

        //- True if list of ranges is empty
        bool empty() const noexcept { return ranges_.empty(); }

        //- Reserve space for at least this size
        void reserve(const label len) { ranges_.reserve(len); }

        //- The linear size (sum of all the element sizes)
        inline label totalSize() const noexcept;

        //- True if the value is contained within any of the sub-ranges
        inline bool contains(const label value) const noexcept;

        //- Add the range to the list
        bool add(const labelRange& range);

        //- Remove the range from the list
        bool remove(const labelRange& range);

        //- Construct a range element at the end of the list,
        //- return reference to the new element.
        template<class... Args>
        inline labelRange& emplace_back(Args&&... args);

        //- Inplace sort of the range elements
        inline void sort();

        //- Return flattened list of all range labels
        List<label> labels() const;


    // Member Operators

       //- Return the value at linear index 'i', -1 for out-of-range
       label operator[](const label i) const;


    // STL iterator

        //- Forward input iterator with const access
        class const_iterator
        {
            //- The list of ranges for which this is an iterator
            const UList<labelRange>* list_;

            //- The index into the list
            label index_;

            //- Index of current element at list-index
            label subIndex_;

        public:

        // Constructors

            //- Construct from range list at given index (and sub-index)
            inline explicit constexpr const_iterator
            (
                const UList<labelRange>* list,
                const label idx = 0,
                const label subIdx = 0
            ) noexcept;


        // Member Operators

            //- Return the current label
            inline label operator*() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline constexpr bool operator==(const const_iterator& iter)
            const noexcept;

            inline constexpr bool operator!=(const const_iterator& iter)
            const noexcept;
        };


        //- A const_iterator set to the beginning of the list
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator cend() const noexcept;

        //- A const_iterator set to the beginning of the list
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to beyond the end of the list
        inline const const_iterator end() const noexcept;


        //- Return const_iterator at linear offset i from begin,
        //- clamped to [0,size] range
        inline const_iterator cbegin(const label i) const;

        //- Return const_iterator at linear offset i from begin,
        //- clamped to [0,size] range
        inline const_iterator begin(const label i) const;


    // Reading/writing

        //- Read List of labelRange from Istream, discarding contents
        Istream& readList(Istream& is);

        //- Write List of labelRange, with line-breaks in ASCII
        //- when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;


    // IOstream Operators

        //- Use the readList() method to read contents from Istream.
        friend Istream& operator>>(Istream& is, labelRanges& list);

        //- Write to Ostream. Uses the writeList() method
        friend Ostream& operator<<(Ostream& os, const labelRanges& list);


    // Housekeeping

        //- Same as contains()
        bool found(const label value) const { return contains(value); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
