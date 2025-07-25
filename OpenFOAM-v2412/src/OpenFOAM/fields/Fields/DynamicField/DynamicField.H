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
    Foam::DynamicField

Description
    Dynamically sized Field.

SourceFiles
    DynamicFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_DynamicField_H
#define Foam_DynamicField_H

#include "Field.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T, int SizeMin> class DynamicField;

template<class T, int SizeMin>
inline Istream& operator>>(Istream&, DynamicField<T, SizeMin>&);

template<class T, int SizeMin>
inline Ostream& operator<<(Ostream&, const DynamicField<T, SizeMin>&);

/*---------------------------------------------------------------------------*\
                        Class DynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin=64>
class DynamicField
:
    public Field<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private Data

        //- The capacity (allocated size) of the underlying field.
        label capacity_;


    // Private Member Functions

        //- Copy assignment from another list
        template<class ListType>
        inline void doAssignDynList(const ListType& list);

        //- Alter the size of the underlying storage
        //  The 'nocopy' option will not attempt to recover old content
        inline void doCapacity(const bool nocopy, const label len);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        //  The 'nocopy' option will not attempt to recover old content
        inline void doReserve(const bool nocopy, const label len);

        //- Reserve allocation space for at least this size.
        //  Never shrinks the allocated size, use setCapacity() for that.
        //  The 'nocopy' option will not attempt to recover old content
        inline void doResize(const bool nocopy, const label len);


public:

    // Static Member Functions

        //- Return a null DynamicField (reference to a nullObject).
        static const DynamicField<T, SizeMin>& null() noexcept
        {
            return NullObjectRef<DynamicField<T, SizeMin>>();
        }


    // Constructors

        //- Default construct, an empty field without allocation.
        inline constexpr DynamicField() noexcept;

        //- Construct empty field with given initial capacity
        inline explicit DynamicField(const label initialCapacity);

        //- Construct given size and initial value
        inline DynamicField(const label len, const T& val);

        //- Construct given size and initial value of zero
        inline DynamicField(const label len, const Foam::zero);

        //- Copy construct
        inline DynamicField(const DynamicField<T, SizeMin>& list);

        //- Copy construct with different sizing parameters
        template<int AnySizeMin>
        inline DynamicField(const DynamicField<T, AnySizeMin>& list);

        //- Copy construct from UList. Size set to UList size.
        //  Also constructs from DynamicField with different sizing parameters.
        inline explicit DynamicField(const UList<T>& list);

        //- Copy construct from IndirectList
        template<class Addr>
        inline explicit DynamicField(const IndirectListBase<T, Addr>& list);

        //- Move construct from List contents
        inline explicit DynamicField(List<T>&& content) noexcept;

        //- Move construct from dynamic Field contents
        inline DynamicField(DynamicField<T, SizeMin>&& content) noexcept;

        //- Move construct with different sizing parameters
        template<int AnySizeMin>
        inline DynamicField(DynamicField<T, AnySizeMin>&& content) noexcept;

        //- Move construct from DynamicList
        template<int AnySizeMin>
        inline DynamicField(DynamicList<T, AnySizeMin>&& content) noexcept;

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelUList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const FieldMapper& map
        );

        //- Construct from Istream. Size set to size of list read.
        inline explicit DynamicField(Istream& is);

        //- Clone
        inline tmp<DynamicField<T, SizeMin>> clone() const;


    // Member Functions

    // Capacity

        //- Normal lower capacity limit - the SizeMin template parameter
        static constexpr label min_size() noexcept { return SizeMin; }

        //- Size of the underlying storage.
        label capacity() const noexcept { return capacity_; }

        //- Number of contiguous bytes of the underlying storage.
        //  \note Only meaningful for contiguous data
        inline std::streamsize capacity_bytes() const noexcept;


    // Sizing

        //- Alter the size of the underlying storage.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        //  Use this or reserve() in combination with push_back().
        inline void setCapacity(const label len);

        //- Alter the size of the underlying storage,
        //- \em without retaining old content.
        //  The addressed size will be truncated if needed to fit, but will
        //  remain otherwise untouched.
        inline void setCapacity_nocopy(const label len);

        //- Change the value for the list capacity directly (ADVANCED, UNSAFE)
        //- Does not perform any memory management or resizing.
        inline void setCapacity_unsafe(const label len) noexcept;

        //- Reserve allocation space for at least this size, allocating new
        //- space if required and \em retaining old content.
        //  Never shrinks the allocated size, use setCapacity() for that.
        inline void reserve(const label len);

        //- Reserve allocation space for at least this size, allocating new
        //- space if required \em without retaining old content.
        //  Never shrinks the allocated size, use setCapacity() for that.
        inline void reserve_nocopy(const label len);

        //- Alter addressable list size, allocating new space if required
        //- while \em recovering old content.
        //  If no reallocation is required, the contents remain untouched.
        //  Otherwise new entries will be uninitialized.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        inline void resize(const label len);

        //- Alter addressable size and fill new space with constant value
        inline void resize(const label len, const T& val);

        //- Alter addressable list size, allocating new space if required
        //- \em without necessarily recovering old content.
        //  If no reallocation is required, the contents remain untouched.
        //  Otherwise all entries will be uninitialized.
        inline void resize_nocopy(const label len);

        //- Alias for resize()
        void setSize(const label n) { this->resize(n); }

        //- Alias for resize()
        void setSize(const label n, const T& val) { this->resize(n, val); }

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear() noexcept;

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to the number of elements used.
        inline void shrink_to_fit();

        //- Shrink the internal bookkeeping of the allocated space to the
        //- number of addressed elements without affecting allocation.
        //  \note when empty() it will delete any allocated memory.
        inline void shrink_unsafe();

        //- Calls shrink_to_fit() and returns a reference to the DynamicField.
        inline DynamicField<T, SizeMin>& shrink();


    // Edit

        //- Swap with plain List content. Implies shrink_to_fit().
        inline void swap(List<T>& list);

        //- Swap content, independent of sizing parameter
        template<int AnySizeMin>
        inline void swap(DynamicField<T, AnySizeMin>& other);

        //- Swap content with DynamicList, independent of sizing parameter
        template<int AnySizeMin>
        inline void swap(DynamicList<T, AnySizeMin>& other);

        //- Transfer the parameter contents into this
        inline void transfer(List<T>& list);

        //- Transfer the parameter contents into this
        template<int AnySizeMin>
        inline void transfer(DynamicList<T, AnySizeMin>& list);

        //- Transfer the parameter contents into this
        template<int AnySizeMin>
        inline void transfer(DynamicField<T, AnySizeMin>& list);

        //- Construct an element at the end of the list,
        //- return reference to the new list element
        template<class... Args>
        inline T& emplace_back(Args&&... args);

        //- Copy append an element at the end of the list
        inline void push_back(const T& val);

        //- Move append an element
        inline void push_back(T&& val);

        //- Copy append another list to the end of this list
        inline void push_back(const UList<T>& list);

        //- Move append another list to the end of this list
        inline void push_back(List<T>&& list);

        //- Reduce size by 1 or more elements. Can be called on an empty list.
        inline void pop_back(label n = 1);


    // Reading/writing

        //- Read from Istream, discarding existing contents
        //  Uses a DynamicList::readList internally
        inline Istream& readList(Istream& is);


    // Assignment

        //- Assign addressed entries to the given value
        inline void operator=(const T& val);

        //- Assign addressed entries to zero
        inline void operator=(const Foam::zero);

        //- Copy assignment
        inline void operator=(const UList<T>& list);

        //- Copy assignment
        inline void operator=(const DynamicField<T, SizeMin>& list);

        //- Copy assign from IndirectList
        template<class Addr>
        inline void operator=(const IndirectListBase<T, Addr>& rhs);

        //- Move assignment
        inline void operator=(List<T>&& list);

        //- Move assignment
        template<int AnySizeMin>
        inline void operator=(DynamicList<T, AnySizeMin>&& list);

        //- Move assignment
        inline void operator=(DynamicField<T, SizeMin>&& list);

        //- Move assignment
        template<int AnySizeMin>
        inline void operator=(DynamicField<T, AnySizeMin>&& list);


    // Member Operators

        //- Return non-const access to an element, resizing list if needed
        inline T& operator()(const label i);


    // IOstream Operators

        //- Use the readList() method to read contents from Istream.
        friend Istream& operator>> <T, SizeMin>
        (
            Istream& is,
            DynamicField<T, SizeMin>& rhs
        );

        //- Write to Ostream
        friend Ostream& operator<< <T, SizeMin>
        (
            Ostream& os,
            const DynamicField<T, SizeMin>& rhs
        );


    // Housekeeping

        //- Append an element at the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const T& val) { this->push_back(val); }

        //- Move append an element
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T&& val) { this->push_back(std::move(val)); }

        //- Append a List at the end of this list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const UList<T>& list) { this->push_back(list); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
