/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::Detail::PtrListDetail

Description
    A rudimentary list of pointers used for PtrList, UPtrList, etc.
    This class is considered implementation detail and should not normally
    be used other than by OpenFOAM container classes.

    It stores a list of pointers, but makes leaves memory management
    to the caller or sub-class.
    The free() method can be used explicitly as required.

SourceFiles
    PtrListDetail.C
    PtrListDetailI.H
    PtrListDetailIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrListDetail_H
#define Foam_PtrListDetail_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                    Class Detail::PtrListDetail Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrListDetail
:
    public List<T*>
{
public:

    // Constructors

        //- Default construct
        inline constexpr PtrListDetail() noexcept;

        //- Construct with specified size, each element initialized to nullptr
        inline explicit PtrListDetail(const label len);

        //- Copy a list of pointers.
        //  The caller is responsible for memory management.
        inline explicit PtrListDetail(const UList<T*>& list);

        //- Copy construct (shallow copies addresses)
        inline PtrListDetail(const PtrListDetail<T>& list);

        //- Move construct
        inline PtrListDetail(PtrListDetail<T>&& list) noexcept;

        //- Copy or move (reuse) construct as specified
        inline PtrListDetail(PtrListDetail<T>& list, bool reuse);


    // Member Functions

        //- Return const pointer to element or nullptr for out-of-range access.
        inline const T* get(const label i) const;

        //- Return pointer to element or nullptr for out-of-range access.
        inline T* get(const label i);

        //- The number of non-nullptr entries in the list
        label count_nonnull() const noexcept;

        //- FatalError if any null exists in the list
        inline void checkNonNull() const;

        //- Locate the first entry that is non-null.
        //  \return the location or -1 if there are no bits set.
        //
        //  \note Method name as per bitSet
        label find_first() const;

        //- Locate the first entry that is null, -1 if there are none (or empty list)
        //  \return the location or -1 if the list is empty or all entries
        //  are non-null.
        //
        //  \note Method name as per bitSet, symmetry with find_first()
        label find_first_not() const;

        //- Locate the next non-null entry, starting one beyond the specified
        //- position
        //  \return the location or -1 if remaining entries are null.
        //
        //  \note Method name as per bitSet
        label find_next(label pos) const;

        //- Locate the next null entry, starting one beyond the specified
        //- position
        //  \return the location or -1 if remaining entries are non-null.
        //
        //  \note Method name as per bitSet
        label find_next_not(label pos) const;

        //- Delete allocated entries and reassign to nullptr.
        //- Does not affect the list size.
        void free();

        //- Make a copy by cloning each of the list pointers.
        template<class... Args>
        PtrListDetail<T> clone(Args&&... args) const;

        //- Reset size of list.
        //  New entries are initialized to nullptr.
        inline void resize(const label newLen);

        //- Set the list to the given size and set all entries to nullptr.
        inline void resize_null(const label newLen);

        //- Set addressed size to be inconsistent with allocated storage.
        //  Use with care
        inline void setAddressableSize(const label n) noexcept;

        //- Write output, optionally silently trimming nullptrs
        Ostream& write(Ostream& os, const bool trimNull=false) const;


    // Member Operators

        //- Copy assignment (shallow copies addresses)
        inline void operator=(const PtrListDetail<T>& list);

        //- Move assignment
        inline void operator=(PtrListDetail<T>&& list);

        //- Assign all entries to nullptr (without deleting)
        inline void operator=(std::nullptr_t);


    // Housekeeping

        // Just use resize().
        void setSize(const label) = delete;
        void setSize(const label, const T&) = delete;
        void setSize(const label, const T*) = delete;

        // Too fragile or dangerous
        void resize_nocopy(const label) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListDetailI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrListDetail.C"
    #include "PtrListDetailIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
