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
    Foam::PtrDynList

Description
    A dynamically resizable PtrList with allocation management.

See Also
    Foam::UPtrList
    Foam::PtrList

SourceFiles
    PtrDynListI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrDynList_H
#define Foam_PtrDynList_H

#include "PtrList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T, int SizeMin> class PtrDynList;

/*---------------------------------------------------------------------------*\
                         Class PtrDynList Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin=64>
class PtrDynList
:
    public PtrList<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private Data

        //- The capacity (allocated size) of the list.
        label capacity_;


public:

    // Constructors

        //- Default construct
        inline constexpr PtrDynList() noexcept;

        //- Construct with given capacity.
        inline explicit PtrDynList(const label len);

        //- Copy construct using 'clone()' method on each element
        inline PtrDynList(const PtrDynList<T, SizeMin>& list);

        //- Move construct
        inline PtrDynList(PtrDynList<T, SizeMin>&& list);

        //- Move construct with different sizing parameters
        template<int AnySizeMin>
        inline PtrDynList(PtrDynList<T, AnySizeMin>&& list);

        //- Move construct from PtrList
        inline PtrDynList(PtrList<T>&& list) noexcept;

        //- Take ownership of pointers in the list, set old pointers to null.
        inline explicit PtrDynList(UList<T*>& list);


    //- Destructor
    ~PtrDynList() = default;


    // Member Functions

    // Sizing

        //- Size of the underlying storage.
        label capacity() const noexcept { return capacity_; }

        //- Reserve allocation space for at least this size.
        inline void reserve(const label len);

        //- Alter the addressed list size.
        inline void resize(const label newLen);

        //- Set the addressed list to the given size,
        //- deleting all existing entries.
        //- Afterwards the list contains all \c nullptr entries.
        inline void resize_null(const label newLen);

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Shrink the allocated space to the number of elements used.
        inline void shrink_to_fit();

        //- Shrink the internal bookkeeping of the allocated space to the
        //- number of addressed elements without affecting allocation.
        //  \note when empty() it will delete any allocated memory.
        inline void shrink_unsafe();

        //- Calls shrink_to_fit()
        void shrink() { shrink_to_fit(); }


    // Edit

        //- Squeeze out intermediate nullptr entries in the list of pointers
        //- and adjust the addressable size accordingly.
        //  \return the number of non-null entries
        inline label squeezeNull();

        //- Swap with plain PtrList content. Implies shrink_to_fit().
        inline void swap(PtrList<T>& list);

        //- Swap content, independent of sizing parameter
        template<int AnySizeMin>
        inline void swap(PtrDynList<T, AnySizeMin>& other) noexcept;

        //- Transfer contents of the argument PtrList into this.
        inline void transfer(PtrList<T>& list);

        //- Transfer contents of any sized PtrDynList into this.
        template<int AnySizeMin>
        inline void transfer(PtrDynList<T, AnySizeMin>& list);

        //- Construct an element at the end of the list,
        //- return reference to the new list element
        template<class... Args>
        inline T& emplace_back(Args&&... args);

        //- Append an element to the end of the list
        inline void push_back(T* ptr);

        //- Move append an element to the end of the list
        inline void push_back(std::unique_ptr<T>&& ptr);

        //- Move append an element to the end of the list
        inline void push_back(autoPtr<T>&& ptr);

        //- Move or clone append a tmp to the end of the list
        inline void push_back(const refPtr<T>& ptr);

        //- Move or clone append a tmp to the end of the list
        inline void push_back(const tmp<T>& ptr);

        //- Move append another list to the end of this list.
        inline void push_back(PtrList<T>&& other);

        //- Move append another list to the end of this list.
        template<int AnySizeMin>
        inline void push_back(PtrDynList<T, AnySizeMin>&& other);

        //- Reduce size by 1 or more elements. Can be called on an empty list.
        inline void pop_back(label n = 1);

        //- Construct and set a new element at given position,
        //- (discard old element at that location).
        //- Auto-sizes list as required.
        //  \param i - the location to set
        //  \param args arguments to forward to the constructor of the element
        //  \return reference to the new list element.
        template<class... Args>
        inline T& emplace_set(const label i, Args&&... args);

        //- Same as emplace_set()
        template<class... Args>
        inline T& emplace(const label i, Args&&... args);

        //- Like emplace_set() but will not overwrite an occupied location.
        //  \param i - the location to set (unless already defined)
        //  \param args arguments to forward to the constructor of the element
        //  \return reference to the existing or the new list element.
        template<class... Args>
        inline T& try_emplace(const label i, Args&&... args);

        //- Set element to given pointer and return old element (can be null).
        //- Auto-sizes list as required.
        inline autoPtr<T> set(const label i, T* ptr);

        //- Set element to given pointer and return old element
        //- Auto-sizes list as required.
        inline autoPtr<T> set(const label i, std::unique_ptr<T>&& ptr);

        //- Set element to given autoPtr and return old element
        //- Auto-sizes list as required.
        inline autoPtr<T> set(const label i, autoPtr<T>&& ptr);

        //- Set element to given refPtr and return old element
        //- Auto-sizes list as required.
        inline autoPtr<T> set(const label i, const refPtr<T>& ptr);

        //- Set element to given tmp and return old element
        //- Auto-sizes list as required.
        inline autoPtr<T> set(const label i, const tmp<T>& ptr);

        //- Reorder elements. Reordering must be unique (ie, shuffle).
        inline void reorder(const labelUList& oldToNew);


    // Member Operators

        //- Copy (clone) assignment
        inline void operator=(const PtrList<T>& list);

        //- Copy (clone) assignment
        inline void operator=(const PtrDynList<T, SizeMin>& list);

        //- Copy (clone) assignment with different sizing parameters
        template<int AnySizeMin>
        inline void operator=(const PtrDynList<T, AnySizeMin>& list);

        //- Move assignment
        inline void operator=(PtrList<T>&& list);

        //- Move assignment
        inline void operator=(PtrDynList<T, SizeMin>&& list);

        //- Move assignment with different sizing parameters
        template<int AnySizeMin>
        inline void operator=(PtrDynList<T, AnySizeMin>&& list);


    // Housekeeping

        //- Disallow push_back with autoPtr without std::move
        void push_back(autoPtr<T>& ptr) = delete;

        //- Set element to given autoPtr and return old element
        //FOAM_DEPRECATED_FOR(2022-10, "set(autoPtr&&))")
        autoPtr<T> set(const label i, autoPtr<T>& ptr)
        {
            return this->set(i, ptr.release());
        }

        //- Same as resize()
        void setSize(const label n) { this->resize(n); }

        //- Move append an element to the end of the list
        void append(autoPtr<T>& ptr) { this->push_back(std::move(ptr)); }

        //- Append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T* ptr) { this->push_back(ptr); }

        //- Move append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(std::unique_ptr<T>&& ptr)
        {
            this->push_back(std::move(ptr));
        }

        //- Move append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(autoPtr<T>&& ptr) { this->push_back(std::move(ptr)); }

        //- Move or clone append a tmp to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const refPtr<T>& ptr) { this->push_back(ptr); }

        //- Move or clone append a tmp to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const tmp<T>& ptr) { this->push_back(ptr); }

        //- Move append another list to the end of this list.
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(PtrList<T>&& other) { this->push_back(std::move(other)); }

        //- Move append another list to the end of this list.
        template<int AnySizeMin>
        void append(PtrDynList<T, AnySizeMin>&& other)
        {
            this->push_back(std::move(other));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrDynListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
