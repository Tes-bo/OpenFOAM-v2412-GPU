/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PtrList

Description
    A list of pointers to objects of type \<T\>, with allocation/deallocation
    management of the pointers.
    The operator[] returns a reference to the object, not the pointer.

See Also
    Foam::UPtrList
    Foam::PtrDynList

SourceFiles
    PtrListI.H
    PtrList.C
    PtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrList_H
#define Foam_PtrList_H

#include "UPtrList.H"
#include "SLPtrListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class autoPtr;
template<class T> class refPtr;
template<class T> class tmp;
template<class T> class PtrList;
template<class T> Istream& operator>>(Istream& is, PtrList<T>& list);

/*---------------------------------------------------------------------------*\
                           Class PtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrList
:
    public UPtrList<T>
{
protected:

    // Protected Member Functions

        //- Read from Istream using Istream constructor class
        template<class INew>
        void readIstream(Istream& is, const INew& inew);


public:

    // Constructors

        //- Default construct
        inline constexpr PtrList() noexcept;

        //- Construct with specified size, each element initialized to nullptr
        inline explicit PtrList(const label len);

        //- Copy construct using 'clone()' method on each element
        inline PtrList(const PtrList<T>& list);

        //- Move construct
        inline PtrList(PtrList<T>&& list) noexcept;

        //- Take ownership of pointers in the list, set old pointers to null.
        inline explicit PtrList(UList<T*>& list);

        //- Copy construct using 'clone()' method on each element
        template<class CloneArg>
        inline PtrList(const PtrList<T>& list, const CloneArg& cloneArgs);

        //- Construct as copy or re-use as specified
        PtrList(PtrList<T>& list, bool reuse);

        //- Copy construct using 'clone()' on each element of SLPtrList\<T\>
        explicit PtrList(const SLPtrList<T>& list);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream& is);


    //- Destructor. Frees all pointers
    ~PtrList();


    // Member Functions

        //- Make a copy by cloning each of the list elements.
        template<class... Args>
        PtrList<T> clone(Args&&... args) const;


    // Access

        //- Return const pointer to element (can be nullptr),
        //- or nullptr for out-of-range access (ie, \em with bounds checking).
        //  The return value can be tested as a bool.
        const T* set(const label i) const { return UPtrList<T>::set(i); }


    // Edit

        //- Clear the PtrList. Delete allocated entries and set size to zero.
        inline void clear();

        //- Free memory and nullify all entries. Does not change the list size.
        inline void free();

        //- Adjust size of PtrList.
        //  New entries are initialized to nullptr, removed entries are deleted
        void resize(const label newLen);

        //- Set the addressed list to the given size,
        //- deleting all existing entries.
        //- Afterwards the list contains all \c nullptr entries.
        inline void resize_null(const label newLen);

        //- Construct and append an element to the end of the list,
        //- return reference to the new list element
        template<class... Args>
        inline T& emplace_back(Args&&... args);

        //- Append an element to the end of the list
        inline void push_back(T* ptr);

        //- Move append an element to the end of the list
        inline void push_back(std::unique_ptr<T>&& ptr);

        //- Move append an element to the end of the list
        inline void push_back(autoPtr<T>&& ptr);

        //- Move or clone append a refPtr to the end of the list
        inline void push_back(const refPtr<T>& ptr);

        //- Move or clone append a tmp to the end of the list
        inline void push_back(const tmp<T>& ptr);

        //- Move append another list to the end of this list.
        inline void push_back(PtrList<T>&& other);

        //- Construct and set a new element at given position,
        //- (discard old element at that location).
        //  \param i - the location to set
        //  \param args arguments to forward to the constructor of the element
        //  \return reference to the new list element.
        template<class... Args>
        inline T& emplace_set(const label i, Args&&... args);

        //- Same as emplace_set()
        template<class... Args>
        inline T& emplace(const label i, Args&&... args);

        //- Like emplace_set() but will not overwrite an occupied (non-null)
        //- location.
        //  \param i - the location to set (unless already defined)
        //  \param args arguments to forward to the constructor of the element
        //  \return reference to the existing or the new list element.
        template<class... Args>
        inline T& try_emplace(const label i, Args&&... args);

        //- Set element to given pointer and return old element (can be null)
        //  No-op if the new pointer value is identical to the current content.
        inline autoPtr<T> set(const label i, T* ptr);

        //- Set element to given unique_ptr and return old element
        inline autoPtr<T> set(const label i, std::unique_ptr<T>&& ptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, autoPtr<T>&& ptr);

        //- Set element to given refPtr and return old element
        inline autoPtr<T> set(const label i, const refPtr<T>& ptr);

        //- Set element to given tmp and return old element
        inline autoPtr<T> set(const label i, const tmp<T>& ptr);

        //- Release ownership of the pointer at the given position.
        //  Out of bounds addressing is a no-op and returns nullptr.
        inline autoPtr<T> release(const label i);

        //- Transfer into this list and annul the argument list
        inline void transfer(PtrList<T>& list);


    // Member Operators

        //- Copy assignment.
        //  For existing list entries, values are copied from the list.
        //  For new list entries, pointers are cloned from the list.
        void operator=(const PtrList<T>& list);

        //- Move assignment
        inline void operator=(PtrList<T>&& list);


    // IOstream Operators

        //- Read from Istream, discarding contents of existing list
        friend Istream& operator>> <T>(Istream& is, PtrList<T>& list);


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
        void setSize(const label newLen) { this->resize(newLen); }

        //- Move append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(autoPtr<T>& ptr) { this->push_back(ptr.release()); }

        //- Append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T* ptr) { this->push_back(ptr); }

        //- Move append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(std::unique_ptr<T>&& ptr)
        {
            this->push_back(ptr.release());
        }

        //- Move append an element to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(autoPtr<T>&& ptr) { this->push_back(ptr.release()); }

        //- Move or clone append a tmp to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const refPtr<T>& ptr) { this->push_back(ptr); }

        //- Move or clone append a tmp to the end of the list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const tmp<T>& ptr) { this->push_back(ptr); }

        //- Move append another list to the end of this list.
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(PtrList<T>&& other) { this->push_back(std::move(other)); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrList.C"
    #include "PtrListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
