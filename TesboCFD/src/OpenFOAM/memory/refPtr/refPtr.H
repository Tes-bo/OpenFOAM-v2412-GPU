/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::refPtr

Description
    A class for managing references or pointers (no reference counting)

SourceFiles
    refPtrI.H

See also
    Foam::autoPtr
    Foam::tmp

\*---------------------------------------------------------------------------*/

#ifndef Foam_refPtr_H
#define Foam_refPtr_H

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class refPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class refPtr
{
    // Private Data

        //- The stored reference type
        enum refType
        {
            PTR,        //!< A managed pointer (not ref-counted)
            REF_Types,  //!< Magic value (ptr: less-than, ref: greater than)
            CREF,       //!< A const reference to an object
            REF         //!< A non-const reference to an object
        };

        //- The managed pointer or address of the object (reference)
        mutable T* ptr_;

        //- The type (managed pointer | object reference)
        mutable refType type_;


public:

    // STL type definitions

        //- Type of object being managed or referenced
        typedef T element_type;

        //- Pointer to type of object being managed or referenced
        typedef T* pointer;


    //- Null reference counter class
    typedef Foam::refCount::zero refCount;


    // Constructors

        //- Construct with no managed pointer.
        inline constexpr refPtr() noexcept;

        //- Implicit construct from literal nullptr: no managed pointer
        inline constexpr refPtr(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit constexpr refPtr(T* p) noexcept;

        //- Implicit construct for a const reference to an object.
        inline constexpr refPtr(const T& obj) noexcept;

        //- Move construct, transferring ownership.
        inline refPtr(refPtr<T>&& rhs) noexcept;

        //- Copy construct (shallow copy)
        inline refPtr(const refPtr<T>& rhs);

        //- Copy/move construct. Optionally reusing pointer.
        inline refPtr(const refPtr<T>& rhs, bool reuse);

        //- Move construct from unique_ptr, transferring ownership.
        inline explicit refPtr(std::unique_ptr<T>&& rhs) noexcept;

        //- No copy construct from autoPtr,
        //- also avoids implicit cast to object or pointer
        refPtr(const autoPtr<T>&) = delete;

        //- Move construct from autoPtr, transferring ownership.
        inline explicit refPtr(autoPtr<T>&& rhs) noexcept;

        //- Reference tmp contents or transfer ownership if requested/possible
        inline refPtr(const tmp<T>& rhs, bool reuse);

        //- Reference the tmp contents
        inline explicit refPtr(const tmp<T>& rhs);

        //- Move construct from tmp, transfer ownership if possible.
        inline explicit refPtr(tmp<T>&& rhs);


    //- Destructor: deletes managed pointer
    inline ~refPtr() noexcept;


    // Factory Methods

        //- Construct refPtr with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_shared, but the overload for
        //      array types is not disabled.
        template<class... Args>
        static refPtr<T> New(Args&&... args)
        {
            return refPtr<T>(new T(std::forward<Args>(args)...));
        }

        //- Construct refPtr from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types.
        //  Future check std::enable_if + std::is_convertible ?
        template<class U, class... Args>
        static refPtr<T> NewFrom(Args&&... args)
        {
            return refPtr<T>(new U(std::forward<Args>(args)...));
        }


    // Static Member Functions

        //- The type-name, constructed from type-name of T
        inline static word typeName();


    // Query

        //- True if pointer/reference is non-null.
        bool good() const noexcept { return bool(ptr_); }

        //- If the stored/referenced content is const
        inline bool is_const() const noexcept;

        //- True if this is a managed pointer (not a reference)
        inline bool is_pointer() const noexcept;

        //- True if this is a reference (not a pointer)
        inline bool is_reference() const noexcept;

        //- True if this is a non-null managed pointer
        inline bool movable() const noexcept;


    // Access

        //- Return pointer without nullptr checking.
        T* get() noexcept { return ptr_; }

        //- Return const pointer without nullptr checking.
        const T* get() const noexcept { return ptr_; }

        //- Return const reference to the object or to the contents
        //- of a (non-null) managed pointer.
        //  Fatal for a null managed pointer
        inline const T& cref() const;

        //- Return non-const reference to the contents of a non-null
        //- managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline T& ref() const;

        //- Return non-const reference to the object or to the contents
        //- of a (non-null) managed pointer, with an additional const_cast.
        //  Fatal for a null managed pointer.
        T& constCast() const { return const_cast<T&>(cref()); }

        //- Return a shallow copy as a wrapped reference, preserving the
        //- const/non-const status.
        inline refPtr<T> shallowClone() const noexcept;


    // Edit

        //- Release ownership and return the pointer.
        //- A no-op for reference objects (returns nullptr).
        //  \remark Method naming consistent with std::unique_ptr
        inline T* release() noexcept;

        //- Return managed pointer for reuse, or clone() the object reference.
        //  \remark Method naming consistent with Foam::tmp
        inline T* ptr() const;

        //- If object pointer points to valid object:
        //- delete object and set pointer to nullptr
        inline void clear() const noexcept;

        //- Delete managed pointer and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Clear existing and transfer ownership.
        inline void reset(refPtr<T>&& other) noexcept;

        //- No reset from autoPtr reference (potentially confusing)
        void reset(const autoPtr<T>&) = delete;

        //- Clear existing and transfer ownership from autoPtr.
        void reset(autoPtr<T>&& other) noexcept { reset(other.release()); }

        //- Clear existing and transfer ownership from unique_ptr
        void reset(std::unique_ptr<T>&& other) { reset(other.release()); }

        //- Reference tmp contents or transfer pointer ownership if possible
        inline void reset(tmp<T>& rhs, bool reuse);

        //- Reset with emplace construction.
        //- Return reference to the new content.
        template<class... Args>
        inline T& emplace(Args&&... args);

        //- Clear existing and set (const) reference from other
        inline void cref(const refPtr<T>& other) noexcept;

        //- Clear existing and set (const) reference
        inline void cref(const T& obj) noexcept;

        //- Clear existing and set (const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void cref(const T* p) noexcept;

        //- Avoid inadvertent casting (to object or pointer)
        void cref(const autoPtr<T>&) = delete;

        //- Avoid inadvertent casting (to object)
        void cref(const tmp<T>&) = delete;

        //- Clear existing and set (non-const) reference
        inline void ref(T& obj) noexcept;

        //- Clear existing and set (non-const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void ref(T* p) noexcept;

        //- Avoid inadvertent casting (to object or pointer)
        void ref(const autoPtr<T>&) = delete;

        //- Avoid inadvertent casting (object)
        void ref(const tmp<T>&) = delete;

        //- Swaps the managed object with other.
        inline void swap(refPtr<T>& other) noexcept;


    // Member Operators

        //- Return const reference to the object.
        //  Fatal if nothing is managed
        inline const T& operator*() const;

        //- Return reference to the managed object.
        //  Fatal if nothing is managed or if the object is const.
        inline T& operator*();

        //- Dereferences (const) pointer to the managed object.
        //  Fatal if nothing is managed.
        inline const T* operator->() const;

        //- Dereferences (non-const) pointer to the managed object.
        //  Fatal if nothing is managed or if the object is const.
        inline T* operator->();

        //- Return const reference to the object - same as cref() method.
        const T& operator()() const { return cref(); }


    // Casting

        //- True if pointer/reference is non-null. Same as good()
        explicit operator bool() const noexcept { return bool(ptr_); }

        //- Cast to underlying data type, using the cref() method.
        operator const T&() const { return cref(); }


    // Assignment

        //- Transfer ownership of managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const refPtr<T>& other);

        //- Clear existing and transfer ownership.
        void operator=(refPtr<T>&& other) noexcept { reset(std::move(other)); }

        //- No copy assignment from plain pointer (uncontrolled access)
        void operator=(T* p) = delete;

        //- Reset via assignment from literal nullptr
        void operator=(std::nullptr_t) noexcept { reset(nullptr); }

        //- No copy assignment from autoPtr (can have unintended behaviour)
        void operator=(const autoPtr<T>&) = delete;

        //- Transfer ownership by move assignment from autoPtr.
        void operator=(autoPtr<T>&& other) noexcept { reset(other.release()); }

        //- Transfer ownership by move assignment from unique_ptr
        void operator=(std::unique_ptr<T>&& other) { reset(other.release()); }

        //- No copy assignment from tmp
        void operator=(const tmp<T>&) = delete;

        //- Move construct, transferring pointer ownership if possible.
        inline void operator=(tmp<T>&& other);


        //- Conversion to tmp, releases pointer or shallow-copies reference
        inline operator tmp<T>();


    // Housekeeping

        //- Identical to good(), or bool operator
        bool valid() const noexcept { return bool(ptr_); }

        //- Deprecated(2020-07) True if a null managed pointer
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Does not need std::swap or Foam::Swap() specialization
// since refPtr is MoveConstructible and MoveAssignable


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
