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
    Foam::tmp

Description
    A class for managing temporary objects.

    This is a combination of std::shared_ptr (with intrusive ref-counting)
    and a shared_ptr without ref-counting and null deleter.
    This allows the tmp to double as a pointer management and an indirect
    pointer to externally allocated objects.
    In contrast to std::shared_ptr, only a limited number of tmp items
    will ever share a pointer.

SourceFiles
    tmpI.H

See also
    Foam::autoPtr
    Foam::refPtr
    Foam::refCount

\*---------------------------------------------------------------------------*/

#ifndef Foam_tmp_H
#define Foam_tmp_H

#include "autoPtr.H"
#include "refCount.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class refPtr;

/*---------------------------------------------------------------------------*\
                             Class tmp Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmp
{
    // Private Data

        //- The stored reference type
        enum refType
        {
            PTR,        //!< A managed pointer (ref-counted)
            CACHE_PTR,  //!< A managed pointer (ref-counted), but immovable
            REF_Types,  //!< Magic value (ptr: less-than, ref: greater than)
            CREF,       //!< A const reference to an object
            REF         //!< A non-const reference to an object
        };

        //- The managed pointer or address of the object (reference)
        mutable T* ptr_;

        //- The type (managed pointer | object reference)
        mutable refType type_;


    // Private Member Functions

        //- Fatal if ref-count of a managed pointer is oversubscribed
        inline void checkUseCount() const;


public:

    // STL type definitions

        //- Type of object being managed or referenced
        typedef T element_type;

        //- Pointer to type of object being managed or referenced
        typedef T* pointer;


    //- Reference counter class
    typedef Foam::refCount refCount;


    // Constructors

        //- Construct with no managed pointer.
        inline constexpr tmp() noexcept;

        //- Implicit construct from literal nullptr: no managed pointer
        inline constexpr tmp(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit tmp(T* p);

        //- Implicit construct for a const reference to an object.
        inline constexpr tmp(const T& obj) noexcept;

        //- Move construct, transferring ownership.
        //  Does not affect ref-count
        inline tmp(tmp<T>&& rhs) noexcept;

        //- Move construct, transferring ownership.
        //  Does not affect ref-count
        //  \note Non-standard definition - should be non-const
        inline tmp(const tmp<T>&& rhs) noexcept;

        //- Copy construct, incrementing ref-count of managed pointer.
        //  \note Non-standard definition - should be non-const
        inline tmp(const tmp<T>& rhs);

        //- Copy/move construct. Optionally reusing ref-counted pointer.
        inline tmp(const tmp<T>& rhs, bool reuse);

        //- No copy construct from autoPtr,
        //- also avoids implicit cast to object or pointer
        tmp(const autoPtr<T>&) = delete;

        //- No copy construct from refPtr,
        //- also avoids implicit cast to object
        tmp(const refPtr<T>&) = delete;

        //- Move construct from autoPtr, transferring ownership.
        inline explicit tmp(autoPtr<T>&& rhs) noexcept;


    //- Destructor: deletes managed pointer when the ref-count is 0
    inline ~tmp() noexcept;


    // Factory Methods

        //- Construct tmp with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_shared, but the overload for
        //      array types is not disabled.
        template<class... Args>
        static tmp<T> New(Args&&... args)
        {
            return tmp<T>(new T(std::forward<Args>(args)...));
        }

        //- Construct tmp from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types.
        //  Future check std::enable_if + std::is_convertible ?
        template<class U, class... Args>
        static tmp<T> NewFrom(Args&&... args)
        {
            return tmp<T>(new U(std::forward<Args>(args)...));
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

        //- True if this is a non-null managed pointer with a unique ref-count
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
        //  Fatal for a null pointer.
        T& constCast() const { return const_cast<T&>(cref()); }


    // Edit

        //- Return managed pointer for reuse, or clone() the object reference.
        inline T* ptr() const;

        //- If object pointer points to valid object:
        //- delete object and set pointer to nullptr
        inline void clear() const noexcept;

        //- Use specified protection against moving for the managed pointer.
        //- No-op for references.
        inline void protect(bool on) noexcept;

        //- Clear existing and transfer ownership.
        inline void reset(tmp<T>&& other) noexcept;

        //- Delete managed temporary object and set to new given pointer.
        inline void reset(T* p = nullptr) noexcept;

        //- Avoid inadvertent casting (to object or pointer)
        void reset(const autoPtr<T>&) = delete;

        //- Avoid inadvertent casting (to object)
        void reset(const refPtr<T>&) = delete;

        //- Reset with emplace construction.
        //- Return reference to the new content.
        template<class... Args>
        inline T& emplace(Args&&... args);

        //- Clear existing and set (const) reference from other
        inline void cref(const tmp<T>& other) noexcept;

        //- Clear existing and set (const) reference
        inline void cref(const T& obj) noexcept;

        //- Clear existing and set (const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void cref(const T* p) noexcept;

        //- Avoid inadvertent casting (to object or pointer)
        void cref(const autoPtr<T>&) = delete;

        //- Avoid inadvertent casting (to object)
        void cref(const refPtr<T>&) = delete;

        //- Clear existing and set to (non-const) reference
        inline void ref(T& obj) noexcept;

        //- Clear existing and set (non-const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void ref(T* p) noexcept;

        //- Avoid inadvertent casting (to object or pointer)
        void ref(const autoPtr<T>&) = delete;

        //- Avoid inadvertent casting (to object)
        void ref(const refPtr<T>&) = delete;

        //- Swaps the managed object with other.
        inline void swap(tmp<T>& other) noexcept;


    // Member Operators

        // Note: no 'operator*()' types since there are large chunks
        // of code that use tmp\<Field\> and Field interchangeable
        // Eg, \code (a * b) - and the '*' would be misinterpreted

        //- Dereferences (const) pointer to the managed object.
        //  Fatal for a null managed pointer.
        inline const T* operator->() const;

        //- Dereferences (non-const) pointer to the managed object.
        //  Fatal for a null managed pointer or if the object is const.
        inline T* operator->();

        //- Return const reference to the object - same as cref() method.
        const T& operator()() const { return cref(); }


    // Casting

        //- True if pointer/reference is non-null. Same as good()
        explicit operator bool() const noexcept { return bool(ptr_); }

        //- Cast to underlying data type, using the cref() method.
        operator const T&() const { return cref(); }


    // Assignment

        //- Transfer ownership of the managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const tmp<T>& other);

        //- Clear existing and transfer ownership.
        inline void operator=(tmp<T>&& other) noexcept;

        //- Take ownership of the pointer.
        //  Fatal for a null pointer, or when the pointer is non-unique.
        inline void operator=(T* p);

        //- Reset via assignment from literal nullptr
        void operator=(std::nullptr_t) noexcept { reset(nullptr); }


    // Housekeeping

        //- Identical to good(), or bool operator
        bool valid() const noexcept { return bool(ptr_); }

        //- Identical to is_pointer(). Prefer is_pointer() or movable().
        //
        //  \deprecated(2020-07) - prefer is_pointer() or movable()
        bool isTmp() const noexcept { return is_pointer(); }

        //- Deprecated(2020-07) True if a null managed pointer
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Does not need std::swap or Foam::Swap() specialization
// since tmp is MoveConstructible and MoveAssignable


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
