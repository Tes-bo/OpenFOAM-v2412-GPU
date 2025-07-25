/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::autoPtr

Description
    Pointer management similar to std::unique_ptr, with some additional
    methods and type checking.

Note
    Parts of the interface now mirror std::unique_ptr, but since it pre-dates
    both C++11 and std::unique_ptr, it has some additional idiosyncrasies.
    The const-reference constructors and assignment operators
    actually use move semantics to allow their participation in
    default constructible, default assignable classes.

SourceFiles
    autoPtrI.H

See also
    Foam::refPtr

\*---------------------------------------------------------------------------*/

#ifndef Foam_autoPtr_H
#define Foam_autoPtr_H

// Transitional features/misfeatures
#define Foam_autoPtr_copyConstruct
#define Foam_autoPtr_castOperator
// #define Foam_autoPtr_nodeprecate_setMethod

#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtr
{
    //- The managed pointer
    T* ptr_;


public:

    // STL type definitions

        //- Type of object being managed
        typedef T element_type;

        //- Pointer to type of object being managed
        typedef T* pointer;


    // Constructors

        //- Construct with no managed pointer.
        constexpr autoPtr() noexcept
        :
            ptr_(nullptr)
        {}

        //- Implicit construct from literal nullptr: no managed pointer
        constexpr autoPtr(std::nullptr_t) noexcept
        :
            ptr_(nullptr)
        {}

        //- Construct, taking ownership of the pointer.
        explicit autoPtr(T* p) noexcept
        :
            ptr_(p)
        {}

        //- Move construct, transferring ownership.
        autoPtr(autoPtr<T>&& rhs) noexcept
        :
            ptr_(rhs.release())
        {}

        //- Move construct, transferring ownership from derived type.
        //  U must be derivable from T
        //  \note Future check std::enable_if + std::is_convertible ?
        template<class U>
        explicit autoPtr(autoPtr<U>&& rhs)
        :
            ptr_(rhs.release())
        {}

        //- Move construct from unique_ptr, transferring ownership.
        explicit autoPtr(std::unique_ptr<T>&& rhs) noexcept
        :
            ptr_(rhs.release())
        {}

        //- A move construct disguised as a copy construct (transfers ownership)
        //  \remark This should ideally be deleted - pending cleanup of code
        //      currently relying on this behaviour.
        #ifdef Foam_autoPtr_copyConstruct
        autoPtr(const autoPtr<T>& rhs) noexcept
        :
            ptr_(const_cast<autoPtr<T>&>(rhs).release())
        {}
        #else
        autoPtr(const autoPtr<T>&) = delete;
        #endif


    //- Destructor: deletes managed pointer
    ~autoPtr() noexcept
    {
        delete ptr_;
    }


    // Factory Methods

        //- Construct autoPtr with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_unique, but the overload for
        //      array types is not disabled.
        template<class... Args>
        static autoPtr<T> New(Args&&... args)
        {
            return autoPtr<T>(new T(std::forward<Args>(args)...));
        }

        //- Construct autoPtr from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types.
        //  Future check std::enable_if + std::is_convertible ?
        template<class U, class... Args>
        static autoPtr<T> NewFrom(Args&&... args)
        {
            return autoPtr<T>(new U(std::forward<Args>(args)...));
        }


    // Member Functions

    // Query

        //- True if the managed pointer is non-null
        bool good() const noexcept { return bool(ptr_); }


    // Access

        //- Return pointer to managed object without nullptr checking.
        //  Pointer remains under autoPtr management.
        T* get() noexcept { return ptr_; }

        //- Return const pointer to managed object without nullptr checking.
        //  Pointer remains under autoPtr management.
        const T* get() const noexcept { return ptr_; }

        //- Return reference to the managed object without nullptr checking.
        //  When get() == nullptr, additional guards may be required to avoid
        //  inadvertent access of a nullptr.
        T& ref() { return *ptr_; }


    // Edit

        //- Release ownership and return the pointer.
        //  \remark Method naming consistent with std::unique_ptr
        inline T* release() noexcept;

        //- Same as \c release().
        //  \remark Method naming consistent with Foam::tmp
        T* ptr() noexcept { return release(); }

        //- Same as \c reset(nullptr)
        //  \remark Method naming consistent with Foam::tmp
        void clear() noexcept { reset(nullptr); }

        //- Delete managed object and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Delete managed object and set to new given pointer
        //  \remark Same as move assign, but better for code documentation
        inline void reset(autoPtr<T>&& other) noexcept;

        //- Reset with emplace construction.
        //- Return reference to the new content.
        template<class... Args>
        inline T& emplace(Args&&... args);

        //- Swaps the managed object with other autoPtr.
        inline void swap(autoPtr<T>& other) noexcept;


    // Other

        //- Copy construct by invoking clone on underlying managed object
        //  A no-op if no pointer is managed
        //  \param args list of arguments for clone
        template<class... Args>
        inline autoPtr<T> clone(Args&&... args) const;


    // Member Operators

        //- Return reference to the managed object.
        //  FatalError if no pointer is managed
        inline T& operator*();

        //- Return const reference to the object.
        //  FatalError if no pointer is managed
        inline const T& operator*() const;

        //- Dereferences (non-const) pointer to the managed object
        //  FatalError if no pointer is managed
        inline T* operator->();

        //- Dereferences (const) pointer to the managed object
        //  FatalError if no pointer is managed
        inline const T* operator->() const;

        //- Return reference to the object data.
        //  FatalError if no pointer is managed
        inline T& operator()();

        //- Return const reference to the object data
        //  FatalError if no pointer is managed
        inline const T& operator()() const;


    // Casting

        //- True if pointer/reference is non-null. Same as good()
        explicit operator bool() const noexcept { return bool(ptr_); }

        //- Implicit cast to \em const pointer type.
        //  \note no cast to non-const pointer,
        //  since this would add uncertainty of ownership.
        operator const T*() const noexcept { return ptr_; }

        //- Deprecated(2019-01) Automatic cast conversion to underlying type
        //  FatalError if no pointer is managed
        //  \deprecated(2019-01) Can result in inadvertent conversions
        //      where the user should really know or check if the pointer
        //      is valid prior to using.
        #ifdef Foam_autoPtr_castOperator
        operator const T&() const { return operator*(); }
        #else
        operator const T&() const = delete;
        #endif


    // Assignment

        //- No copy assignment from plain pointer (uncontrolled access)
        void operator=(T* p) = delete;

        //- No move assignment disguised as a copy assignment
        //  \deprecated(2018-02) can have unintended behaviour
        void operator=(const autoPtr<T>&) = delete;

        //- Transfer object ownership from parameter
        void operator=(autoPtr<T>&& other) noexcept { reset(std::move(other)); }

        //- Transfer object ownership from parameter
        //  \note Future check std::enable_if + std::is_convertible ?
        template<class U>
        void operator=(autoPtr<U>&& other) noexcept
        {
            reset(other.release());
        }

        //- Transfer ownership by move assignment from unique_ptr
        void operator=(std::unique_ptr<T>&& other) { reset(other.release()); }

        //- Reset via assignment from literal nullptr
        void operator=(std::nullptr_t) noexcept { reset(nullptr); }


    // Housekeeping

        //- Identical to good(), or bool operator
        bool valid() const noexcept { return bool(ptr_); }

        //- Deprecated(2020-07) True if the managed pointer is null
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }

        //- Deprecated(2018-02) Identical to reset().
        //  \note Provided for backward compatibility - the older version
        //      enforced a run-time check (Fatal if pointer was already set)
        //      but this was rarely used.
        //  \deprecated(2018-02) Identical to reset().
        #ifndef Foam_autoPtr_nodeprecate_setMethod
        FOAM_DEPRECATED_FOR(2022-01, "reset() - same behaviour")
        #endif
        void set(T* p) noexcept { reset(p); }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Does not need std::swap or Foam::Swap() specialization
// since autoPtr is MoveConstructible and MoveAssignable


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
