/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Includes some standard C++ headers, defines global macros and templates
    used in multiple places by OpenFOAM.

\*---------------------------------------------------------------------------*/

#ifndef Foam_stdFoam_H
#define Foam_stdFoam_H

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <utility>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compile-time warning for use of deprecated methods (compiler-dependent).
// Use within the class declaration.

#if (__cplusplus >= 201402L)
# define FOAM_DEPRECATED(since) [[deprecated("Since " #since)]]
# define FOAM_DEPRECATED_FOR(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
# define FOAM_DEPRECATED_STRICT(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
#elif defined(__GNUC__)
# define FOAM_DEPRECATED(since) __attribute__((__deprecated__("Since " #since)))
# define FOAM_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))
# define FOAM_DEPRECATED_STRICT(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))
#else
# define FOAM_DEPRECATED(since)
# define FOAM_DEPRECATED_FOR(since, replacement)
# define FOAM_DEPRECATED_STRICT(since, replacement)
#endif
#ifndef FOAM_COMPILE_STRICT
# undef  FOAM_DEPRECATED_STRICT
# define FOAM_DEPRECATED_STRICT(since, replacement)
#endif

// Compile-time warning about unused result
// FUTURE: check __has_cpp_attribute(nodiscard) and define with [[nodiscard]]
#if defined(__GNUC__)
# define FOAM_NODISCARD __attribute__((warn_unused_result))
#else
# define FOAM_NODISCARD
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Namespace for OpenFOAM
namespace Foam
{

//- Implementation details for various OpenFOAM classes
namespace Detail {}

//- Additional OpenFOAM modules
namespace Module {}


//- A functor that returns its argument unchanged (cf. C++20 std::identity)
//- Should \em never be specialized.
struct identityOp
{
    using is_transparent = void;

    template<class T>
    constexpr T&& operator()(T&& val) const noexcept
    {
        return std::forward<T>(val);
    }

    // Allow use as an identity array/map
    template<class T>
    constexpr T&& operator[](T&& val) const noexcept
    {
        return std::forward<T>(val);
    }
};


//- Swap non-array types as per std::swap, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
//
//  \note For complex structures, it is usual to provide a swap() member
//  function and specialize Swap()
template<class T>
void Swap(T& a, T& b)
{
    std::swap(a, b);
}


//- Swap array types as per std::swap example, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
template<class T, size_t N>
void Swap(T (&a)[N], T (&b)[N])
{
    for (size_t i = 0; i < N; ++i)
    {
        Foam::Swap(a[i], b[i]);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Namespace for std templates that are are part of future C++ standards
//- or that are in a state of change.
//
// SeeAlso
//     - https://en.cppreference.com/w/cpp/iterator/begin
//     - https://en.cppreference.com/w/cpp/iterator/end
//     - https://en.cppreference.com/w/cpp/iterator/rbegin
//     - https://en.cppreference.com/w/cpp/iterator/rend

namespace stdFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward iteration

//- Return iterator to the beginning of the container \a c.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c.
//  Definition as per std::cbegin C++17
template<class C>
constexpr auto cbegin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return iterator to the end of the container \a c.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(const C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c.
//  Definition as per std::cend C++17
template<class C>
constexpr auto cend(const C& c) -> decltype(c.end())
{
    return c.end();
}


// Reverse iteration

//- Return reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::rbegin C++17
template<class C>
constexpr auto rbegin(C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return const_reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::rbegin C++17
template<class C>
constexpr auto rbegin(const C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return const_reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::crbegin C++17
template<class C>
constexpr auto crbegin(const C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return reverse_iterator to reverse-end of container \a c.
//  Definition as per std::rend C++17
template<class C>
constexpr auto rend(C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return const_reverse_iterator to reverse-end of container \a c.
//  Definition as per std::rend C++17
template<class C>
constexpr auto rend(const C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return const_reverse_iterator to reverse-end of container \a c.
//  Definition as per std::crend C++17
template<class C>
constexpr auto crend(const C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return the lesser of the parameters.
//  Definition as per std::min C++14
template<class T>
constexpr inline const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}

//- Return the greater of the parameters.
//  Definition as per std::max C++14
template<class T>
constexpr inline const T& max(const T& a, const T& b)
{
    return (a < b) ? b : a;
}

} // End namespace stdFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Iterate across all elements in the \a container object.
//  \par Usage
//  \code
//  forAllIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIters, forAllIter, forAllConstIters
#define forAllIters(container,iter)                                            \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::begin(container);                                 \
        iter != stdFoam::end(container);                                       \
        ++iter                                                                 \
    )


//- Iterate across all elements of the \a container object with const access.
//  \par Usage
//  \code
//  forAllConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIters, forAllIter, forAllConstIter
#define forAllConstIters(container,iter)                                       \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::cbegin(container);                                \
        iter != stdFoam::cend(container);                                      \
        ++iter                                                                 \
    )


//- Reverse iterate across elements in the \a container object of type
//  \a Container.
//  \par Usage
//  \code
//  forAllReverseIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstReverseIters
#define forAllReverseIters(container,iter)                                     \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::rbegin(container);                                \
        iter != stdFoam::rend(container);                                      \
        ++iter                                                                 \
    )


//- Reverse iterate across elements of \a container object with const access.
//  \par Usage
//  \code
//  forAllReverseConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllReverseIters
#define forAllConstReverseIters(container,iter)                                \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::crbegin(container);                               \
        iter != stdFoam::crend(container);                                     \
        ++iter                                                                 \
    )


//- Loop across all elements in \a list
// \par Usage
// \code
// forAll(anyList, i)
// {
//      statements;
// }
// \endcode
// \sa forAllReverse
#define forAll(list, i) \
    for (Foam::label i=0; i<(list).size(); ++i)


//- Reverse loop across all elements in \a list
//  \par Usage
//  \code
//  forAllReverse(anyList, i)
//  {
//       statements;
//  }
//  \endcode
//  \sa forAll
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; --i)


// Compatibility macros for pre C++11

//- Iterate across all elements in the \a container object
//  of type \a Container.
//  \par Usage
//  \code
//  forAllIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIter
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin();                        \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )


//- Iterate across all elements in the \a container object
//  of type \a Container with const access.
//  \par Usage
//  \code
//  forAllConstIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIter
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).cbegin();                 \
        iter != (container).cend();                                            \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace stdFoam
{

/*---------------------------------------------------------------------------*\
                        Class stdFoam::span Declaration
\*---------------------------------------------------------------------------*/

//- Rudimentary functionality similar to std::span for holding memory view
template<class Type>
class span
{
    // Private Data

        //- The data pointer
        Type* data_;

        //- The size of the span
        std::size_t size_;


public:

    // STL type definitions

        using element_type = Type;
        using value_type = std::remove_cv<Type>;
        using size_type = std::size_t;
        using pointer = Type*;
        using const_pointer = const Type*;
        using reference = Type&;
        using const_reference = const Type&;
        using difference_type = std::ptrdiff_t;
        using iterator = Type*;


    // Generated Methods

        //- Copy construct
        span(const span& other) noexcept = default;

        //- Copy assignment
        span& operator=(const span& other) noexcept = default;


    // Constructors

        //- Default construct
        constexpr span() noexcept
        :
            data_(nullptr),
            size_(0)
        {}

        //- Construct from pointer and size
        constexpr span(pointer ptr, size_type count) noexcept
        :
            data_(ptr),
            size_(count)
        {}

        //- Construct from begin/end pointers
        constexpr span(pointer first, pointer last) noexcept
        :
            data_(first),
            size_(last - first)
        {}


    //- Destructor
    ~span() noexcept = default;


    // Member Functions

    // Iterators

        //- Iterator to begin of span
        constexpr iterator begin() const noexcept { return data_; }

        //- Iterator to one-past end of span
        constexpr iterator end() const noexcept { return (data_ + size_); }


    // Element access

        //- Access the first element. Undefined if span is empty
        constexpr reference front() const { return *(data_); }

        //- Access the last element. Undefined if span is empty
        constexpr reference back() const { return *(data_ + size_ - 1); }

        //- Access an element of the sequence
        constexpr reference operator[](size_type idx) const
        {
           return *(data_ + idx);
        }

        //- Return a pointer to the beginning of the sequence
        constexpr pointer data() const noexcept { return data_; }


    // Observers

        //- Number of elements in the sequence
        constexpr size_type size() const noexcept { return size_; }

        //- The size of the sequence in bytes
        constexpr size_type size_bytes() const noexcept
        {
            return (size_*sizeof(Type));
        }

        //- True if the sequence is empty
        constexpr bool empty() const noexcept { return !size_; }


    // Subview

        //- Obtains a span of the first count elements
        span<Type> first(size_type count) const noexcept
        {
            return span<Type>(data_, count);
        }

        //- Obtains a span of the last count elements
        span<Type> last(size_type count) const noexcept
        {
            return span<Type>(data_ + (size_ - count), count);
        }

        //- Obtains a sub-span starting at pos until end of the sequence
        span<Type> subspan(size_type pos) const noexcept
        {
            return span<Type>(data_ + pos, size_ - pos);
        }

        //- Obtains a sub-span of length len starting at pos.
        //  Graciously handles excess lengths.
        span<Type> subspan(size_type pos, size_type len) const noexcept
        {
            return span<Type>(data_ + pos, std::min(size_ - pos, len));
        }


    // Additional members, similar to UList etc.
    // The std::span has as_bytes() and as_writeable_bytes() as free functions

        //- A readonly view as byte content
        constexpr const char* cdata_bytes() const noexcept
        {
            return reinterpret_cast<const char*>(data_);
        }

        //- A writable view as byte content (if the pointer type is non-const).
        //- Like data(), the const access itself is const.
        template<class TypeT = Type>
        typename std::enable_if<!std::is_const<TypeT>::value, char*>::type
        data_bytes() const noexcept
        {
            return reinterpret_cast<char*>(data_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stdFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
