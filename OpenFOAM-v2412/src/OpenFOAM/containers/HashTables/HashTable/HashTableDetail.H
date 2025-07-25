/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::Detail::HashTablePair

Description
    Internal storage type for HashTable

SourceFiles
    HashTableDetail.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_HashTableDetail_H
#define Foam_HashTableDetail_H

#include "zero.H"
#include <memory>
#include <utility>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Ostream;
template<class T> class autoPtr;

namespace Detail
{

/*---------------------------------------------------------------------------*\
                        Class isPointerLike Declaration
\*---------------------------------------------------------------------------*/

//- Pointer-like behaviour
template<class T> struct isPointerLike : std::false_type {};

//- Pointer-like behaviour for autoPtr
template<class T> struct isPointerLike<autoPtr<T>> : std::true_type {};

//- Pointer-like behaviour for std::unique_ptr
template<class T> struct isPointerLike<std::unique_ptr<T>> : std::true_type {};


/*---------------------------------------------------------------------------*\
                        Class HashTablePair Declaration
\*---------------------------------------------------------------------------*/

//- Internal storage type for HashTable entries
//  Comprises a (K,V) tuple and a linked-list entry for hash collisions.
//  Could store key/val as std::pair, but no particular advantage
//  unless the iterator dereference type changes.
template<class K, class V>
struct HashTablePair
{
    // Data Types

        //- Type of key
        typedef K key_type;

        //- Type of content
        typedef V mapped_type;

        //- This class stores a value
        static constexpr bool stores_value() noexcept
        {
            return true;
        }


    // Data

        //- The lookup key
        const K key_;

        //- The mapped data
        V val_;

        //- Addressing (next in collision list)
        HashTablePair* next_;


    // Generated Methods

        //- No copy construct
        HashTablePair(const HashTablePair&) = delete;

        //- No copy assignment
        void operator=(const HashTablePair&) = delete;


    // Constructors

        //- Construct from next pointer, key, contents
        template<class... Args>
        HashTablePair
        (
            HashTablePair* next,
            const key_type& key,
            Args&&... args
        )
        :
            key_(key),
            val_(std::forward<Args>(args)...),
            next_(next)
        {}


    // Member Functions

        //- The key
        const key_type& key() const noexcept
        {
            return key_;
        }

        //- Const access to the mapped value
        const mapped_type& cval() const noexcept
        {
            return val_;
        }

        //- Const access to the mapped value
        const mapped_type& val() const noexcept
        {
            return val_;
        }

        //- Non-const access to the mapped value
        mapped_type& val() noexcept
        {
            return val_;
        }


        //- Compare keys
        bool operator<(const HashTablePair& rhs) const
        {
            return key_ < rhs.key_;
        }

        //- Write (key, val) pair - for pointer types
        template<class TypeT = V>
        typename std::enable_if
        <
            (
                std::is_pointer<TypeT>::value
             || Detail::isPointerLike<TypeT>::value
            ),
            void
        >::type
        print(Ostream& os) const
        {
            os << key_;

            if (val_)
            {
                os << ' ' << *val_;
            }
        }

        //- Write (key, val) pair - for non-pointer types
        template<class TypeT = V>
        typename std::enable_if
        <
            (
                !std::is_pointer<TypeT>::value
             && !Detail::isPointerLike<TypeT>::value
            ),
            void
        >::type
        print(Ostream& os) const
        {
            os << key_ << ' ' << val_;
        }
};


/*---------------------------------------------------------------------------*\
                       Class HashTableSingle Declaration
\*---------------------------------------------------------------------------*/

//- Internal storage type for HashSet entries
//  Comprises a single (K) value and a linked-list entry for hash collisions
template<class K>
struct HashTableSingle
{
    // Data Types

        //- Type of key
        typedef K key_type;

        //- Type of content (no content: placeholder)
        typedef Foam::zero mapped_type;

        //- Content storage type to the entry
        typedef key_type value_type;

        //- This class does not store any value
        static constexpr bool stores_value() noexcept
        {
            return false;
        }


    // Data

        //- The lookup key == content
        const K key_;

        //- Addressing (next in collision list)
        HashTableSingle* next_;


    // Generated Methods

        //- No copy construct
        HashTableSingle(const HashTableSingle&) = delete;

        //- No copy assignment
        void operator=(const HashTableSingle&) = delete;


    // Constructors

        //- Construct from next pointer, key, (unused) contents
        template<class... Args>
        HashTableSingle
        (
            HashTableSingle* next,
            const key_type& key,
            Args&&...
        )
        :
            key_(key),
            next_(next)
        {}


    // Member Functions

        //- The key
        const key_type& key() const noexcept
        {
            return key_;
        }

        //- Const access to the (dummy) mapped value (for API only)
        const mapped_type& cval() const noexcept
        {
            return Foam::zero::dummy;
        }

        //- Const access to the (dummy) mapped value (for API only)
        const mapped_type& val() const noexcept
        {
            return Foam::zero::dummy;
        }

        //- Non-const access to the (dummy) mapped value (for API only)
        mapped_type& val() noexcept
        {
            return Foam::zero::dummy;
        }


        //- Compare keys
        bool operator<(const HashTableSingle& rhs) const
        {
            return key_ < rhs.key_;
        }

        //- Write the key. There is no value to write.
        void print(Ostream& os) const
        {
            os << key_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
