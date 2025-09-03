/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::NullObject

Description
    Singleton null-object class and instance.

    Its contents occupy enough space to also be reinterpreted
    as another class with a null pointer or zero long for its first
    member, with additional zero parameters for safe casting to List etc.

SourceFiles
    nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_nullObject_H
#define Foam_nullObject_H

#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class NullObject;

/*---------------------------------------------------------------------------*\
                         Class NullObject Declaration
\*---------------------------------------------------------------------------*/

class NullObject
{
    //- A %union of zero data types
    union zeros
    {
        void* ptr;
        unsigned long val;
    };


    // Private Data

        //- The zero data content
        zeros data_[4];


    // Constructors

        //- Private constructor for singleton only
        //  Could also rely on bit-wise zero initialization for union content
        NullObject()
        :
           data_{{nullptr}, {nullptr}, {nullptr}, {nullptr}}
        {}

        //- No copy construct
        NullObject(const NullObject&) = delete;

        //- No copy assignment
        void operator=(const NullObject&) = delete;


public:

    // Static Data

        //- A unique null object
        static const NullObject nullObject;


    // Member Functions

        //- A nullptr pointer content
        const void* pointer() const noexcept
        {
            return data_[0].ptr;
        }

        //- Zero valued integer content
        unsigned long value() const noexcept
        {
            return 0;
        }

        //- No elements
        bool empty() const noexcept
        {
            return true;
        }

        //- Zero elements
        label size() const noexcept
        {
            return 0;
        }

        //- No-op method (for HashTable replacement)
        const NullObject& toc() const noexcept
        {
            return *this;
        }

        //- No-op method (for HashTable replacement)
        const NullObject& sortedToc() const noexcept
        {
            return *this;
        }


    // Member Operators

        //- Swallow assignment (cf, std::ignore)
        template<class T>
        const NullObject& operator=(const T&) const noexcept
        {
            return *this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Globals

//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


// IOstream Operators

//- Read from Istream consumes no content, does not change NullObject
inline Istream& operator>>(Istream& is, const NullObject&) noexcept
{
    return is;
}

//- Write to Ostream emits no content
inline Ostream& operator<<(Ostream& os, const NullObject&) noexcept
{
    return os;
}


// Global Functions

//- Const pointer (of type T) to the nullObject
template<class T>
inline const T* NullObjectPtr() noexcept
{
    return reinterpret_cast<const T*>(nullObjectPtr);
}

//- Non-const pointer (of type T) to the nullObject.
//- Only use when nothing will be written into it!
template<class T>
inline T* NullObjectPtr_constCast() noexcept
{
    return reinterpret_cast<T*>(const_cast<NullObject*>(nullObjectPtr));
}


//- Const reference (of type T) to the nullObject
template<class T>
inline const T& NullObjectRef() noexcept
{
    return *reinterpret_cast<const T*>(nullObjectPtr);
}

//- Non-const reference (of type T) to the nullObject
//- Only use when nothing will be written into it!
template<class T>
inline T& NullObjectRef_constCast() noexcept
{
    return *reinterpret_cast<T*>(const_cast<NullObject*>(nullObjectPtr));
}


//- True if ptr is a pointer (of type T) to the nullObject
template<class T>
inline bool isNull(const T* ptr) noexcept
{
    return ptr == NullObjectPtr<T>();
}

//- True if obj is a reference (of type T) to the nullObject
template<class T>
inline bool isNull(const T& obj) noexcept
{
    return &obj == NullObjectPtr<T>();
}


//- True if ptr is not a pointer (of type T) to the nullObject
template<class T>
inline bool notNull(const T* ptr) noexcept
{
    return ptr != NullObjectPtr<T>();
}

//- True if obj is not a reference (of type T) to the nullObject
template<class T>
inline bool notNull(const T& obj) noexcept
{
    return &obj != NullObjectPtr<T>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
