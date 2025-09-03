/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::Hash

Description
    Hash function class.
    The default definition is for primitives.
    Non-primitives used to hash entries on hash tables will need
    a specialized version.

\*---------------------------------------------------------------------------*/

#ifndef Hash_H
#define Hash_H

#include "Hasher.H"
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class T>
struct Hash
{
    unsigned operator()(const T& obj, unsigned seed=0) const
    {
        return Foam::Hasher(&obj, sizeof(T), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing of nullptr, always 0
template<>
struct Hash<std::nullptr_t>
{
    unsigned operator()(std::nullptr_t, unsigned seed=0) const noexcept
    {
        return seed;
    }
};

//- Hashing of pointers, treat as unsigned integer
template<>
struct Hash<void*>
{
    unsigned operator()(const void* const ptr, unsigned seed=0) const
    {
        const uintptr_t addr = uintptr_t(ptr);
        return Foam::Hasher(&addr, sizeof(addr), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialization for common integral types

#undef  FOAM_HASH_SPECIALIZATION
#define FOAM_HASH_SPECIALIZATION(Type)                                        \
                                                                              \
    /*! \brief Hashing of integral type: Type */                              \
    /*! Unseeded (single value) uses natural order, otherwise incremental */  \
    template<>                                                                \
    struct Hash<Type>                                                         \
    {                                                                         \
        unsigned operator()(const Type val) const                             \
        {                                                                     \
            return static_cast<unsigned>(val);                                \
        }                                                                     \
        unsigned operator()(const Type val, unsigned seed) const              \
        {                                                                     \
            return Foam::Hasher(&val, sizeof(Type), seed);                    \
        }                                                                     \
    }

FOAM_HASH_SPECIALIZATION(bool);
FOAM_HASH_SPECIALIZATION(char);
FOAM_HASH_SPECIALIZATION(int32_t);
FOAM_HASH_SPECIALIZATION(int64_t);
FOAM_HASH_SPECIALIZATION(uint32_t);
#undef FOAM_HASH_SPECIALIZATION


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
