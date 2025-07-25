/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::HashFunction

Description
    Hash function class.

    Verify that template overloads are properly resolved

Note
    The second template parameter (bool) is used for SFINAE overloading,

\*---------------------------------------------------------------------------*/

#ifndef HashFunction_H
#define HashFunction_H

#include "Hash.H"

#ifdef FULLDEBUG
#define HashTypeInfo(Args)  void info() { std::cerr<< "" Args << "\n"; }
#else
#define HashTypeInfo(Args)  void info() {}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class T, class SFINAEType=bool>
struct HashFun
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "default"; }
    #endif
    HashTypeInfo("plain hash")

    unsigned operator()(const T& obj, unsigned seed=0) const
    {
        return Foam::Hasher(&obj, sizeof(obj), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing for label
template<> struct HashFun<Foam::label> : Hash<label>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "label"; }
    #endif
    HashTypeInfo("hash label")
};


//- Hashing for pointers, interpret pointer as a integer type
template<> struct HashFun<void*> : Hash<void *>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "pointer"; }
    #endif
    HashTypeInfo("hash ptr")
};


//- Hashing for string types
template<class StringType>
struct HashFun
<
    StringType,
    typename std::enable_if
    <
        std::is_base_of<std::string, StringType>::value, bool
    >::type
>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "string"; }
    #endif
    HashTypeInfo("hash string")

    unsigned operator()(const std::string& obj, unsigned seed=0) const
    {
        return Foam::Hasher(obj.data(), obj.size(), seed);
    }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Various
#include "edge.H"
#include "face.H"
#include "triFace.H"
#include "Pair.H"
#include "Tuple2.H"
#include "DynamicList.H"
#include "FixedList.H"

namespace Foam
{

template<> struct HashFun<edge> : Hash<edge>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "edge"; }
    #endif
    HashTypeInfo("hash edge")
};


template<> struct HashFun<face> : Hash<face>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "face"; }
    #endif
    HashTypeInfo("hash face")
};


template<> struct HashFun<triFace> : Hash<triFace>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "triFace"; }
    #endif
    HashTypeInfo("hash triFace")
};


template<class T>
struct HashFun<Pair<T>> : Hash<Pair<T>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "Pair"; }
    #endif
    HashTypeInfo("hash Pair")
};

template<class T1, class T2>
struct HashFun<Tuple2<T1, T2>> : Hash<Tuple2<T1, T2>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "Tuple2"; }
    #endif
    HashTypeInfo("hash Tuple2")
};


template<class T>
struct HashFun<List<T>> : Hash<List<T>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "List"; }
    #endif
    HashTypeInfo("hash List")
};

template<class T> struct HashFun<UList<T>> : Hash<UList<T>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "UList"; }
    #endif
    HashTypeInfo("hash UList")
};

template<class T, int SizeMin>
struct HashFun<DynamicList<T, SizeMin>> : Hash<DynamicList<T, SizeMin>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "DynamicList"; }
    #endif
    HashTypeInfo("hash DynamicList")
};

template<class T, unsigned N>
struct HashFun<FixedList<T, N>> : Hash<FixedList<T, N>>
{
    #ifdef FULLDEBUG
    static constexpr const char* name() noexcept { return "FixedList"; }
    #endif
    HashTypeInfo("hash FixedList")
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
