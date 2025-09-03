/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#ifdef __GNUC__
    #define USE_BUILTIN_BYTESWAP
#else
    #undef USE_BUILTIN_BYTESWAP
#endif

// for Debugging:
// #undef USE_BUILTIN_BYTESWAP


// Some ideas about a templated approach for swapping bytes:
//     http://www.cplusplus.com/forum/general/27544/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::endian::isBig() noexcept
{
    const short testBig = 0x0100;

    // Yields 0x01 for big endian
    return *(reinterpret_cast<const char*>(&testBig));
}


inline bool Foam::endian::isLittle() noexcept
{
    const short testLittle = 0x0001;

    // Yields 0x01 for little endian
    return *(reinterpret_cast<const char*>(&testLittle));
}


inline uint32_t Foam::endian::swap32(uint32_t u)
{
#ifdef USE_BUILTIN_BYTESWAP
    return __builtin_bswap32(u);
#else
    return
    (
        (((u) & 0xff000000) >> 24)      // 3 -> 0
      | (((u) & 0x00ff0000) >>  8)      // 2 -> 1
      | (((u) & 0x0000ff00) <<  8)      // 2 <- 1
      | (((u) & 0x000000ff) << 24)      // 3 <- 0
    );

    // Alternative formulation
    //
    // u = ((u<<8) & 0xFF00FF00) | ((u>>8) & 0x00FF00FF);
    // u = (u>>16) | (u<<16);
    // return u;
#endif
}


inline uint64_t Foam::endian::swap64(uint64_t u)
{
#ifdef USE_BUILTIN_BYTESWAP
    return __builtin_bswap64(u);
#else
    return
    (
        (((u) & 0xff00000000000000ull) >> 56)   // 7 -> 0
      | (((u) & 0x00ff000000000000ull) >> 40)   // 6 -> 1
      | (((u) & 0x0000ff0000000000ull) >> 24)   // 5 -> 2
      | (((u) & 0x000000ff00000000ull) >>  8)   // 4 -> 3
      | (((u) & 0x00000000ff000000ull) <<  8)   // 4 <- 3
      | (((u) & 0x0000000000ff0000ull) << 24)   // 5 <- 2
      | (((u) & 0x000000000000ff00ull) << 40)   // 6 <- 1
      | (((u) & 0x00000000000000ffull) << 56)   // 7 <- 0
    );

    // Alternative formulation
    /*
     u = ((u<< 8) & 0xFF00FF00FF00FF00ull) | ((u>> 8) & 0x00FF00FF00FF00FFull);
     u = ((u<<16) & 0xFFFF0000FFFF0000ull) | ((u>>16) & 0x0000FFFF0000FFFFull);
     return (u >> 32) | (u << 32);
    */
#endif
}


// ************************************************************************* //
