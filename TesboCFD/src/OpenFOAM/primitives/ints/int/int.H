/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

Primitive
    int

Description
    System signed integer

SourceFiles
    intIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_int_H
#define Foam_primitives_int_H

#include "int8.H"
#include "int16.H"
#include "int32.H"
#include "int64.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define MAXMIN(RetType, Type1, Type2)                                          \
                                                                               \
/*! \brief Floating/integral min. Use stdFoam::min() to preserve references */ \
inline RetType min(const Type1 s1, const Type2 s2)                             \
{                                                                              \
    return (s1 < s2) ? s1 : s2;                                                \
}                                                                              \
                                                                               \
/*! \brief Floating integral max. Use stdFoam::max() to preserve references */ \
inline RetType max(const Type1 s1, const Type2 s2)                             \
{                                                                              \
    return (s2 < s1) ? s1 : s2;                                                \
}


MAXMIN(int8_t, int8_t, int8_t)
MAXMIN(int16_t, int16_t, int16_t)

MAXMIN(int32_t, int32_t, int32_t)
MAXMIN(int64_t, int64_t, int32_t)
MAXMIN(int64_t, int32_t, int64_t)
MAXMIN(int64_t, int64_t, int64_t)


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read int from stream
int readInt(Istream& is);

//- Parse entire buffer as an int, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
int readInt(const char* buf);

//- Parse entire string as an int, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline int readInt(const std::string& str)
{
    return readInt(str.c_str());
}

//- Read entire buffer as an int, skipping leading/trailing whitespace.
//  \return True if successful.
bool readInt(const char* buf, int& val);

//- Read entire string as an int32_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readInt(const std::string& str, int& val)
{
    return readInt(str.c_str(), val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
