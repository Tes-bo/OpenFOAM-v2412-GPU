/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool& value, const label comm)
{}

void Foam::UPstream::reduceOr(bool& value, const label comm)
{}


void Foam::reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag,
    const label comm
)
{}

void Foam::reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag,
    const label comm
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Common reductions

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native)                                      \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const minOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const maxOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Floating-point reductions

#undef  Pstream_FloatReductions
#define Pstream_FloatReductions(Native)                                       \
                                                                              \
Pstream_CommonReductions(Native);                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    UPstream::Request& req                                                    \
)                                                                             \
{}                                                                            \
                                                                              \
/* Deprecated: prefer version with UPstream::Request */                       \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    UPstream::Request& req                                                    \
)                                                                             \
{}                                                                            \
                                                                              \
/* Deprecated: prefer version with UPstream::Request */                       \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,                                                            \
    const label comm,                                                         \
    label& requestID                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::sumReduce                                                          \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Bitwise reductions

#undef  Pstream_BitwiseReductions
#define Pstream_BitwiseReductions(Native)                                     \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const bitOrOp<Native>&,                                                   \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const bitOrOp<Native>&,                                                   \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Pstream_CommonReductions(int32_t);
Pstream_CommonReductions(int64_t);
Pstream_CommonReductions(uint32_t);
Pstream_CommonReductions(uint64_t);

Pstream_FloatReductions(float);
Pstream_FloatReductions(double);

Pstream_BitwiseReductions(unsigned char);
Pstream_BitwiseReductions(unsigned int);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef Pstream_CommonReductions
#undef Pstream_FloatReductions
#undef Pstream_BitwiseReductions


// ************************************************************************* //
