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

#include "UPstream.H"
#include <cstring>  // memmove

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native)                                        \
                                                                              \
void Foam::UPstream::mpiGather                                                \
(                                                                             \
    const Native* sendData,                                                   \
    Native* recvData,                                                         \
    int count,                                                                \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    std::memmove(recvData, sendData, count*sizeof(Native));                   \
}                                                                             \
                                                                              \
                                                                              \
void Foam::UPstream::mpiScatter                                               \
(                                                                             \
    const Native* sendData,                                                   \
    Native* recvData,                                                         \
    int count,                                                                \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    std::memmove(recvData, sendData, count*sizeof(Native));                   \
}                                                                             \
                                                                              \
                                                                              \
void Foam::UPstream::mpiAllGather                                             \
(                                                                             \
    Native* allData,                                                          \
    int count,                                                                \
    const label comm                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
                                                                              \
void Foam::UPstream::gather                                                   \
(                                                                             \
    const Native* sendData,                                                   \
    int sendCount,                                                            \
                                                                              \
    Native* recvData,                                                         \
    const UList<int>& recvCounts,                                             \
    const UList<int>& recvOffsets,                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    /* recvCounts[0] may be invalid - use sendCount instead */                \
    std::memmove(recvData, sendData, sendCount*sizeof(Native));               \
}                                                                             \
                                                                              \
void Foam::UPstream::scatter                                                  \
(                                                                             \
    const Native* sendData,                                                   \
    const UList<int>& sendCounts,                                             \
    const UList<int>& sendOffsets,                                            \
                                                                              \
    Native* recvData,                                                         \
    int recvCount,                                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    std::memmove(recvData, sendData, recvCount*sizeof(Native));               \
}


//TDB: Pstream_CommonRoutines(bool);
Pstream_CommonRoutines(char);
Pstream_CommonRoutines(int32_t);
Pstream_CommonRoutines(int64_t);
Pstream_CommonRoutines(uint32_t);
Pstream_CommonRoutines(uint64_t);
Pstream_CommonRoutines(float);
Pstream_CommonRoutines(double);

#undef Pstream_CommonRoutines

// ************************************************************************* //
