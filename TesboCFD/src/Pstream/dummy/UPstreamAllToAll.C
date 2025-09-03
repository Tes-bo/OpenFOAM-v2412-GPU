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
#include "Map.H"

#include <cinttypes>
#include <cstring>  // memmove

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native)                                        \
void Foam::UPstream::allToAll                                                 \
(                                                                             \
    const UList<Native>& sendData,                                            \
    UList<Native>& recvData,                                                  \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    recvData.deepCopy(sendData);                                              \
}


Pstream_CommonRoutines(int32_t);
Pstream_CommonRoutines(int64_t);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native)                                        \
void Foam::UPstream::allToAllConsensus                                        \
(                                                                             \
    const UList<Native>& sendData,                                            \
    UList<Native>& recvData,                                                  \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    recvData.deepCopy(sendData);                                              \
}                                                                             \
void Foam::UPstream::allToAllConsensus                                        \
(                                                                             \
    const Map<Native>& sendData,                                              \
    Map<Native>& recvData,                                                    \
    const int tag,                                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    recvData = sendData;                                                      \
}


Pstream_CommonRoutines(int32_t);
Pstream_CommonRoutines(int64_t);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native)                                        \
void Foam::UPstream::allToAll                                                 \
(                                                                             \
    const Native* sendData,                                                   \
    const UList<int>& sendCounts,                                             \
    const UList<int>& sendOffsets,                                            \
    Native* recvData,                                                         \
    const UList<int>& recvCounts,                                             \
    const UList<int>& recvOffsets,                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    if (recvCounts[0] != sendCounts[0])                                       \
    {                                                                         \
        FatalErrorInFunction                                                  \
            << "Number to send " << sendCounts[0]                             \
            << " does not equal number to receive " << recvCounts[0]          \
            << Foam::abort(FatalError);                                       \
    }                                                                         \
    std::memmove(recvData, sendData, recvCounts[0]*sizeof(Native));           \
}


// Unused: Pstream_CommonRoutines(char);

#undef Pstream_CommonRoutines

// ************************************************************************* //
