/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "UIPstream.H"
#include "IPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPBstream::UIPBstream
(
    const UPstream::commsTypes commsType,
    const int rootProcNo,
    DynamicList<char>& receiveBuf,
    label& receiveBufPosition,
    const int tag,
    const label comm,
    const bool clearAtEnd,
    IOstreamOption::streamFormat fmt
)
:
    UIPstreamBase
    (
        commsType,              // irrelevant
        rootProcNo,             // normally UPstream::masterNo()
        receiveBuf,
        receiveBufPosition,
        tag,                    // irrelevant
        comm,
        clearAtEnd,
        fmt
    )
{
    bufferIPCrecv();
}


Foam::IPBstream::IPBstream
(
    const UPstream::commsTypes commsType,
    const int rootProcNo,
    const label bufSize,
    const int tag,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    Pstream(commsType, bufSize),
    UIPBstream
    (
        commsType,              // irrelevant
        rootProcNo,             // normally UPstream::masterNo()
        Pstream::transferBuf_,
        UIPstreamBase::storedRecvBufPos_,   // Internal only
        tag,                    // irrelevant
        comm,
        false,  // Do not clear Pstream::transferBuf_ if at end
        fmt
    )
{}


Foam::IPBstream::IPBstream
(
    const int rootProcNo,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    IPBstream
    (
        UPstream::commsTypes::scheduled,    // irrelevant
        rootProcNo,
        label(0),  // bufSize
        UPstream::msgType(),                // irrelevant
        comm,
        fmt
    )
{}


Foam::IPBstream::IPBstream
(
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    IPBstream
    (
        UPstream::commsTypes::scheduled,    // irrelevant
        UPstream::masterNo(),  // rootProcNo
        label(0),  // bufSize
        UPstream::msgType(),                // irrelevant
        comm,
        fmt
    )
{}


// ************************************************************************* //
