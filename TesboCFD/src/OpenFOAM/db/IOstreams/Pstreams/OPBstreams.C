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

#include "UOPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPBstream::UOPBstream
(
    const UPstream::commsTypes commsType,
    const int rootProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const label comm,
    const bool sendAtDestruct,
    IOstreamOption::streamFormat fmt
)
:
    UOPstreamBase
    (
        commsType,              // irrelevant
        rootProcNo,             // normally UPstream::masterNo()
        sendBuf,
        tag,                    // irrelevant
        comm,
        sendAtDestruct,
        fmt
    )
{}


Foam::OPBstream::OPBstream
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
    UOPBstream
    (
        commsType,              // irrelevant
        rootProcNo,             // normally UPstream::masterNo()
        Pstream::transferBuf_,
        tag,                    // irrelevant
        comm,
        true,  // sendAtDestruct
        fmt
    )
{}


Foam::OPBstream::OPBstream
(
    const int rootProcNo,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    OPBstream
    (
        UPstream::commsTypes::scheduled,    // irrelevant
        rootProcNo,
        label(0),  // bufSize
        UPstream::msgType(),                // irrelevant
        comm,
        fmt
    )
{}


Foam::OPBstream::OPBstream
(
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    OPBstream
    (
        UPstream::commsTypes::scheduled,    // irrelevant
        UPstream::masterNo(),  // rootProcNo
        label(0),  // bufSize
        UPstream::msgType(),                // irrelevant
        comm,
        fmt
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPBstream::~UOPBstream()
{
    if (sendAtDestruct_)
    {
        if (!bufferIPCsend())
        {
            FatalErrorInFunction
                << "Failed broadcast message of size "
                << sendBuf_.size() << " root: " << toProcNo_
                << Foam::abort(FatalError);
        }
    }
}


// ************************************************************************* //
