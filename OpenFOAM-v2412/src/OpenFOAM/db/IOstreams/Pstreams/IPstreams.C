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

#include "UIPstream.H"
#include "IPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
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
        commsType,
        fromProcNo,
        receiveBuf,
        receiveBufPosition,
        tag,
        comm,
        clearAtEnd,
        fmt
    )
{
    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Message is already received into buffer
    }
    else
    {
        bufferIPCrecv();
    }
}


Foam::UIPstream::UIPstream(const int fromProcNo, PstreamBuffers& buffers)
:
    UIPstreamBase(fromProcNo, buffers)
{
    if (commsType() == UPstream::commsTypes::nonBlocking)
    {
        // Message is already received into buffer
        messageSize_ = recvBuf_.size();

        if (debug)
        {
            Perr<< "UIPstream::UIPstream PstreamBuffers :"
                << " fromProcNo:" << fromProcNo_
                << " tag:" << tag_ << " comm:" << comm_
                << " receive buffer size:" << messageSize_
                << Foam::endl;
        }
    }
    else
    {
        bufferIPCrecv();
    }
}


Foam::UIPstream::UIPstream
(
    const DynamicList<char>& recvBuf,
    IOstreamOption::streamFormat fmt
)
:
    UIPstreamBase(recvBuf, fmt)
{}


Foam::IPstream::IPstream
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    const label bufSize,
    const int tag,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    Pstream(commsType, bufSize),
    UIPstream
    (
        commsType,
        fromProcNo,
        Pstream::transferBuf_,
        UIPstreamBase::storedRecvBufPos_,   // Internal only
        tag,
        comm,
        false,  // Do not clear Pstream::transferBuf_ if at end
        fmt
    )
{}


// ************************************************************************* //
