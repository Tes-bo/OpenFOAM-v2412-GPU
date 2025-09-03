/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPstream::UOPstream
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const label comm,
    const bool sendAtDestruct,
    IOstreamOption::streamFormat fmt
)
:
    UOPstreamBase(commsType, toProcNo, sendBuf, tag, comm, sendAtDestruct, fmt)
{}


Foam::UOPstream::UOPstream(const int toProcNo, PstreamBuffers& buffers)
:
    UOPstreamBase(toProcNo, buffers)
{}


Foam::UOPstream::UOPstream
(
    DynamicList<char>& sendBuf,
    IOstreamOption::streamFormat fmt
)
:
    UOPstreamBase(sendBuf, fmt)
{}


Foam::OPstream::OPstream
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const label bufSize,
    const int tag,
    const label comm,
    IOstreamOption::streamFormat fmt
)
:
    Pstream(commsType, bufSize),
    UOPstream
    (
        commsType,
        toProcNo,
        Pstream::transferBuf_,
        tag,
        comm,
        true,  // sendAtDestruct
        fmt
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPstream::send()
{
    sendAtDestruct_ = false;
    return bufferIPCsend();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPstream::~UOPstream()
{
    if (sendAtDestruct_)
    {
        // Note: sendAtDestruct_ and nonBlocking is a questionable combination
        // since the transfer buffer will be destroyed before
        // the non-blocking send completes!
        //
        // Could flag as an error, but not actually used anywhere.

        if (!bufferIPCsend())
        {
            FatalErrorInFunction
                << "Failed sending outgoing message of size "
                << sendBuf_.size() << " to processor " << toProcNo_
                << Foam::abort(FatalError);
        }
    }
}


// ************************************************************************* //
