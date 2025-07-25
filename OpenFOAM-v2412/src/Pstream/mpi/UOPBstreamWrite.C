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

#include "UOPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::UOPBstream::bufferIPCsend()
{
    // Uses double broadcast
    // 1. for the data size
    // 2. for the data itself
    // With this information, can determine and resize receive buffer

    PstreamGlobals::checkCommunicator(comm_, toProcNo_);

    // Same type must be expected in UIPBstream::bufferIPCrecv()
    std::streamsize bufSize(sendBuf_.size());

    // TODO: some corrective action
    #if 0
    if (bufSize > std::streamsize(INT_MAX))
    {
        Perr<< "UOPBstream::write() :"
            << " exceeds INT_MAX bytes" << Foam::endl;
        error::printStack(Perr);
    }
    #endif

    // Broadcast #1 - data size
    if
    (
        !UPstream::broadcast
        (
            reinterpret_cast<char*>(&bufSize),
            sizeof(std::streamsize),
            comm_,
            toProcNo_  //< is actually rootProcNo
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Bcast failure sending buffer size:" << bufSize << nl
            << Foam::abort(FatalError);
        return false;
    }


    // Broadcast #2 - data content
    // - skip if there is no data to send
    if
    (
        (bufSize > 0)
     && !UPstream::broadcast
        (
            sendBuf_.data(),
            sendBuf_.size(),  // same as bufSize
            comm_,
            toProcNo_  //< is actually rootProcNo
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Bcast failure sending buffer data:"
            << sendBuf_.size() << nl
            << Foam::abort(FatalError);
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPBstream::write
(
    const int rootProcNo,
    const char* buf,
    const std::streamsize bufSize,
    const label comm
)
{
    if
    (
        !UPstream::broadcast(const_cast<char*>(buf), bufSize, comm, rootProcNo)
    )
    {
        FatalErrorInFunction
            << "MPI_Bcast failure sending buffer data:" << label(bufSize) << nl
            << Foam::abort(FatalError);
        return false;
    }

    return true;
}


// ************************************************************************* //
