/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
#include "profilingPstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::UOPstream::bufferIPCsend()
{
    return UOPstream::write
    (
        commsType(),
        toProcNo_,
        sendBuf_.cdata(),
        sendBuf_.size(),
        tag_,
        comm_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPstream::write
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator,
    UPstream::Request* req,
    const UPstream::sendModes sendMode
)
{
    PstreamGlobals::reset_request(req);

    // TODO: some corrective action, at least when not nonBlocking
    #if 0
    if (bufSize > std::streamsize(INT_MAX))
    {
        Perr<< "UOPstream::write() : to rank " << toProcNo
            << " exceeds INT_MAX bytes" << Foam::endl;
        error::printStack(Perr);
    }
    #endif

    if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
    {
        Perr<< "UOPstream::write : starting write to:" << toProcNo
            << " size:" << label(bufSize)
            << " tag:" << tag << " comm:" << communicator
            << " commType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Perr);
    }
    else if (UPstream::debug)
    {
        Perr<< "UOPstream::write : starting write to:" << toProcNo
            << " size:" << label(bufSize)
            << " tag:" << tag << " comm:" << communicator
            << " commType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }

    PstreamGlobals::checkCommunicator(communicator, toProcNo);

    int returnCode = MPI_ERR_UNKNOWN;

    profilingPstream::beginTiming();

    if (commsType == UPstream::commsTypes::buffered)
    {
        returnCode = MPI_Bsend
        (
            const_cast<char*>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,
            tag,
            PstreamGlobals::MPICommunicators_[communicator]
        );

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (UPstream::debug)
        {
            Perr<< "UOPstream::write : finished buffered send to:"
                << toProcNo
                << " size:" << label(bufSize) << " tag:" << tag
                << Foam::endl;
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::sendModes::sync == sendMode)
        {
            returnCode = MPI_Ssend
            (
                const_cast<char*>(buf),
                bufSize,
                MPI_BYTE,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }
        else
        {
            returnCode = MPI_Send
            (
                const_cast<char*>(buf),
                bufSize,
                MPI_BYTE,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (UPstream::debug)
        {
            Perr<< "UOPstream::write : finished send to:"
                << toProcNo
                << " size:" << label(bufSize) << " tag:" << tag
                << Foam::endl;
        }
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        MPI_Request request;

        if (UPstream::sendModes::sync == sendMode)
        {
            returnCode = MPI_Issend
            (
                const_cast<char*>(buf),
                bufSize,
                MPI_BYTE,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );
        }
        else
        {
            returnCode = MPI_Isend
            (
                const_cast<char*>(buf),
                bufSize,
                MPI_BYTE,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );
        }

        if (UPstream::debug)
        {
            Perr<< "UOPstream::write : started non-blocking send to:"
                << toProcNo
                << " size:" << label(bufSize) << " tag:" << tag
                << " request:" <<
                (req ? label(-1) : PstreamGlobals::outstandingRequests_.size())
                << Foam::endl;
        }

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << Foam::abort(FatalError);
    }

    return (returnCode == MPI_SUCCESS);
}


// ************************************************************************* //
