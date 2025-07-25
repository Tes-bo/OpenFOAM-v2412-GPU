/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "IOstreams.H"

// FUTURE? probe and receive message
// - as of 2023-06 appears to be broken with INTELMPI + PMI-2 (slurm)
//   and perhaps other places so currently avoid

#undef Pstream_use_MPI_Get_count

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// General blocking/non-blocking MPI receive
static std::streamsize UPstream_mpi_receive
(
    const Foam::UPstream::commsTypes commsType,
    char* buf,
    const std::streamsize bufSize,
    const int fromProcNo,
    const int tag,
    const Foam::label communicator,
    Foam::UPstream::Request* req
)
{
    using namespace Foam;

    PstreamGlobals::reset_request(req);

    // TODO: some corrective action, at least when not nonBlocking
    #if 0
    // No warnings here, just on the sender side.
    if (bufSize > std::streamsize(INT_MAX))
    {
        Perr<< "UIPstream::read() : from rank " << fromProcNo
            << " exceeds INT_MAX bytes" << Foam::endl;
        error::printStack(Perr);
    }
    #endif

    if (UPstream::warnComm >= 0 && communicator != UPstream::warnComm)
    {
        Perr<< "UIPstream::read : starting read from:" << fromProcNo
            << " size:" << label(bufSize)
            << " tag:" << tag << " comm:" << communicator
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Perr);
    }
    else if (UPstream::debug)
    {
        Perr<< "UIPstream::read : starting read from:" << fromProcNo
            << " size:" << label(bufSize)
            << " tag:" << tag << " comm:" << communicator
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }

    int returnCode = MPI_ERR_UNKNOWN;

    profilingPstream::beginTiming();

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        // Not UPstream::commsTypes::nonBlocking

        MPI_Status status;

        {
            returnCode = MPI_Recv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            );
        }

        profilingPstream::addGatherTime();

        if (returnCode != MPI_SUCCESS)
        {
            FatalErrorInFunction
                << "MPI_Recv cannot receive incoming message"
                << Foam::abort(FatalError);
            return 0;
        }
        else if (UPstream::debug)
        {
            Perr<< "UIPstream::read : finished recv from:"
                << fromProcNo
                << " size:" << label(bufSize) << " tag:" << tag
                << Foam::endl;
        }

        // Check size of message read
        #ifdef Pstream_use_MPI_Get_count
        int count(0);
        MPI_Get_count(&status, MPI_BYTE, &count);
        #else
        MPI_Count count(0);
        MPI_Get_elements_x(&status, MPI_BYTE, &count);
        #endif

        // Errors
        if (count == MPI_UNDEFINED || int64_t(count) < 0)
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "returned undefined or negative value"
                << Foam::abort(FatalError);
        }
        else if (int64_t(count) > int64_t(UList<char>::max_size()))
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "count is larger than UList<char>::max_size() bytes"
                << Foam::abort(FatalError);
        }


        if (bufSize < std::streamsize(count))
        {
            FatalErrorInFunction
                << "buffer (" << label(bufSize)
                << ") not large enough for incoming message ("
                << label(count) << ')'
                << Foam::abort(FatalError);
        }

        return std::streamsize(count);
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        MPI_Request request;

        {
            returnCode = MPI_Irecv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &request
            );
        }

        if (returnCode != MPI_SUCCESS)
        {
            FatalErrorInFunction
                << "MPI_Irecv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();


        if (UPstream::debug)
        {
            Perr<< "UIPstream::read : started non-blocking recv from:"
                << fromProcNo
                << " size:" << label(bufSize) << " tag:" << tag
                << " request:" <<
                (req ? label(-1) : PstreamGlobals::outstandingRequests_.size())
                << Foam::endl;
        }

        // Assume the message will be completely received.
        return bufSize;
    }

    FatalErrorInFunction
        << "Unsupported communications type " << int(commsType)
        << Foam::abort(FatalError);

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UIPstream::bufferIPCrecv()
{
    // Called by constructor
    if (UPstream::debug)
    {
        Perr<< "UIPstream IPC read buffer :"
            << " from:" << fromProcNo_
            << " tag:" << tag_ << " comm:" << comm_
            << " wanted size:" << recvBuf_.capacity()
            << Foam::endl;
    }

    // Fallback value
    messageSize_ = recvBuf_.capacity();

    if (commsType() == UPstream::commsTypes::nonBlocking)
    {
        // Non-blocking
        // ~~~~~~~~~~~~
        // No chance of probing for size nor relying on the returned message
        // size (since it returns immediately without any further checks)
        //
        // Fortunately there are not many (any?) places that are using
        // a non-blocking IPstream with streaming anyhow.

        messageSize_ = recvBuf_.size();
    }
    else if (!recvBuf_.capacity())
    {
        // No buffer size allocated/specified - probe size of incoming message
        profilingPstream::beginTiming();

        MPI_Status status;

        MPI_Probe
        (
            fromProcNo_,
            tag_,
            PstreamGlobals::MPICommunicators_[comm_],
           &status
        );

        profilingPstream::addProbeTime();


        #ifdef Pstream_use_MPI_Get_count
        int count(0);
        MPI_Get_count(&status, MPI_BYTE, &count);
        #else
        MPI_Count count(0);
        MPI_Get_elements_x(&status, MPI_BYTE, &count);
        #endif

        // Errors
        if (count == MPI_UNDEFINED || int64_t(count) < 0)
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "returned undefined or negative value"
                << Foam::abort(FatalError);
        }
        else if (int64_t(count) > int64_t(UList<char>::max_size()))
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "count is larger than UList<char>::max_size() bytes"
                << Foam::abort(FatalError);
        }

        if (UPstream::debug)
        {
            Perr<< "UIPstream::UIPstream : probed size:"
                << label(count) << Foam::endl;
        }

        recvBuf_.resize(label(count));
        messageSize_ = label(count);
    }

    std::streamsize count = UPstream_mpi_receive
    (
        commsType(),
        recvBuf_.data(),
        messageSize_,   // The expected size
        fromProcNo_,
        tag_,
        comm_,
        nullptr   // UPstream::Request
    );

    if (count < 0)
    {
        FatalErrorInFunction
            << "MPI_recv() with negative size??"
            << Foam::abort(FatalError);
    }
    else if (int64_t(count) > int64_t(UList<char>::max_size()))
    {
        FatalErrorInFunction
            << "MPI_recv() larger than "
                "UList<char>::max_size() bytes"
            << Foam::abort(FatalError);
    }

    // Set addressed size. Leave actual allocated memory intact.
    recvBuf_.resize(label(count));
    messageSize_ = label(count);

    if (recvBuf_.empty())
    {
        setEof();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::streamsize Foam::UIPstream::read
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator,
    UPstream::Request* req
)
{
    return UPstream_mpi_receive
    (
        commsType,
        buf,
        bufSize,
        fromProcNo,
        tag,
        communicator,
        req
    );
}


// ************************************************************************* //
