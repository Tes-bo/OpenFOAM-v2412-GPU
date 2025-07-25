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

Application
    Test-parallel-chunks

Description
    Test for sending contiguous data in chunk-wise.
    Largely mirrors Pstream::exchange or vice versa

\*---------------------------------------------------------------------------*/

#define Foam_PstreamExchange_debug_chunks

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

using namespace Foam;


//- Number of elements corresponding to max byte transfer.
//  Normal upper limit is INT_MAX since MPI sizes are limited to <int>.
template<class Type>
inline std::size_t maxTransferCount
(
    const std::size_t max_bytes = std::size_t(0)
) noexcept
{
    return
    (
        (max_bytes == 0)                        // ie, unlimited
      ? (std::size_t(0))                        //
      : (max_bytes > std::size_t(INT_MAX))      // MPI limit is <int>
      ? (std::size_t(INT_MAX) / sizeof(Type))   //
      : (max_bytes > sizeof(Type))              // require an integral number
      ? (max_bytes / sizeof(Type))              //
      : (std::size_t(1))                        // min of one element
    );
}


//- Upper limit on number of transfer bytes.
//  Max bytes is normally INT_MAX since MPI sizes are limited to <int>.
//  Negative values indicate a subtraction from INT_MAX.
inline std::size_t PstreamDetail_maxTransferBytes
(
    const int64_t max_bytes
) noexcept
{
    return
    (
        (max_bytes < 0)  // (numBytes fewer than INT_MAX)
      ? std::size_t(INT_MAX + max_bytes)
      : std::size_t(max_bytes)
    );
}


template<class Container, class Type>
void broadcast_chunks
(
    Container& sendData,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
    const int64_t maxComms_bytes = UPstream::maxCommsSize
)
{
    // OR  static_assert(is_contiguous<T>::value, "Contiguous data only!")
    if (!is_contiguous<Type>::value)
    {
        FatalErrorInFunction
            << "Contiguous data only." << sizeof(Type)
            << Foam::abort(FatalError);
    }

    if (maxComms_bytes == 0)
    {
        // Do in one go
        Info<< "send " << sendData.size() << " elements in one go" << endl;
        Pstream::broadcast(sendData, comm);
        return;
    }

    label sendSize(sendData.size());
    Pstream::broadcast(sendSize, comm);

    label recvSize(sendSize);

    sendData.resize_nocopy(recvSize);  // A no-op on master


    // The chunk size (number of elements) corresponding to max byte transfer
    // Is zero for non-chunked exchanges.
    const std::size_t chunkSize
    (
        PstreamDetail_maxTransferCount<Type>
        (
            PstreamDetail_maxTransferBytes(maxComms_bytes)
        )
    );


    if (chunkSize)
    {
        // Convert from send count (elements) to number of chunks.
        // Can normally calculate with (count-1), but add some safety
        label nChunks = 1 + (sendSize/label(chunkSize));

        Info
            << "send " << sendSize << " elements ("
            << (sendSize*sizeof(Type)) << " bytes) in " << nChunks
            << " chunks of " << label(chunkSize) << " elements ("
            << label(chunkSize*sizeof(Type)) << " bytes) for maxCommsSize:"
            << label(maxComms_bytes)
            << endl;
    }


    // stress-test with shortened sendSize
    // will produce useless loops, but no calls
    // sendSize /= 2;

    typedef stdFoam::span<Type> sendType;

    do
    {
        sendType payload(sendData.data(), sendData.size());

        if (!chunkSize)
        {
            UPstream::broadcast
            (
                payload.data_bytes(),
                payload.size_bytes(),
                comm
            );
            break;
        }

        // Dispatch chunk-wise until there is nothing left
        for (int iter = 0; /*true*/; ++iter)
        {
            // The begin/end for the data window
            const std::size_t beg = (std::size_t(iter)*chunkSize);
            const std::size_t end = (std::size_t(iter+1)*chunkSize);

            if (payload.size() <= beg)
            {
                // No more data windows
                break;
            }

            sendType window
            (
                (end < payload.size())
              ? payload.subspan(beg, end - beg)
              : payload.subspan(beg)
            );

            Info<< "iter " << iter
                << ": beg=" << label(beg) << " len=" << label(window.size())
                << " (" << label(window.size_bytes()) << " bytes)" << endl;

            UPstream::broadcast
            (
                window.data_bytes(),
                window.size_bytes(),
                comm
            );
        }
    }
    while (false);

    Info<< "final" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addOption("comms-size", "int", "override Pstream::maxCommsSize");

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    labelList input1;
    if (UPstream::master())
    {
        input1 = identity(500);
    }
    broadcast_chunks<labelList, label>(input1);

    Pstream::maxCommsSize = 33;

    args.readIfPresent("comms-size", Pstream::maxCommsSize);

    broadcast_chunks<labelList, label>(input1);

    // Mostly the same with PstreamBuffers
    if (false)
    {
        PstreamBuffers pBufs;

        labelList sendData;
        if (Pstream::master())
        {
            sendData = identity(500);

            for (const int proci : Pstream::subProcs())
            {
                UOPstream os(proci, pBufs);
                os << sendData;
            }
        }

        Info<< "call finishedSends()" << endl;
        pBufs.finishedScatters();

        if (!Pstream::master())
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> sendData;
        }
    }

    // Manually
    Info<< "perform list exchange" << endl;
    {
        labelListList sendBufs(UPstream::nProcs());
        labelListList recvBufs(UPstream::nProcs());
        labelList recvSizes;

        if (Pstream::master())
        {
            for (const int proci : Pstream::allProcs())
            {
                if (proci != Pstream::myProcNo())
                {
                    sendBufs[proci] = identity(500);
                }
            }
        }

        Pstream::exchangeSizes(sendBufs, recvSizes);

        Pstream::exchange<labelList, label>
        (
            sendBufs,
            recvSizes,
            recvBufs
        );
    }


    Info<< "perform Map exchange" << endl;
    {
        Map<labelList> sendBufs;
        Map<labelList> recvBufs;
        Map<label> recvSizes;

        if (Pstream::master())
        {
            for (const int proci : Pstream::allProcs())
            {
                if (proci != Pstream::myProcNo())
                {
                    sendBufs(proci) = identity(500);
                }
            }
        }

        Pstream::exchangeSizes(sendBufs, recvSizes);

        Pstream::exchange<labelList, label>
        (
            sendBufs,
            recvSizes,
            recvBufs
        );
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
