/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Note
    The send/recv windows for chunk-wise transfers:

        iter    data window
        ----    -----------
        0       [0, 1*chunk]
        1       [1*chunk, 2*chunk]
        2       [2*chunk, 3*chunk]
        ...

    In older versions (v2312 and earlier) the number of chunks was
    determined by the sender sizes and used an extra MPI_Allreduce.
    However instead rely on the send/recv buffers having been consistently
    sized, which avoids the additional reduction.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "contiguous.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * Details * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

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
inline std::size_t maxTransferBytes(const int64_t max_bytes) noexcept
{
    return
    (
        (max_bytes < 0)  // (numBytes fewer than INT_MAX)
      ? std::size_t(INT_MAX + max_bytes)
      : std::size_t(max_bytes)
    );
}


//- Exchange of \em contiguous data, with or without chunking.
//- Setup sends and receives, each specified as [rank, span] tuple.
//  The serial list of tuples can be populated from other lists, from
//  maps of data or subsets of lists/maps etc.
//
//  Any waiting for requests is done outside by the caller
template<class Type>
void exchangeBuffers
(
    const UList<std::pair<int, stdFoam::span<const Type>>>& sends,
    const UList<std::pair<int, stdFoam::span<Type>>>& recvs,
    const int tag,
    const label comm,
    const int64_t maxComms_bytes = UPstream::maxCommsSize
)
{
    typedef stdFoam::span<const Type> sendType;
    typedef stdFoam::span<Type> recvType;

    // Caller already checked for parRun

    if (sends.empty() && recvs.empty())
    {
        // Nothing to do
        return;
    }

    const int myProci = UPstream::myProcNo(comm);


    // The chunk size (number of elements) corresponding to max byte transfer.
    // Is zero for non-chunked exchanges.
    const std::size_t chunkSize
    (
        PstreamDetail::maxTransferCount<Type>
        (
            PstreamDetail::maxTransferBytes(maxComms_bytes)
        )
    );


    // Set up receives
    // ~~~~~~~~~~~~~~~

    for (auto& slot : recvs)
    {
        // [rank, span]
        const auto proci = slot.first;
        auto& payload = slot.second;

        // No self-communication or zero-size payload
        if (proci == myProci || payload.empty())
        {
            continue;
        }
        else if (!chunkSize)
        {
            // Dispatch without chunks
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                payload.data_bytes(),
                payload.size_bytes(),
                tag,
                comm
            );
            continue;
        }

        // Dispatch chunk-wise until there is nothing left
        for (int iter = 0; /*true*/; ++iter)
        {
            // The begin/end for the data window
            const std::size_t beg = (std::size_t(iter)*chunkSize);
            const std::size_t end = (std::size_t(iter+1)*chunkSize);

            // The message tag augmented by the iteration number
            // - avoids message collisions between different chunks
            const int msgTagChunk = (tag + iter);

            if (payload.size() <= beg)
            {
                // No more data windows
                break;
            }

            recvType window
            (
                (end < payload.size())
              ? payload.subspan(beg, end - beg)
              : payload.subspan(beg)
            );

            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                window.data_bytes(),
                window.size_bytes(),
                msgTagChunk,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    for (const auto& slot : sends)
    {
        // [rank, span]
        const auto proci = slot.first;
        const auto& payload = slot.second;

        // No self-communication or zero-size payload
        if (proci == myProci || payload.empty())
        {
            continue;
        }
        else if (!chunkSize)
        {
            // Dispatch without chunks
            bool ok = UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                payload.cdata_bytes(),
                payload.size_bytes(),
                tag,
                comm
            );

            if (!ok)
            {
                FatalErrorInFunction
                    << "Failure sending message to:" << proci
                    << " nBytes:" << label(payload.size_bytes()) << nl
                    << Foam::abort(FatalError);
            }
            continue;
        }

        // Dispatch chunk-wise until there is nothing left
        for (int iter = 0; /*true*/; ++iter)
        {
            // The begin/end for the data window
            const std::size_t beg = (std::size_t(iter)*chunkSize);
            const std::size_t end = (std::size_t(iter+1)*chunkSize);

            // The message tag augmented by the iteration number
            // - avoids message collisions between different chunks
            const int msgTagChunk = (tag + iter);

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

            bool ok = UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                window.cdata_bytes(),
                window.size_bytes(),
                msgTagChunk,
                comm
            );

            if (!ok)
            {
                FatalErrorInFunction
                    << "Failure sending message to:" << proci
                    << " nBytes:" << label(window.size_bytes()) << nl
                    << Foam::abort(FatalError);
            }
        }
    }
}


//- Exchange \em contiguous data using point-to-point communication.
//- Sends sendBufs, receives into recvBufs.
//  Data provided and received as container all of which have been
//  properly sized before calling
//
// No internal guards or resizing.
template<class Container, class Type>
void exchangeContainer
(
    const UList<Container>& sendBufs,
    UList<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait,        //!< Wait for requests to complete
    const int64_t maxComms_bytes = UPstream::maxCommsSize
)
{
    typedef stdFoam::span<const Type> sendType;
    typedef stdFoam::span<Type> recvType;

    // Caller already checked for parRun

    if (sendBufs.empty() && recvBufs.empty())
    {
        // Nothing to do
        return;
    }

    const label startOfRequests = UPstream::nRequests();
    const int myProci = UPstream::myProcNo(comm);

    // The chunk size (number of elements) corresponding to max byte transfer
    const std::size_t chunkSize
    (
        PstreamDetail::maxTransferCount<Type>
        (
            PstreamDetail::maxTransferBytes(maxComms_bytes)
        )
    );


    // Set up receives
    // ~~~~~~~~~~~~~~~

    forAll(recvBufs, proci)
    {
        // [rank, span]
        recvType payload
        (
            recvBufs[proci].data(),
            std::size_t(recvBufs[proci].size())
        );

        // No self-communication or zero-size payload
        if (proci == myProci || payload.empty())
        {
            continue;
        }
        else if (!chunkSize)
        {
            // Dispatch without chunks
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                payload.data_bytes(),
                payload.size_bytes(),
                tag,
                comm
            );
            continue;
        }

        // Dispatch chunk-wise until there is nothing left
        for (int iter = 0; /*true*/; ++iter)
        {
            // The begin/end for the data window
            const std::size_t beg = (std::size_t(iter)*chunkSize);
            const std::size_t end = (std::size_t(iter+1)*chunkSize);

            // The message tag augmented by the iteration number
            // - avoids message collisions between different chunks
            const int msgTagChunk = (tag + iter);

            if (payload.size() <= beg)
            {
                // No more data windows
                break;
            }

            recvType window
            (
                (end < payload.size())
              ? payload.subspan(beg, end - beg)
              : payload.subspan(beg)
            );

            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                window.data_bytes(),
                window.size_bytes(),
                msgTagChunk,
                comm
            );
        }
    }


    // Set up sends
    // ~~~~~~~~~~~~

    forAll(sendBufs, proci)
    {
        // [rank, span]
        sendType payload
        (
            sendBufs[proci].cdata(),
            std::size_t(sendBufs[proci].size())
        );

        // No self-communication or zero-size payload
        if (proci == myProci || payload.empty())
        {
            continue;
        }
        else if (!chunkSize)
        {
            // Dispatch without chunks
            bool ok = UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                payload.cdata_bytes(),
                payload.size_bytes(),
                tag,
                comm
            );

            if (!ok)
            {
                FatalErrorInFunction
                    << "Fallure sending message to:" << proci
                    << " nBytes:" << label(payload.size_bytes()) << nl
                    << Foam::abort(FatalError);
            }
            continue;
        }

        // Dispatch chunk-wise until there is nothing left
        for (int iter = 0; /*true*/; ++iter)
        {
            // The begin/end for the data window
            const std::size_t beg = (std::size_t(iter)*chunkSize);
            const std::size_t end = (std::size_t(iter+1)*chunkSize);

            // The message tag augmented by the iteration number
            // - avoids message collisions between different chunks
            const int msgTagChunk = (tag + iter);

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

            bool ok = UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                window.cdata_bytes(),
                window.size_bytes(),
                msgTagChunk,
                comm
            );

            if (!ok)
            {
                FatalErrorInFunction
                    << "Failure sending message to:" << proci
                    << " nBytes:" << label(window.size_bytes()) << nl
                    << Foam::abort(FatalError);
            }
        }
    }


    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (UPstream::debug)
    {
        Perr<< "Pstream::exchange with "
            << (UPstream::nRequests() - startOfRequests)
            << " requests" << nl;
    }

    if (wait)
    {
        UPstream::waitRequests(startOfRequests);
    }
}


//- Exchange \em contiguous data using point-to-point communication.
//- Sends sendBufs, receives into recvBufs.
//  Data provided and received as container all of which have been
//  properly sized before calling
//
// No internal guards or resizing.
template<class Container, class Type>
void exchangeContainer
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait,        //!< Wait for requests to complete
    const int64_t maxComms_bytes = UPstream::maxCommsSize
)
{
    typedef stdFoam::span<const Type> sendType;
    typedef stdFoam::span<Type> recvType;

    typedef std::pair<int, sendType> sendTuple;
    typedef std::pair<int, recvType> recvTuple;

    const label startOfRequests = UPstream::nRequests();
    const label myProci = UPstream::myProcNo(comm);

    // Serialize recv sequences
    DynamicList<recvTuple> recvs(recvBufs.size());
    {
        forAllIters(recvBufs, iter)
        {
            const auto proci = iter.key();
            auto& recvData = recvBufs[proci];

            recvType payload
            (
                recvData.data(),
                std::size_t(recvData.size())
            );

            // No self-communication or zero-size payload
            if (proci != myProci && !payload.empty())
            {
                recvs.emplace_back(proci, payload);
            }
        }

        std::sort
        (
            recvs.begin(),
            recvs.end(),
            [=](const recvTuple& a, const recvTuple& b)
            {
                // Sorted processor order
                return (a.first < b.first);

                // OR: // Shorter messages first
                // return (a.second.size() < b.second.size());
            }
        );
    }

    // Serialize send sequences
    DynamicList<sendTuple> sends(sendBufs.size());
    {
        forAllConstIters(sendBufs, iter)
        {
            const auto proci = iter.key();
            const auto& sendData = iter.val();

            sendType payload
            (
                sendData.cdata(),
                std::size_t(sendData.size())
            );

            // No self-communication or zero-size payload
            if (proci != myProci && !payload.empty())
            {
                sends.emplace_back(proci, payload);
            }
        }

        std::sort
        (
            sends.begin(),
            sends.end(),
            [=](const sendTuple& a, const sendTuple& b)
            {
                // Sorted processor order
                return (a.first < b.first);

                // OR: // Shorter messages first
                // return (a.second.size() < b.second.size());
            }
        );
    }

    // Exchange buffers in chunk-wise transfers
    PstreamDetail::exchangeBuffers<Type>
    (
        sends,
        recvs,
        tag,
        comm,
        maxComms_bytes
    );

    // Wait for all to finish
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (UPstream::debug)
    {
        Perr<< "Pstream::exchange with "
            << (UPstream::nRequests() - startOfRequests)
            << " requests" << nl;
    }

    if (wait)
    {
        UPstream::waitRequests(startOfRequests);
    }
}

} // namespace PstreamDetail
} // namespace Foam

#include "PstreamExchangeConsensus.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class Type>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(comm);
    const label numProc = UPstream::nProcs(comm);

    if (sendBufs.size() != numProc)
    {
        FatalErrorInFunction
            << "List size " << sendBufs.size()
            << " != number of ranks " << numProc << nl
            << Foam::abort(FatalError);
    }

    recvBufs.resize_nocopy(numProc);

    if (UPstream::is_parallel(comm))
    {
        // Presize all receive buffers
        forAll(recvSizes, proci)
        {
            const label count = recvSizes[proci];

            if (proci != myProci && count > 0)
            {
                recvBufs[proci].resize_nocopy(count);
            }
            else
            {
                recvBufs[proci].clear();
            }
        }

        PstreamDetail::exchangeContainer<Container, Type>
        (
            sendBufs,
            recvBufs,
            tag,
            comm,
            wait
            // (default: UPstream::maxCommsSize)
        );
    }

    // Do myself. Already checked if in communicator
    recvBufs[myProci] = sendBufs[myProci];
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const Map<Container>& sendBufs,
    const Map<label>& recvSizes,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    static_assert(is_contiguous<Type>::value, "Contiguous data only!");

    const int myProci = UPstream::myProcNo(comm);

    // Initial: clear out receive 'slots'
    // Preferrable to clear out the map entries instead of the map itself
    // since this can potentially preserve allocated space
    // (eg DynamicList entries) between calls

    forAllIters(recvBufs, iter)
    {
        iter.val().clear();
    }

    if (UPstream::is_parallel(comm))
    {
        // Presize all receive buffers
        forAllIters(recvSizes, iter)
        {
            const label proci = iter.key();
            const label count = iter.val();

            if (proci != myProci && count > 0)
            {
                recvBufs(proci).resize_nocopy(count);
            }
        }

        PstreamDetail::exchangeContainer<Container, Type>
        (
            sendBufs,
            recvBufs,
            tag,
            comm,
            wait
            // (default: UPstream::maxCommsSize)
        );
    }

    // Do myself (if actually in the communicator)
    if (UPstream::is_rank(comm))
    {
        const auto iter = sendBufs.find(myProci);

        bool needsCopy = iter.good();

        if (needsCopy)
        {
            const auto& sendData = iter.val();

            needsCopy = !sendData.empty();
            if (needsCopy)
            {
                // insert_or_assign
                recvBufs(myProci) = sendData;
            }
        }

        if (!needsCopy)
        {
            recvBufs.erase(myProci);
        }
    }
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const Container& sendBufs,
    labelList& recvSizes,
    const label tag,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(comm);
    const label numProc = UPstream::nProcs(comm);

    if (sendBufs.size() != numProc)
    {
        FatalErrorInFunction
            << "Container size " << sendBufs.size()
            << " != number of ranks " << numProc << nl
            << Foam::abort(FatalError);
    }

    labelList sendSizes(numProc);
    for (label proci = 0; proci < numProc; ++proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }

    recvSizes.resize_nocopy(numProc);
    recvSizes = 0;  // Ensure non-received entries are properly zeroed

    // Preserve self-send, even if not described by neighbourhood
    recvSizes[myProci] = sendSizes[myProci];

    const label startOfRequests = UPstream::nRequests();

    for (const label proci : recvProcs)
    {
        if (proci != myProci)
        {
            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                reinterpret_cast<char*>(&recvSizes[proci]),
                sizeof(label),
                tag,
                comm
            );
        }
    }

    for (const label proci : sendProcs)
    {
        if (proci != myProci)
        {
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                reinterpret_cast<char*>(&sendSizes[proci]),
                sizeof(label),
                tag,
                comm
            );
        }
    }

    UPstream::waitRequests(startOfRequests);
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const labelUList& neighProcs,
    const Container& sendBufs,
    labelList& recvSizes,
    const label tag,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    Pstream::exchangeSizes<Container>
    (
        neighProcs,  // send
        neighProcs,  // recv
        sendBufs,
        recvSizes,
        tag,
        comm
    );
}


// Sparse sending
template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Map<Container>& sendBufs,
    Map<label>& recvSizes,
    const label tag,
    const label comm
)
{
    recvSizes.clear();  // Done in allToAllConsensus too, but be explicit here

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    Map<label> sendSizes;
    sendSizes.reserve(sendBufs.size());

    forAllConstIters(sendBufs, iter)
    {
        const label proci = iter.key();
        const label count = iter.val().size();

        if (count)
        {
            sendSizes.emplace(proci, count);
        }
    }

    UPstream::allToAllConsensus
    (
        sendSizes,
        recvSizes,
        (tag + 314159),  // some unique tag?
        comm
    );
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Container& sendBufs,
    labelList& recvSizes,
    const label comm
)
{
    if (!UPstream::is_rank(comm))
    {
        recvSizes.clear();
        return;  // Process not in communicator
    }

    const label numProc = UPstream::nProcs(comm);

    if (sendBufs.size() != numProc)
    {
        FatalErrorInFunction
            << "Container size " << sendBufs.size()
            << " != number of ranks " << numProc << nl
            << Foam::abort(FatalError);
    }

    labelList sendSizes(numProc);
    forAll(sendBufs, proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }
    recvSizes.resize_nocopy(sendSizes.size());

    if
    (
        UPstream::nProcsNonblockingExchange > 1
     && UPstream::nProcsNonblockingExchange <= numProc
    )
    {
        // Use algorithm NBX: Nonblocking Consensus Exchange

        UPstream::allToAllConsensus
        (
            sendSizes,
            recvSizes,
            (UPstream::msgType() + 314159),  // some unique tag?
            comm
        );
        return;
    }

    UPstream::allToAll(sendSizes, recvSizes, comm);
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    // Algorithm PEX: Personalized Exchange
    // - Step 1: each process writes the data sizes to each peer and
    //   redistributes the vector (eg, MPI_Alltoall or non-blocking consensus)
    // - Step 2: size receive buffers and setup receives for all
    //   non-zero sendcounts. Post all sends and wait.

    labelList recvSizes;
    Pstream::exchangeSizes(sendBufs, recvSizes, comm);

    Pstream::exchange<Container, Type>
    (
        sendBufs,
        recvSizes,
        recvBufs,
        tag,
        comm,
        wait
    );
}


template<class Container, class Type>
void Foam::Pstream::exchange
(
    const Map<Container>& sendBufs,
    Map<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool wait
)
{
    // Algorithm PEX: Personalized Exchange
    // but using nonblocking consensus exchange for the sizes

    Map<label> recvSizes;
    Pstream::exchangeSizes(sendBufs, recvSizes, tag, comm);

    Pstream::exchange<Container, Type>
    (
        sendBufs,
        recvSizes,
        recvBufs,
        tag,
        comm,
        wait
    );
}


// ************************************************************************* //
