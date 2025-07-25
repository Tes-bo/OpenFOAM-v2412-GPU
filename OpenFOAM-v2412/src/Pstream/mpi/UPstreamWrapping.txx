/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "UPstreamWrapping.H"
#include "profilingPstream.H"
#include "PstreamGlobals.H"
#include "Map.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PstreamDetail::broadcast0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const label comm
)
{
    if (!UPstream::is_parallel(comm))
    {
        return;
    }

    profilingPstream::beginTiming();

    // const int returnCode =
    MPI_Bcast
    (
        values,
        count,
        datatype,
        0,  // (root rank) == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[comm]
    );

    profilingPstream::addBroadcastTime();
}


template<class Type>
void Foam::PstreamDetail::reduce0
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label comm
)
{
    if (!UPstream::is_parallel(comm))
    {
        return;
    }

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        Perr<< "** MPI_Reduce (blocking):";
        if (count == 1)
        {
            Perr<< (*values);
        }
        else
        {
            Perr<< UList<Type>(values, count);
        }
        Perr<< " with comm:" << comm
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Perr);
    }

    profilingPstream::beginTiming();

    // const int returnCode =
    MPI_Reduce
    (
        MPI_IN_PLACE,  // recv is also send
        values,
        count,
        datatype,
        optype,
        0,  // (root rank) == UPstream::masterNo()
        PstreamGlobals::MPICommunicators_[comm]
    );

    profilingPstream::addReduceTime();
}


template<class Type>
void Foam::PstreamDetail::allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_parallel(comm))
    {
        return;
    }

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Iallreduce (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Allreduce (blocking):";
        }
        if (count == 1)
        {
            Perr<< (*values);
        }
        else
        {
            Perr<< UList<Type>(values, count);
        }
        Perr<< " with comm:" << comm
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Perr);
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Iallreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iallreduce failed for "
                << UList<Type>(values, count)
                << Foam::abort(FatalError);
        }


        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Allreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Allreduce failed for "
                << UList<Type>(values, count)
                << Foam::abort(FatalError);
        }

        profilingPstream::addReduceTime();
    }
}


template<class Type>
void Foam::PstreamDetail::allToAll
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm))
    {
        return;
    }

    const label numProc = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Ialltoall (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Alltoall (blocking):";
        }
        Perr<< " numProc:" << numProc
            << " sendData:" << sendData.size()
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        (sendData.size() != numProc || recvData.size() != numProc)
    )
    {
        FatalErrorInFunction
            << "Have " << numProc << " ranks, but size of sendData:"
            << sendData.size() << " or recvData:" << recvData.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::is_parallel(comm))
    {
        recvData.deepCopy(sendData);
        return;
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Ialltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ialltoall [comm: " << comm << "] failed."
                << " For " << sendData
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall [comm: " << comm << "] failed."
                << " For " << sendData
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }
}


template<class Type>
void Foam::PstreamDetail::allToAllv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm))
    {
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Ialltoallv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Alltoallv (blocking):";
        }
        Perr<< " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        (sendCounts.size() != np || sendOffsets.size() < np)
     || (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << ", sendOffsets:" << sendOffsets.size()
            << ", recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::is_parallel(comm))
    {
        if (recvCounts[0] != sendCounts[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendCounts[0]
                << " does not equal bytes to receive " << recvCounts[0]
                << Foam::abort(FatalError);
        }
        std::memmove
        (
            recvData,
            (sendData + sendOffsets[0]),
            recvCounts[0]*sizeof(Type)
        );

        return;
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Ialltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ialltoallv [comm: " << comm << "] failed."
                << " For sendCounts " << sendCounts
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Alltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv [comm: " << comm << "] failed."
                << " For sendCounts " << sendCounts
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }

}


template<class Type>
void Foam::PstreamDetail::allToAllConsensus
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const int tag,
    const label comm
)
{
    const bool initialBarrier = (UPstream::tuning_NBX_ > 0);

    if (!UPstream::is_rank(comm))
    {
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(comm);
    const label numProc = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        Perr<< "** non-blocking consensus Alltoall (list):";
        Perr<< " numProc:" << numProc
            << " sendData:" << sendData.size()
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if (sendData.size() != numProc || recvData.size() != numProc)
    {
        FatalErrorInFunction
            << "Have " << numProc << " ranks, but size of sendData:"
            << sendData.size() << " or recvData:" << recvData.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    // Initial: assign zero everywhere. Values of zero are never transmitted
    const Type zeroValue = pTraits<Type>::zero;
    recvData = zeroValue;

    if (!UPstream::is_parallel(comm))
    {
        // Non-parallel : deep copy
        recvData.deepCopy(sendData);
        return;
    }

    // Fake send/recv for myself
    {
        recvData[myProci] = sendData[myProci];
    }

    // Implementation description
    // --------------------------
    // "Scalable Communication Protocols for Dynamic Sparse Data Exchange",
    // Hoeffler, Siebert, Lumsdaine
    // May 2010 ACM SIGPLAN Notices 45(5):159-168
    // https://doi.org/10.1145/1837853.1693476
    //
    // - http://unixer.de/publications/img/hoefler-dsde-protocols.pdf
    //
    // Algorithm NBX: Nonblocking consensus

    // This specific specialization is largely just for integer data
    // so we initialise the receiving data with zero and then
    // do not send/recv them.
    // This is because we are dealing with a flat list of entries to
    // send and not a sparse Map etc.


    // ------------------------------------------------------------------------
    // Setup sends
    // ------------------------------------------------------------------------

    profilingPstream::beginTiming();

    // An initial barrier may help to avoid synchronisation problems
    // caused elsewhere
    if (initialBarrier)
    {
        MPI_Barrier(PstreamGlobals::MPICommunicators_[comm]);
    }

    DynamicList<MPI_Request> sendRequests(sendData.size());

    // Start nonblocking synchronous send to destination rank
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (sendData[proci] == zeroValue)
        {
            // Do not send/recv empty data
        }
        else if (proci != myProci)
        {
            // Has data to send

            MPI_Issend
            (
               &sendData[proci],
                1,              // one element per rank
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
               &sendRequests.emplace_back()
            );
        }
    }


    // ------------------------------------------------------------------------
    // Probe and receive
    // ------------------------------------------------------------------------

    MPI_Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        int flag = 0;
        MPI_Status status;

        MPI_Iprobe
        (
            MPI_ANY_SOURCE,
            tag,
            PstreamGlobals::MPICommunicators_[comm],
           &flag,
           &status
        );

        if (flag)
        {
            // Message found, receive into dest buffer location
            const label proci = status.MPI_SOURCE;

            // Only send/recv a single (fundamental) data type
            int count(0);
            MPI_Get_count(&status, datatype, &count);

            if (count != 1)
            {
                FatalErrorInFunction
                    << "Incorrect message size from proc=" << proci
                    << ". Expected 1 but had " << count << nl
                    << exit(FatalError);
            }

            // Regular blocking receive [the data are small]

            MPI_Recv
            (
               &recvData[proci],
                count,          // count=1 (see above)
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
                MPI_STATUS_IGNORE
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);

            if (flag)
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            MPI_Testall
            (
                sendRequests.size(),
                sendRequests.data(),
               &flag, MPI_STATUSES_IGNORE
            );

            if (flag)
            {
                MPI_Ibarrier
                (
                    PstreamGlobals::MPICommunicators_[comm],
                   &barrierRequest
                );
                barrier_active = true;
            }
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::allToAllConsensus
(
    const Map<Type>& sendBufs,
    Map<Type>& recvBufs,
    MPI_Datatype datatype,
    const int tag,
    const label comm
)
{
    const bool initialBarrier = (UPstream::tuning_NBX_ > 0);

    const label myProci = UPstream::myProcNo(comm);
    const label numProc = UPstream::nProcs(comm);

    if (!UPstream::is_rank(comm))
    {
        return;
    }

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        Perr<< "** non-blocking consensus Alltoall (map):";
        Perr<< " numProc:" << numProc
            << " sendData:" << sendBufs.size()
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    // Initial: clear out everything
    recvBufs.clear();

    // Fake send/recv for myself - parallel or non-parallel
    {
        const auto iter = sendBufs.find(myProci);
        if (iter.good())
        {
            // Do myself: insert_or_assign
            recvBufs(iter.key()) = iter.val();
        }
    }

    if (!UPstream::is_parallel(comm))
    {
        // Nothing left to do
        return;
    }


    // ------------------------------------------------------------------------
    // Setup sends
    // ------------------------------------------------------------------------

    // Algorithm NBX: Nonblocking consensus
    // Implementation like above, but sending map data.

    DynamicList<MPI_Request> sendRequests(sendBufs.size());

    profilingPstream::beginTiming();

    // An initial barrier may help to avoid synchronisation problems
    // caused elsewhere
    if (initialBarrier)
    {
        MPI_Barrier(PstreamGlobals::MPICommunicators_[comm]);
    }


    // Start nonblocking synchronous send to destination ranks

    // Same as forAllConstIters()
    const auto endIter = sendBufs.cend();
    for (auto iter = sendBufs.cbegin(); iter != endIter; ++iter)
    {
        const label proci = iter.key();
        const auto& sendData = iter.val();

        if (proci != myProci && proci >= 0 && proci < numProc)
        {
            MPI_Issend
            (
               &sendData,
                1,              // one element per rank
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
               &sendRequests.emplace_back()
            );
        }
    }


    // ------------------------------------------------------------------------
    // Probe and receive
    // ------------------------------------------------------------------------

    MPI_Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        int flag = 0;
        MPI_Status status;

        MPI_Iprobe
        (
            MPI_ANY_SOURCE,
            tag,
            PstreamGlobals::MPICommunicators_[comm],
           &flag,
           &status
        );

        if (flag)
        {
            // Message found, receive into dest buffer location
            const label proci = status.MPI_SOURCE;

            // Only send/recv a single (fundamental) data type
            int count(0);
            MPI_Get_count(&status, datatype, &count);

            if (count != 1)
            {
                FatalErrorInFunction
                    << "Incorrect message size from proc=" << proci
                    << ". Expected 1 but had " << count << nl
                    << exit(FatalError);
            }

            auto& recvData = recvBufs(proci);

            // Regular blocking receive [the data are small]

            MPI_Recv
            (
               &recvData,
                count,          // count=1 (see above)
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[comm],
                MPI_STATUS_IGNORE
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);

            if (flag)
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            MPI_Testall
            (
                sendRequests.size(),
                sendRequests.data(),
               &flag, MPI_STATUSES_IGNORE
            );

            if (flag)
            {
                MPI_Ibarrier
                (
                    PstreamGlobals::MPICommunicators_[comm],
                   &barrierRequest
                );
                barrier_active = true;
            }
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::gather
(
    const Type* sendData,
    Type* recvData,

    int count,
    MPI_Datatype datatype,

    const label comm,
    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm) || !count)
    {
        return;
    }
    if (!UPstream::is_parallel(comm))
    {
        if (recvData)
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
        return;
    }

    const label numProc = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Igather (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Gather (blocking):";
        }
        Perr<< " numProc:" << numProc
            << " count:" << count
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Igather
            (
                const_cast<Type*>(sendData), count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Igather [comm: " << comm << "] failed."
                << " count:" << count << nl
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Gather
            (
                const_cast<Type*>(sendData), count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gather [comm: " << comm << "] failed."
                << " count:" << count << nl
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


template<class Type>
void Foam::PstreamDetail::scatter
(
    const Type* sendData,
    Type* recvData,

    int count,
    MPI_Datatype datatype,

    const label comm,
    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm) || !count)
    {
        return;
    }
    if (!UPstream::is_parallel(comm))
    {
        if (recvData)
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
        return;
    }

    const label numProc = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Iscatter (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Scatter (blocking):";
        }
        Perr<< " numProc:" << numProc
            << " count:" << count
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Iscatter
            (
                const_cast<Type*>(sendData), count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iscatter [comm: " << comm << "] failed."
                << " count:" << count << nl
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Scatter
            (
                const_cast<Type*>(sendData), count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatter [comm: " << comm << "] failed."
                << " count:" << count << nl
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


template<class Type>
void Foam::PstreamDetail::gatherv
(
    const Type* sendData,
    int sendCount,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm))
    {
        return;
    }
    if (!UPstream::is_parallel(comm))
    {
        // recvCounts[0] may be invalid - use sendCount instead
        std::memmove(recvData, sendData, sendCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Igatherv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Gatherv (blocking):";
        }
        Perr<< " np:" << np
            << " recvCounts:" << recvCounts
            << " recvOffsets:" << recvOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        UPstream::master(comm)
     && (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }

    // Ensure send/recv consistency on master
    if (UPstream::master(comm) && !recvCounts[0])
    {
        sendCount = 0;
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Igatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Igatherv failed [comm: " << comm << ']'
                << " sendCount " << sendCount
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Gatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed [comm: " << comm << ']'
                << " sendCount " << sendCount
                << " recvCounts " << recvCounts
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


template<class Type>
void Foam::PstreamDetail::scatterv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_rank(comm))
    {
        return;
    }
    if (!UPstream::is_parallel(comm))
    {
        std::memmove(recvData, sendData, recvCount*sizeof(Type));
        return;
    }

    const label np = UPstream::nProcs(comm);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Iscatterv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Scatterv (blocking):";
        }
        Perr<< " np:" << np
            << " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        UPstream::master(comm)
     && (sendCounts.size() != np || sendOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << " or sendOffsets:" << sendOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Iscatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iscatterv [comm: " << comm << "] failed."
                << " sendCounts " << sendCounts
                << " sendOffsets " << sendOffsets
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Scatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv [comm: " << comm << "] failed."
                << " sendCounts " << sendCounts
                << " sendOffsets " << sendOffsets
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


template<class Type>
void Foam::PstreamDetail::allGather
(
    Type* allData,
    int count,

    MPI_Datatype datatype,
    const label comm,

    UPstream::Request* req,
    label* requestID
)
{
    PstreamGlobals::reset_request(req, requestID);

    const bool immediate = (req || requestID);

    if (!UPstream::is_parallel(comm))
    {
        // Nothing to do - ignore
        return;
    }

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        if (immediate)
        {
            Perr<< "** MPI_Iallgather (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Allgather (blocking):";
        }
        Perr<< " numProc:" << UPstream::nProcs(comm)
            << " with comm:" << comm
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }


#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        if
        (
            MPI_Iallgather
            (
                MPI_IN_PLACE, count, datatype,
                allData, count, datatype,
                PstreamGlobals::MPICommunicators_[comm],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iallgather [comm: " << comm << "] failed."
                << Foam::abort(FatalError);
        }

        PstreamGlobals::push_request(request, req, requestID);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Allgather
            (
                MPI_IN_PLACE, count, datatype,
                allData, count, datatype,
                PstreamGlobals::MPICommunicators_[comm]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Allgather [comm: " << comm << "] failed."
                << Foam::abort(FatalError);
        }

        // Is actually gather/scatter but we can't split it apart
        profilingPstream::addGatherTime();
    }
}


// ************************************************************************* //
