/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "int.H"
#include "UPstreamWrapping.H"
#include "collatedFileOperation.H"

#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>

#undef Pstream_use_MPI_Get_count

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The min value and default for MPI buffer length
constexpr int minBufLen = 20000000;

// Track size of attached MPI buffer
static int attachedBufLen = 0;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Attach user-defined send buffer
static void attachOurBuffers()
{
#ifndef SGIMPI
    if (attachedBufLen)
    {
        return;  // Already attached
    }

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str, len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    char* buf = new char[len];

    if (MPI_SUCCESS == MPI_Buffer_attach(buf, len))
    {
        // Properly attached
        attachedBufLen = len;

        if (Foam::UPstream::debug)
        {
            Foam::Perr<< "UPstream::init : buffer-size " << len << '\n';
        }
    }
    else
    {
        delete[] buf;
        Foam::Perr<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


// Remove an existing user-defined send buffer
// IMPORTANT:
//     This operation will block until all messages currently in the
//     buffer have been transmitted.
static void detachOurBuffers()
{
#ifndef SGIMPI
    if (!attachedBufLen)
    {
        return;  // Nothing to detach
    }

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

    char* buf = nullptr;
    int len = 0;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        // This was presumably the buffer that we attached
        // and not someone else.
        delete[] buf;
    }

    // Nothing attached
    attachedBufLen = 0;
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        if (UPstream::debug)
        {
            Perr<< "UPstream::initNull : was already initialized\n";
        }
    }
    else
    {
        // Not already initialized

        MPI_Init_thread
        (
            nullptr,    // argc
            nullptr,    // argv
            MPI_THREAD_SINGLE,
            &flag       // provided_thread_support
        );

        ourMpi = true;
    }

    // Could also attach buffers etc.

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int numprocs = 0, myRank = 0;
    int provided_thread_support = 0;
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized.
        // Warn if we've called twice, but skip if initialized externally

        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already initialized - cannot perform MPI_Init" << nl
                << "This could indicate an application programming error!"
                << endl;

            return true;
        }
        else if (UPstream::debug)
        {
            Perr<< "UPstream::init : was already initialized\n";
        }
    }
    else
    {
        MPI_Init_thread
        (
            &argc,
            &argv,
            (
                needsThread
              ? MPI_THREAD_MULTIPLE
              : MPI_THREAD_SINGLE
            ),
            &provided_thread_support
        );

        ourMpi = true;
    }

    // Check argument list for local world
    label worldIndex = -1;
    word world;
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-world") == 0)
        {
            worldIndex = argi++;
            if (argi >= argc)
            {
                FatalErrorInFunction
                    << "Missing world name to argument \"world\""
                    << Foam::abort(FatalError);
            }
            world = argv[argi];
            break;
        }
    }

    // Filter 'world' option
    if (worldIndex != -1)
    {
        for (label i = worldIndex+2; i < argc; i++)
        {
            argv[i-2] = argv[i];
        }
        argc -= 2;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (UPstream::debug)
    {
        Perr<< "UPstream::init :"
            << " thread-support : requested:" << needsThread
            << " obtained:"
            << (
                   (provided_thread_support == MPI_THREAD_SINGLE)
                 ? "SINGLE"
                 : (provided_thread_support == MPI_THREAD_SERIALIZED)
                 ? "SERIALIZED"
                 : (provided_thread_support == MPI_THREAD_MULTIPLE)
                 ? "MULTIPLE"
                 : "other"
               )
            << " procs:" << numprocs
            << " rank:" << myRank
            << " world:" << world << endl;
    }

    if (worldIndex == -1 && numprocs <= 1)
    {
        FatalErrorInFunction
            << "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    if (worldIndex != -1)
    {
        // During startup, so commWorld() == commGlobal()

        wordList worlds(numprocs);
        worlds[UPstream::myProcNo(UPstream::commGlobal())] = world;
        Pstream::gatherList
        (
            worlds,
            UPstream::msgType(),
            UPstream::commGlobal()
        );

        // Compact
        if (UPstream::master(UPstream::commGlobal()))
        {
            DynamicList<word> worldNames(numprocs);
            worldIDs_.resize_nocopy(numprocs);

            forAll(worlds, proci)
            {
                const word& world = worlds[proci];

                worldIDs_[proci] = worldNames.find(world);

                if (worldIDs_[proci] == -1)
                {
                    worldIDs_[proci] = worldNames.size();
                    worldNames.push_back(world);
                }
            }

            allWorlds_.transfer(worldNames);
        }
        Pstream::broadcasts(UPstream::commGlobal(), allWorlds_, worldIDs_);

        const label myWorldId =
            worldIDs_[UPstream::myProcNo(UPstream::commGlobal())];

        DynamicList<label> subRanks;
        forAll(worldIDs_, proci)
        {
            if (worldIDs_[proci] == myWorldId)
            {
                subRanks.push_back(proci);
            }
        }

        // Allocate new communicator with comm-global as its parent
        const label subComm =
            UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);


        // Override worldComm
        UPstream::worldComm = subComm;
        // For testing: warn use of non-worldComm
        UPstream::warnComm = UPstream::worldComm;

        // MPI_COMM_SELF : the processor number wrt the new world communicator
        if (procIDs_[UPstream::commSelf()].size())
        {
            procIDs_[UPstream::commSelf()].front() =
                UPstream::myProcNo(subComm);
        }

        if (UPstream::debug)
        {
            // Check
            int subNumProcs, subRank;
            MPI_Comm_size
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subNumProcs
            );
            MPI_Comm_rank
            (
                PstreamGlobals::MPICommunicators_[subComm],
                &subRank
            );

            Perr<< "UPstream::init : in world:" << world
                << " using local communicator:" << subComm
                << " rank " << subRank
                << " of " << subNumProcs
                << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + world + '/' +  name(myProcNo(subComm)) + "] ";
        Perr.prefix() = Pout.prefix();
    }
    else
    {
        // All processors use world 0
        worldIDs_.resize_nocopy(numprocs);
        worldIDs_ = 0;
    }

    attachOurBuffers();

    return true;
}


void Foam::UPstream::shutdown(int errNo)
{
    int flag = 0;

    MPI_Initialized(&flag);
    if (!flag)
    {
        // MPI not initialized - we have nothing to do
        return;
    }

    MPI_Finalized(&flag);
    if (flag)
    {
        // MPI already finalized - we have nothing to do
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (UPstream::debug && errNo == 0)
        {
            Perr<< "UPstream::shutdown : was already finalized\n";
        }
        ourMpi = false;
        return;
    }

    if (!ourMpi)
    {
        WarningInFunction
            << "Finalizing MPI, but was initialized elsewhere\n";
    }
    ourMpi = false;


    // Abort - stop now, without any final synchonization steps!
    // -----

    if (errNo != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, errNo);
        return;
    }


    // Regular cleanup
    // ---------------

    if (UPstream::debug)
    {
        Perr<< "UPstream::shutdown\n";
    }

    // Check for any outstanding requests
    {
        label nOutstanding = 0;

        for (MPI_Request request : PstreamGlobals::outstandingRequests_)
        {
            if (MPI_REQUEST_NULL != request)
            {
                // TBD: MPI_Cancel(&request); MPI_Request_free(&request);
                ++nOutstanding;
            }
        }

        if (nOutstanding)
        {
            WarningInFunction
                << "Still have " << nOutstanding
                << " outstanding MPI requests."
                << " Should not happen for a normal code exit."
                << endl;
        }

        PstreamGlobals::outstandingRequests_.clear();
    }


    {
        detachOurBuffers();

        forAllReverse(myProcNo_, communicator)
        {
            freeCommunicatorComponents(communicator);
        }
    }


    MPI_Finalize();
}


void Foam::UPstream::exit(int errNo)
{
    UPstream::shutdown(errNo);
    std::exit(errNo);
}


void Foam::UPstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::allocateCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPICommunicators_.size())
    {
        // Extend storage with null values
        PstreamGlobals::pendingMPIFree_.emplace_back(false);
        PstreamGlobals::MPICommunicators_.emplace_back(MPI_COMM_NULL);
    }
    else if (index > PstreamGlobals::MPICommunicators_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Global communicator. Same as world communicator for single-world

        if (index != UPstream::commGlobal())
        {
            FatalErrorInFunction
                << "base world communicator should always be index "
                << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }

        PstreamGlobals::pendingMPIFree_[index] = false;
        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;

        // TBD: MPI_Comm_dup(MPI_COMM_WORLD, ...);
        // with pendingMPIFree_[index] = true
        // Note: freeCommunicatorComponents() may need an update

        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of ranks to the actual number
        int numProcs;
        MPI_Comm_size
        (
            PstreamGlobals::MPICommunicators_[index],
           &numProcs
        );

        // identity [0-numProcs], as 'int'
        procIDs_[index].resize_nocopy(numProcs);
        std::iota(procIDs_[index].begin(), procIDs_[index].end(), 0);
    }
    else if (parentIndex == -2)
    {
        // MPI_COMM_SELF

        PstreamGlobals::pendingMPIFree_[index] = false;
        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_SELF;

        MPI_Comm_rank(MPI_COMM_SELF, &myProcNo_[index]);

        // Number of ranks is always 1 (self communicator)

        #ifdef FULLDEBUG
        int numProcs;
        MPI_Comm_size(MPI_COMM_SELF, &numProcs);

        if (numProcs != 1)
        {
            // Already finalized - this is an error
            FatalErrorInFunction
                << "MPI_COMM_SELF had " << numProcs << " != 1 ranks!\n"
                << Foam::abort(FatalError);
        }
        #endif

        // For MPI_COMM_SELF : the process IDs within the world communicator.
        // Uses MPI_COMM_WORLD in case called before UPstream::commGlobal()
        // was initialized

        procIDs_[index].resize_nocopy(1);
        MPI_Comm_rank(MPI_COMM_WORLD, &procIDs_[index].front());
    }
    else
    {
        // General sub-communicator

        PstreamGlobals::pendingMPIFree_[index] = true;

        // Starting from parent
        MPI_Group parent_group;
        MPI_Comm_group
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
           &parent_group
        );

        MPI_Group active_group;
        MPI_Group_incl
        (
            parent_group,
            procIDs_[index].size(),
            procIDs_[index].cdata(),
           &active_group
        );

        #if defined(MSMPI_VER)
        // ms-mpi (10.0 and others?) does not have MPI_Comm_create_group
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            active_group,
           &PstreamGlobals::MPICommunicators_[index]
        );
        #else
        // Create new communicator for this group
        MPI_Comm_create_group
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            active_group,
            UPstream::msgType(),
           &PstreamGlobals::MPICommunicators_[index]
        );
        #endif

        // Groups not needed after this...
        MPI_Group_free(&parent_group);
        MPI_Group_free(&active_group);

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            // No communicator created
            myProcNo_[index] = -1;
            PstreamGlobals::pendingMPIFree_[index] = false;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freeCommunicatorComponents(const label index)
{
    if (UPstream::debug)
    {
        Perr<< "freeCommunicatorComponents: " << index
            << " from " << PstreamGlobals::MPICommunicators_.size() << endl;
    }

    // Only free communicators that we have specifically allocated ourselves
    //
    // Bounds checking needed since there are no UPstream communicator indices
    // when MPI is initialized outside of OpenFOAM

    if
    (
        (index >= 0 && index < PstreamGlobals::MPICommunicators_.size())
     && PstreamGlobals::pendingMPIFree_[index]
    )
    {
        PstreamGlobals::pendingMPIFree_[index] = false;

        // Free communicator. Sets communicator to MPI_COMM_NULL
        if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[index])
        {
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[index]);
        }
    }
}


void Foam::UPstream::barrier(const label communicator, UPstream::Request* req)
{
    // No-op for non-parallel or not on communicator
    if (!UPstream::parRun() || !UPstream::is_rank(communicator))
    {
        PstreamGlobals::reset_request(req);
        return;
    }

    if (req)
    {
        MPI_Request request;

        // Non-blocking
        if
        (
            MPI_Ibarrier
            (
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ibarrier returned with error"
                << Foam::abort(FatalError);
        }

        *req = UPstream::Request(request);
    }
    else
    {
        // Blocking
        if
        (
            MPI_Barrier
            (
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Barrier returned with error"
                << Foam::abort(FatalError);
        }
    }
}


std::pair<int,int64_t>
Foam::UPstream::probeMessage
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    const int tag,
    const label communicator
)
{
    std::pair<int,int64_t> result(-1, 0);

    // No-op for non-parallel or not on communicator
    if (!UPstream::parRun() || !UPstream::is_rank(communicator))
    {
        return result;
    }

    const int source = (fromProcNo < 0) ? MPI_ANY_SOURCE : fromProcNo;
    // Supporting MPI_ANY_TAG is not particularly useful...

    int flag = 0;
    MPI_Status status;

    if (UPstream::commsTypes::buffered == commsType)
    {
        // Blocking
        profilingPstream::beginTiming();

        if
        (
            MPI_Probe
            (
                source,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Probe returned with error"
                << Foam::abort(FatalError);
        }

        profilingPstream::addProbeTime();
        flag = 1;
    }
    else
    {
        // Non-blocking
        profilingPstream::beginTiming();

        if
        (
            MPI_Iprobe
            (
                source,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &flag,
                &status
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iprobe returned with error"
                << Foam::abort(FatalError);
        }

        profilingPstream::addRequestTime();
    }

    if (flag)
    {
        // Unlikely to be used with large amounts of data,
        // but use MPI_Get_elements_x() instead of MPI_Count() anyhow

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
        else if (int64_t(count) > int64_t(INT_MAX))
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "count is larger than INI_MAX bytes"
                << Foam::abort(FatalError);
        }


        result.first = status.MPI_SOURCE;
        result.second = int64_t(count);
    }

    return result;
}


// ************************************************************************* //
