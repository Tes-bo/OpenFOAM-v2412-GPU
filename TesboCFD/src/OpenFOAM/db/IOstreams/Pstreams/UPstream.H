/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

Class
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstreamCommsStruct.C
    UPstreamTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstream_H
#define Foam_UPstream_H

#include "wordList.H"
#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "Map.H"
#include "Enum.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Implementation details for UPstream/Pstream/MPI etc.
namespace PstreamDetail {}

//- Interface handling for UPstream/Pstream/MPI etc.
namespace PstreamUtils {}

/*---------------------------------------------------------------------------*\
                          Class UPstream Declaration
\*---------------------------------------------------------------------------*/

class UPstream
{
public:

    //- Int ranges are used for MPI ranks (processes)
    typedef IntRange<int> rangeType;

    //- Communications types
    enum class commsTypes : char
    {
        buffered,       //!< "buffered"                 : (MPI_Bsend, MPI_Recv)
        scheduled,      //!< "scheduled" (MPI standard) : (MPI_Send, MPI_Recv)
        nonBlocking,    //!< "nonBlocking" (immediate)  : (MPI_Isend, MPI_Irecv)
        // Aliases
        blocking = buffered  //!< compatibility name for buffered
    };

    //- Enumerated names for the communication types
    static const Enum<commsTypes> commsTypeNames;

    //- Different MPI-send modes (ignored for commsTypes::buffered)
    enum class sendModes : char
    {
        normal,         //!< (MPI_Send, MPI_Isend)
        sync            //!< (MPI_Ssend, MPI_Issend)
    };


    // Public Classes

        //- Wrapper for MPI_Comm
        class Communicator;  // Forward Declaration

        //- Wrapper for MPI_Request
        class Request;  // Forward Declaration

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private Data

                //- The procID of the processor \em directly above
                label above_;

                //- The procIDs of processors \em directly below
                labelList below_;

                //- The procIDs of all processors below myProcNo,
                //- not just directly below
                labelList allBelow_;

                //- The procIDs of all processors not below myProcNo
                //  (inverse of allBelow_ without myProcNo)
                labelList allNotBelow_;


        public:

            // Constructors

                //- Default construct with above == -1
                commsStruct() noexcept : above_(-1) {}

                //- Move construct from components
                commsStruct
                (
                    const label above,
                    labelList&& below,
                    labelList&& allBelow,
                    labelList&& allNotBelow
                );

                //- Copy construct from below, allBelow components
                commsStruct
                (
                    const label numProcs,
                    const label myProcID,
                    const label above,
                    const labelUList& below,
                    const labelUList& allBelow
                );


            // Member Functions

                //- Print un-directed graph in graphviz dot format
                static void printGraph
                (
                    Ostream& os,
                    const UList<UPstream::commsStruct>& comms,
                    const label proci = 0  // starting node
                );


            // Access

                //- The number of processors addressed by the structure
                label nProcs() const;

                //- The procID of the processor \em directly above
                label above() const noexcept { return above_; }

                //- The procIDs of the processors \em directly below
                const labelList& below() const noexcept { return below_; }

                //- The procIDs of all processors below
                //- (so not just directly below)
                const labelList& allBelow() const noexcept { return allBelow_; }

                //- The procIDs of all processors not below myProcNo.
                //- The inverse set of allBelow without myProcNo.
                const labelList& allNotBelow() const noexcept
                {
                    return allNotBelow_;
                }


            // Edit

                //- Reset to default constructed state
                void reset();

                //- Reset with automatic linear/tree selection
                void reset(const label procID, const label numProcs);


            // Member / Friend Operators

                bool operator==(const commsStruct&) const;
                bool operator!=(const commsStruct&) const;

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };


private:

    // Private Data

        //- Communications type of this stream
        commsTypes commsType_;


    // Private Static Data

        //- By default this is not a parallel run
        static bool parRun_;

        //- Have support for threads?
        static bool haveThreads_;

        //- Standard transfer message type
        static int msgType_;

        //- Names of all worlds
        static wordList allWorlds_;

        //- Per processor the world index (into allWorlds_)
        static labelList worldIDs_;

        //- Intra-host communicator
        static label intraHostComm_;

        //- Inter-host communicator (between host leaders)
        static label interHostComm_;


    // Communicator specific data

        //- My processor number
        static DynamicList<int> myProcNo_;

        //- List of process IDs
        static DynamicList<List<int>> procIDs_;

        //- Parent communicator
        static DynamicList<label> parentComm_;

        //- Free communicators
        static DynamicList<label> freeComms_;

        //- Linear communication schedule
        static DynamicList<List<commsStruct>> linearCommunication_;

        //- Multi level communication schedule
        static DynamicList<List<commsStruct>> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun(const label nProcs, const bool haveThreads);

        //- Initialise entries for new communicator. Return the index
        static label getAvailableCommIndex(const label parentIndex);

        //- Allocate MPI components of communicator with given index
        static void allocateCommunicatorComponents
        (
            const label parentIndex,
            const label index
        );

        //- Free MPI components of communicator.
        //  Does not touch the first two communicators (SELF, WORLD)
        static void freeCommunicatorComponents(const label index);

        //- Allocate inter-host, intra-host communicators
        //- with comm-world as parent
        static bool allocateHostCommunicatorPairs();


public:

    //- Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static Data

        //- Should compact transfer be used in which floats replace doubles
        //- reducing the bandwidth requirement at the expense of some loss
        //- in accuracy
        static bool floatTransfer;

        //- Number of processors to change from linear to tree communication
        static int nProcsSimpleSum;

        //- Number of processors to change to nonBlocking consensual
        //- exchange (NBX). Ignored for zero or negative values.
        static int nProcsNonblockingExchange;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Optional maximum message size (bytes)
        static int maxCommsSize;

        //- Tuning parameters for non-blocking exchange (NBX)
        static int tuning_NBX_;

        //- MPI buffer-size (bytes)
        static const int mpiBufferSize;


    // Standard Communicators

        //- Communicator for all ranks.
        //- May differ from commGlobal() if local worlds are in use
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;

        //- Communicator for all ranks, irrespective of any local worlds
        static constexpr label commGlobal() noexcept { return 0; }

        //- Communicator within the current rank only
        static constexpr label commSelf() noexcept { return 1; }

        //- Communicator for all ranks (respecting any local worlds)
        static label commWorld() noexcept { return worldComm; }

        //- Set world communicator. Negative values are a no-op.
        //  \returns old world communicator index
        static label commWorld(const label communicator) noexcept
        {
            label old(worldComm);
            if (communicator >= 0) worldComm = communicator;
            return old;
        }

        //- Alter communicator debugging setting.
        //- Warns for use of any communicator differing from specified.
        //  \returns the previous warn index
        static label commWarn(const label communicator) noexcept
        {
            label old(warnComm);
            warnComm = communicator;
            return old;
        }

        //- Number of currently defined communicators
        static label nComms() noexcept { return parentComm_.size(); }

        //- Debugging: print the communication tree
        static void printCommTree(const label communicator);


    // Host Communicators

        //- Demand-driven: Intra-host communicator (respects any local worlds)
        static label commIntraHost();

        //- Demand-driven: Inter-host communicator (respects any local worlds)
        static label commInterHost();

        //- Test for presence of any intra or inter host communicators
        static bool hasHostComms();

        //- Remove any existing intra and inter host communicators
        static void clearHostComms();


    // Constructors

        //- Construct for given communication type
        explicit UPstream(const commsTypes commsType) noexcept
        :
            commsType_(commsType)
        {}


    // Member Functions

        //- Allocate new communicator with contiguous sub-ranks
        //- on the parent communicator.
        static label allocateCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The contiguous sub-ranks of parent to use
            const labelRange& subRanks,

            //! Call allocateCommunicatorComponents()
            const bool withComponents = true
        );

        //- Allocate new communicator with sub-ranks on the parent communicator
        static label allocateCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The sub-ranks of parent to use (ignore negative values)
            const labelUList& subRanks,

            //! Call allocateCommunicatorComponents()
            const bool withComponents = true
        );

        //- Free a previously allocated communicator.
        //  Ignores placeholder (negative) communicators.
        static void freeCommunicator
        (
            const label communicator,
            const bool withComponents = true
        );

        //- Allocate an inter-host communicator
        static label allocateInterHostCommunicator
        (
            const label parentCommunicator = worldComm
        );

        //- Allocate an intra-host communicator
        static label allocateIntraHostCommunicator
        (
            const label parentCommunicator = worldComm
        );


        //- Wrapper class for allocating/freeing communicators. Always invokes
        //- allocateCommunicatorComponents() and freeCommunicatorComponents()
        class communicator
        {
            label comm_;

        public:

            //- No copy construct
            communicator(const communicator&) = delete;

            //- No copy assignment
            void operator=(const communicator&) = delete;

            //- Default construct (a placeholder communicator)
            communicator() noexcept : comm_(-1) {}

            //- Move construct, takes ownership
            communicator(communicator&& c) : comm_(c.comm_) { c.comm_ = -1; }

            //- Allocate communicator for contiguous sub-ranks on given parent
            communicator
            (
                //! The parent communicator
                const label parentComm,

                //! The contiguous sub-ranks of parent to use
                const labelRange& subRanks
            )
            :
                comm_(UPstream::allocateCommunicator(parentComm, subRanks))
            {}

            //- Allocate communicator for sub-ranks on given parent
            communicator
            (
                //! The parent communicator
                const label parentComm,

                //! The sub-ranks of parent to use (negative values ignored)
                const labelUList& subRanks
            )
            :
                comm_(UPstream::allocateCommunicator(parentComm, subRanks))
            {}

            //- Free allocated communicator
            ~communicator() { UPstream::freeCommunicator(comm_); }

            //- True if communicator is non-negative (ie, was allocated)
            bool good() const noexcept { return (comm_ >= 0); }

            //- The communicator label
            label comm() const noexcept { return comm_; }

            //- Release ownership of the communicator, return old value.
            //  Leave further management to the caller
            label release() noexcept { label c(comm_); comm_ = -1; return c; }

            //- Free allocated communicator
            void reset() { UPstream::freeCommunicator(comm_); comm_ = -1; }

            //- Allocate with contiguous sub-ranks of parent communicator
            void reset(label parent, const labelRange& subRanks)
            {
                UPstream::freeCommunicator(comm_);
                comm_ = UPstream::allocateCommunicator(parent, subRanks);
            }

            //- Allocate with sub-ranks of parent communicator
            void reset(label parent, const labelUList& subRanks)
            {
                UPstream::freeCommunicator(comm_);
                comm_ = UPstream::allocateCommunicator(parent, subRanks);
            }

            //- Take ownership, free allocated communicator
            //  \caution do not call as self-assignment
            void reset(communicator&& c)
            {
                if (comm_ != c.comm_) UPstream::freeCommunicator(comm_);
                comm_ = c.comm_;
                c.comm_ = -1;
            }

            //- Move assignment, takes ownership
            //  \caution do not call as self-assignment
            void operator=(communicator&& c) { reset(std::move(c)); }

            //- Cast to label - the same as comm()
            operator label() const noexcept { return comm_; }
        };


        //- Return physical processor number (i.e. processor number in
        //- worldComm) given communicator and processor
        static int baseProcNo(label comm, int procID);

        //- Return processor number in communicator (given physical processor
        //- number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //- and communicator)
        static label procNo
        (
            const label comm,
            const label currentComm,
            const int currentProcID
        );

        //- Add the valid option this type of communications library
        //- adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns sub-processes and initialises inter-communication
        static bool init(int& argc, char**& argv, const bool needsThread);

        //- Special purpose initialisation function.
        //  Performs a basic MPI_Init without any other setup.
        //  Only used for applications that need MPI communication when
        //  OpenFOAM is running in a non-parallel mode.
        //  \note Behaves as a no-op if MPI has already been initialized.
        //      Fatal if MPI has already been finalized.
        static bool initNull();

        //- Impose a synchronisation barrier (optionally non-blocking)
        static void barrier
        (
            const label communicator,
            UPstream::Request* req = nullptr
        );

        //- Probe for an incoming message.
        //
        //  \param commsType    Non-blocking or not
        //  \param fromProcNo   The source rank (negative == ANY_SOURCE)
        //  \param tag          The source message tag
        //  \param communicator The communicator index
        //
        //  \returns source rank and message size (bytes)
        //      and (-1, 0) on failure
        static std::pair<int,int64_t> probeMessage
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            const int tag = UPstream::msgType(),
            const label communicator = worldComm
        );


    // Requests (non-blocking comms).
    // Pending requests are usually handled as an internal (global) list,
    // since this simplifies the overall tracking and provides a convenient
    // wrapping to avoid exposing MPI-specific types, but can also handle
    // with a wrapped UPstream::Request as well.

        //- Number of outstanding requests (on the internal list of requests)
        static label nRequests() noexcept;

        //- Truncate outstanding requests to given length, which is
        //- expected to be in the range [0 to nRequests()].
        //  A no-op for out-of-range values.
        static void resetRequests(const label n);

        //- Transfer the (wrapped) MPI request to the internal global list.
        //  A no-op for non-parallel. No special treatment for null requests.
        static void addRequest(UPstream::Request& req);

        //- Non-blocking comms: cancel and free outstanding request.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static void cancelRequest(const label i);

        //- Non-blocking comms: cancel and free outstanding request.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false
        static void cancelRequest(UPstream::Request& req);

        //- Non-blocking comms: cancel and free outstanding requests.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false or list is empty
        static void cancelRequests(UList<UPstream::Request>& requests);

        //- Non-blocking comms: cancel/free outstanding requests
        //- (from position onwards) and remove from internal list of requests.
        //- Corresponds to MPI_Cancel() + MPI_Request_free()
        //  A no-op if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to remove (negative = until the end)
        static void removeRequests(const label pos, label len = -1);

        //- Non-blocking comms: free outstanding request.
        //- Corresponds to MPI_Request_free()
        //  A no-op if parRun() == false
        static void freeRequest(UPstream::Request& req);

        //- Non-blocking comms: free outstanding requests.
        //- Corresponds to MPI_Request_free()
        //  A no-op if parRun() == false or list is empty
        static void freeRequests(UList<UPstream::Request>& requests);

        //- Wait until all requests (from position onwards) have finished.
        //- Corresponds to MPI_Waitall()
        //  A no-op if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  If checking a trailing portion of the list, it will also trim
        //  the list of outstanding requests as a side-effect.
        //  This is a feature (not a bug) to conveniently manange the list.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        static void waitRequests(const label pos, label len = -1);

        //- Wait until all requests have finished.
        //- Corresponds to MPI_Waitall()
        //  A no-op if parRun() == false, or the list is empty.
        static void waitRequests(UList<UPstream::Request>& requests);

        //- Wait until any request (from position onwards) has finished.
        //- Corresponds to MPI_Waitany()
        //  A no-op and returns false if parRun() == false,
        //  if the position is out-of-range [0 to nRequests()],
        //  or the internal list of requests is empty.
        //
        //  \returns true if any pending request completed.
        //  \returns false if all requests have already been handled.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        static bool waitAnyRequest(const label pos, label len = -1);

        //- Wait until some requests (from position onwards) have finished.
        //- Corresponds to MPI_Waitsome()
        //  A no-op and returns false if parRun() == false,
        //  if the position is out-of-range [0 to nRequests],
        //  or the internal list of requests is empty.
        //
        //  \returns true if some pending requests completed.
        //  \returns false if all requests have already been handled
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        //  \param[out] indices the completed request indices relative to the
        //      starting position. This is an optional parameter that can be
        //      used to recover the indices or simply to avoid reallocations
        //      when calling within a loop.
        static bool waitSomeRequests
        (
            const label pos,
            label len = -1,
            DynamicList<int>* indices = nullptr
        );

        //- Wait until some requests have finished.
        //- Corresponds to MPI_Waitsome()
        //  A no-op and returns false if parRun() == false,
        //  the list is empty,
        //  or if all the requests have already been handled.
        //
        //  \param requests the requests
        //  \param[out] indices the completed request indices relative to the
        //      starting position. This is an optional parameter that can be
        //      used to recover the indices or simply to avoid reallocations
        //      when calling within a loop.
        static bool waitSomeRequests
        (
            UList<UPstream::Request>& requests,
            DynamicList<int>* indices = nullptr
        );

        //- Wait until any request has finished and return its index.
        //- Corresponds to MPI_Waitany()
        //  Returns -1 if parRun() == false, or the list is empty,
        //  or if all the requests have already been handled
        static label waitAnyRequest(UList<UPstream::Request>& requests);

        //- Wait until request i has finished.
        //- Corresponds to MPI_Wait()
        //  A no-op if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static void waitRequest(const label i);

        //- Wait until specified request has finished.
        //- Corresponds to MPI_Wait()
        //  A no-op if parRun() == false or for a null-request
        static void waitRequest(UPstream::Request& req);

        //- Non-blocking comms: has request i finished?
        //- Corresponds to MPI_Test()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        static bool finishedRequest(const label i);

        //- Non-blocking comms: has request finished?
        //- Corresponds to MPI_Test()
        //  A no-op and returns true if parRun() == false
        //  or for a null-request
        static bool finishedRequest(UPstream::Request& req);

        //- Non-blocking comms: have all requests (from position onwards)
        //- finished?
        //- Corresponds to MPI_Testall()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the index is out-of-range (0 to nRequests)
        //  or the addressed range is empty etc.
        //
        //  \param pos starting position within the internal list of requests
        //  \param len length of slice to check (negative = until the end)
        static bool finishedRequests(const label pos, label len = -1);

        //- Non-blocking comms: have all requests finished?
        //- Corresponds to MPI_Testall()
        //  A no-op and returns true if parRun() == false or list is empty
        static bool finishedRequests(UList<UPstream::Request>& requests);

        //- Non-blocking comms: have both requests finished?
        //- Corresponds to pair of MPI_Test()
        //  A no-op and returns true if parRun() == false,
        //  if there are no pending requests,
        //  or if the indices are out-of-range (0 to nRequests)
        //  Each finished request parameter is set to -1 (ie, done).
        static bool finishedRequestPair(label& req0, label& req1);

        //- Non-blocking comms: wait for both requests to finish.
        //- Corresponds to pair of MPI_Wait()
        //  A no-op if parRun() == false,
        //  if there are no pending requests,
        //  or if the indices are out-of-range (0 to nRequests)
        //  Each finished request parameter is set to -1 (ie, done).
        static void waitRequestPair(label& req0, label& req1);


    // General

        //- Set as parallel run on/off.
        //  \return the previous value
        static bool parRun(const bool on) noexcept
        {
            bool old(parRun_);
            parRun_ = on;
            return old;
        }

        //- Test if this a parallel run
        //  Modify access is deprecated
        static bool& parRun() noexcept { return parRun_; }

        //- Have support for threads
        static bool haveThreads() noexcept { return haveThreads_; }

        //- Relative rank for the master process - is always 0.
        static constexpr int masterNo() noexcept { return 0; }

        //- Number of ranks in parallel run (for given communicator).
        //- It is 1 for serial run
        static label nProcs(const label communicator = worldComm)
        {
            return procIDs_[communicator].size();
        }

        //- Rank of this process in the communicator (starting from masterNo()).
        //- Can be negative if the process is not a rank in the communicator
        static int myProcNo(const label communicator = worldComm)
        {
            return myProcNo_[communicator];
        }

        //- True if process corresponds to the master rank in the communicator
        static bool master(const label communicator = worldComm)
        {
            return myProcNo_[communicator] == masterNo();
        }

        //- True if process corresponds to any rank (master or sub-rank)
        //- in the given communicator
        static bool is_rank(const label communicator = worldComm)
        {
            return myProcNo_[communicator] >= 0;
        }

        //- True if process corresponds to a sub-rank in the given communicator
        static bool is_subrank(const label communicator = worldComm)
        {
            return myProcNo_[communicator] > masterNo();
        }

        //- True if parallel algorithm or exchange is required.
        //  This is when parRun() == true, the process corresponds to a rank
        //  in the communicator and there is more than one rank in the
        //  communicator
        static bool is_parallel(const label communicator = worldComm)
        {
            return
            (
                parRun_ && is_rank(communicator) && nProcs(communicator) > 1
            );
        }

        //- The parent communicator
        static label parent(const label communicator)
        {
            return parentComm_(communicator);
        }

        //- The list of ranks within a given communicator
        static List<int>& procID(const label communicator)
        {
            return procIDs_[communicator];
        }


    // Worlds

        //- All worlds
        static const wordList& allWorlds() noexcept
        {
            return allWorlds_;
        }

        //- The indices into allWorlds for all processes
        static const labelList& worldIDs() noexcept
        {
            return worldIDs_;
        }

        //- My worldID
        static label myWorldID()
        {
            return worldIDs_[myProcNo(UPstream::commGlobal())];
        }

        //- My world
        static const word& myWorld()
        {
            return allWorlds_[worldIDs_[myProcNo(UPstream::commGlobal())]];
        }


    // Rank addressing

        //- Range of process indices for all processes
        static rangeType allProcs(const label communicator = worldComm)
        {
            // Proc 0 -> nProcs (int value)
            return rangeType(static_cast<int>(nProcs(communicator)));
        }

        //- Range of process indices for sub-processes
        static rangeType subProcs(const label communicator = worldComm)
        {
            // Proc 1 -> nProcs (int value)
            return rangeType(1, static_cast<int>(nProcs(communicator)-1));
        }

        //- Communication schedule for linear all-to-master (proc 0)
        static const List<commsStruct>&
        linearCommunication
        (
            const label communicator = worldComm
        );

        //- Communication schedule for tree all-to-master (proc 0)
        static const List<commsStruct>&
        treeCommunication
        (
            const label communicator = worldComm
        );

        //- Communication schedule for all-to-master (proc 0) as
        //- linear/tree/none with switching based on UPstream::nProcsSimpleSum
        //- and the is_parallel() state
        static const List<commsStruct>& whichCommunication
        (
            const label communicator = worldComm
        )
        {
            const label np
            (
                parRun_ && is_rank(communicator)  // cf. is_parallel()
              ? nProcs(communicator)
              : 0
            );

            return
            (
                np <= 1
              ? List<commsStruct>::null()
              : np < nProcsSimpleSum
              ? linearCommunication(communicator)
              : treeCommunication(communicator)
            );
        }


        //- Message tag of standard messages
        static int& msgType() noexcept
        {
            return msgType_;
        }

        //- Set the message tag for standard messages
        //  \return the previous value
        static int msgType(int val) noexcept
        {
            int old(msgType_);
            msgType_ = val;
            return old;
        }

        //- Increment the message tag for standard messages
        //  \return the previous value
        static int incrMsgType(int val = 1) noexcept
        {
            int old(msgType_);
            msgType_ += val;
            return old;
        }

        //- Get the communications type of the stream
        commsTypes commsType() const noexcept
        {
            return commsType_;
        }

        //- Set the communications type of the stream
        //  \return the previous value
        commsTypes commsType(const commsTypes ct) noexcept
        {
            commsTypes old(commsType_);
            commsType_ = ct;
            return old;
        }


        //- Shutdown (finalize) MPI as required.
        //  Uses MPI_Abort instead of MPI_Finalize if errNo is non-zero
        static void shutdown(int errNo = 0);

        //- Call MPI_Abort with no other checks or cleanup
        static void abort();

        //- Shutdown (finalize) MPI as required and exit program with errNo.
        static void exit(int errNo = 1);


        #undef  Pstream_CommonRoutines
        #define Pstream_CommonRoutines(Native)                                \
                                                                              \
        /*!\brief Exchange \c Native data with all ranks in communicator   */ \
        /*! \em non-parallel : simple copy of \p sendData to \p recvData   */ \
        static void allToAll                                                  \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci */                 \
            const UList<Native>& sendData,                                    \
            /*! [out] The data received from the other ranks */               \
            UList<Native>& recvData,                                          \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \em non-zero \c Native data between ranks [NBX] */ \
        /*! \p recvData is always initially assigned zero and no non-zero  */ \
        /*! values are sent/received from other ranks.                     */ \
        /*! \em non-parallel : simple copy of \p sendData to \p recvData   */ \
        /*! \note The message tag should be chosen to be a unique value    */ \
        /*! since the implementation uses probing with ANY_SOURCE !!       */ \
        /*! An initial barrier may help to avoid synchronisation problems  */ \
        /*! caused elsewhere (See "nbx.tuning" opt switch)                 */ \
        static void allToAllConsensus                                         \
        (                                                                     \
            /*! [in] The \em non-zero value at [proci] is sent to proci */    \
            const UList<Native>& sendData,                                    \
            /*! [out] The non-zero value received from each rank */           \
            UList<Native>& recvData,                                          \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \c Native data between ranks [NBX]              */ \
        /*! \p recvData map is always cleared initially so a simple check  */ \
        /*! of its keys is sufficient to determine connectivity.           */ \
        /*! \em non-parallel : copy own rank (if it exists)                */ \
        /*! See notes about message tags and "nbx.tuning" opt switch       */ \
        static void allToAllConsensus                                         \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci. */                \
            const Map<Native>& sendData,                                      \
            /*! [out] The values received from given ranks. */                \
            Map<Native>& recvData,                                            \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*!\brief Exchange \c Native data between ranks [NBX]              */ \
        /*! \returns any received data as a Map                            */ \
        static Map<Native> allToAllConsensus                                  \
        (                                                                     \
            /*! [in] The value at [proci] is sent to proci. */                \
            const Map<Native>& sendData,                                      \
            /*! Message tag for the communication */                          \
            const int tag,                                                    \
            const label communicator = worldComm                              \
        )                                                                     \
        {                                                                     \
            Map<Native> recvData;                                             \
            allToAllConsensus(sendData, recvData, tag, communicator);         \
            return recvData;                                                  \
        }

        Pstream_CommonRoutines(int32_t);
        Pstream_CommonRoutines(int64_t);

        #undef Pstream_CommonRoutines


    // Low-level gather/scatter routines

        #undef  Pstream_CommonRoutines
        #define Pstream_CommonRoutines(Native)                                \
                                                                              \
        /*! \brief Receive identically-sized \c Native data from all ranks */ \
        static void mpiGather                                                 \
        (                                                                     \
            /*! On rank: individual value to send */                          \
            const Native* sendData,                                           \
            /*! On master: receive buffer with all values */                  \
            Native* recvData,                                                 \
            /*! Number of send/recv data per rank. Globally consistent! */    \
            int count,                                                        \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Send identically-sized \c Native data to all ranks */      \
        static void mpiScatter                                                \
        (                                                                     \
            /*! On master: send buffer with all values */                     \
            const Native* sendData,                                           \
            /*! On rank: individual value to receive */                       \
            Native* recvData,                                                 \
            /*! Number of send/recv data per rank. Globally consistent! */    \
            int count,                                                        \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Gather/scatter identically-sized \c Native data */         \
        /*! Send data from proc slot, receive into all slots */               \
        static void mpiAllGather                                              \
        (                                                                     \
            /*! On all ranks: the base of the data locations */               \
            Native* allData,                                                  \
            /*! Number of send/recv data per rank. Globally consistent! */    \
            int count,                                                        \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Receive variable length \c Native data from all ranks */   \
        static void gather                                                    \
        (                                                                     \
            const Native* sendData,                                           \
            int sendCount,  /*!< Ignored on master if recvCount[0] == 0 */    \
            Native* recvData,               /*!< Ignored on non-root rank */  \
            const UList<int>& recvCounts,   /*!< Ignored on non-root rank */  \
            const UList<int>& recvOffsets,  /*!< Ignored on non-root rank */  \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Send variable length \c Native data to all ranks */        \
        static void scatter                                                   \
        (                                                                     \
            const Native* sendData,         /*!< Ignored on non-root rank */  \
            const UList<int>& sendCounts,   /*!< Ignored on non-root rank */  \
            const UList<int>& sendOffsets,  /*!< Ignored on non-root rank */  \
            Native* recvData,                                                 \
            int recvCount,                                                    \
            const label communicator = worldComm                              \
        );

        Pstream_CommonRoutines(char);
        Pstream_CommonRoutines(int32_t);
        Pstream_CommonRoutines(int64_t);
        Pstream_CommonRoutines(uint32_t);
        Pstream_CommonRoutines(uint64_t);
        Pstream_CommonRoutines(float);
        Pstream_CommonRoutines(double);

        #undef Pstream_CommonRoutines


    // Gather single, contiguous value(s)

        //- Allgather individual values into list locations.
        //  The returned list has size nProcs, identical on all ranks.
        template<class T>
        static List<T> allGatherValues
        (
            const T& localValue,
            const label communicator = worldComm
        );

        //- Gather individual values into list locations.
        //  On master list length == nProcs, otherwise zero length.
        //  \n
        //  For \b non-parallel :
        //  the returned list length is 1 with localValue.
        template<class T>
        static List<T> listGatherValues
        (
            const T& localValue,
            const label communicator = worldComm
        );

        //- Scatter individual values from list locations.
        //  On master input list length == nProcs, ignored on other procs.
        //  \n
        //  For \b non-parallel :
        //  returns the first list element (or default initialized).
        template<class T>
        static T listScatterValues
        (
            const UList<T>& allValues,
            const label communicator = worldComm
        );


    // Broadcast Functions

        //- Broadcast buffer contents to all processes in given communicator.
        //- The sizes must match on all processes.
        //  For \b non-parallel : do nothing.
        //  \return True on success
        static bool broadcast
        (
            char* buf,
            const std::streamsize bufSize,
            const label communicator,
            const int rootProcNo = masterNo()
        );


    // Logical reductions

        //- Logical (and) reduction (MPI_AllReduce)
        //  For \b non-parallel : do nothing
        static void reduceAnd
        (
            bool& value,
            const label communicator = worldComm
        );

        //- Logical (or) reduction (MPI_AllReduce)
        //  For \b non-parallel : do nothing
        static void reduceOr
        (
            bool& value,
            const label communicator = worldComm
        );


    // Housekeeping

        //- Wait for all requests to finish.
        //  \deprecated(2023-01) Probably not what you want.
        //     Should normally be restricted to a particular starting request.
        FOAM_DEPRECATED_FOR(2023-01, "waitRequests(int) method")
        static void waitRequests() { waitRequests(0); }

        //- Process index of first sub-process
        //  \deprecated(2020-09) use subProcs() method instead
        FOAM_DEPRECATED_FOR(2020-09, "subProcs() method")
        static constexpr int firstSlave() noexcept
        {
            return 1;
        }

        //- Process index of last sub-process
        //  \deprecated(2020-09) use subProcs() method instead
        FOAM_DEPRECATED_FOR(2020-09, "subProcs() or allProcs() method")
        static int lastSlave(const label communicator = worldComm)
        {
            return nProcs(communicator) - 1;
        }
};


/*---------------------------------------------------------------------------*\
                   Class UPstream::Communicator Declaration
\*---------------------------------------------------------------------------*/

//- An opaque wrapper for MPI_Comm with a vendor-independent
//- representation without any \c <mpi.h> header.
//  The MPI standard states that MPI_Comm is always an opaque object.
//  Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).
class UPstream::Communicator
{
public:

    // Public Types

        //- Storage for MPI_Comm (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Comm (as wrapped value)
        value_type value_;

public:

    // Generated Methods

        //- Copy construct
        Communicator(const Communicator&) noexcept = default;

        //- Move construct
        Communicator(Communicator&&) noexcept = default;

        //- Copy assignment
        Communicator& operator=(const Communicator&) noexcept = default;

        //- Move assignment
        Communicator& operator=(Communicator&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Communicator& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Communicator& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_COMM_NULL
        Communicator() noexcept;

        //- Construct from MPI_Comm (as pointer type)
        explicit Communicator(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Comm (as integer type)
        explicit Communicator(value_type val) noexcept
        :
            value_(val)
        {}


    // Factory Methods

        //- Transcribe internally indexed communicator to wrapped value.
        //  Example,
        //  \code
        //  PstreamUtils::Cast::to_mpi
        //  (
        //      UPstream::Communicator::lookup(UPstream::commWorld())
        //  )
        //  \endcode
        static Communicator lookup(const label comm);


    // Member Functions

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }

        //- True if not equal to MPI_COMM_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_COMM_NULL)
        void reset() noexcept;
};


/*---------------------------------------------------------------------------*\
                      Class UPstream::Request Declaration
\*---------------------------------------------------------------------------*/

//- An opaque wrapper for MPI_Request with a vendor-independent
//- representation without any \c <mpi.h> header.
//  The MPI standard states that MPI_Request is always an opaque object.
//  Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).
class UPstream::Request
{
public:

    // Public Types

        //- Storage for MPI_Request (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Request (as wrapped value)
        value_type value_;

public:

    // Generated Methods

        //- Copy construct
        Request(const Request&) noexcept = default;

        //- Move construct
        Request(Request&&) noexcept = default;

        //- Copy assignment
        Request& operator=(const Request&) noexcept = default;

        //- Move assignment
        Request& operator=(Request&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Request& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Request& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_REQUEST_NULL
        Request() noexcept;

        //- Construct from MPI_Request (as pointer type)
        explicit Request(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Request (as integer type)
        explicit Request(value_type val) noexcept
        :
            value_(val)
        {}


    // Member Functions

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }

        //- True if not equal to MPI_REQUEST_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_REQUEST_NULL)
        void reset() noexcept;

        //- Same as calling UPstream::cancelRequest()
        void cancel() { UPstream::cancelRequest(*this); }

        //- Same as calling UPstream::freeRequest()
        void free() { UPstream::freeRequest(*this); }

        //- Same as calling UPstream::finishedRequest()
        bool finished() { return UPstream::finishedRequest(*this); }

        //- Same as calling UPstream::waitRequest()
        void wait() { UPstream::waitRequest(*this); }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream&, const UPstream::commsStruct&);


// * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * //

// Template specialisation for access of commsStruct
template<>
UPstream::commsStruct&
UList<UPstream::commsStruct>::operator[](const label procID);

template<>
const UPstream::commsStruct&
UList<UPstream::commsStruct>::operator[](const label procID) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstreamTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
