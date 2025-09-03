/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "debug.H"
#include "registerSwitch.H"
#include "dictionary.H"
#include "SHA1.H"
#include "OSspecific.H"  // for hostName()
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);
}

const Foam::Enum
<
    Foam::UPstream::commsTypes
>
Foam::UPstream::commsTypeNames
({
    { commsTypes::buffered, "buffered" },   // "buffered"
    { commsTypes::scheduled, "scheduled" },
    { commsTypes::nonBlocking, "nonBlocking" },  // "immediate"
    // compatibility names
    { commsTypes::buffered, "blocking" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Determine host grouping.
// Uses SHA1 of hostname instead of MPI_Comm_split or MPI_Comm_split_type
// for two reasons:
// - Comm_split returns an MPI_COMM_NULL on non-participating process
//   which does not easily fit into the OpenFOAM framework
//
// - use the SHA1 of hostname allows a single MPI_Gather, determination of
//   the inter-host vs intra-host (on the master) followed by a single
//   broadcast of integers.
//
// Returns: the unique host indices with the leading hosts encoded
// with negative values
static List<int> getHostGroupIds(const label parentCommunicator)
{
    const label numProcs = UPstream::nProcs(parentCommunicator);

    List<SHA1Digest> digests;
    if (UPstream::master(parentCommunicator))
    {
        digests.resize(numProcs);
    }

    // Could also add lowercase etc, but since hostName()
    // will be consistent within the same node, there is no need.
    SHA1Digest myDigest(SHA1(hostName()).digest());

    // The fixed-length digest allows use of MPI_Gather
    UPstream::mpiGather
    (
        myDigest.cdata_bytes(),     // Send
        digests.data_bytes(),       // Recv
        SHA1Digest::max_size(),     // Num send/recv data per rank
        parentCommunicator
    );

    List<int> hostIDs(numProcs);

    // Compact numbering of hosts.
    if (UPstream::master(parentCommunicator))
    {
        DynamicList<SHA1Digest> uniqDigests;

        forAll(digests, proci)
        {
            const SHA1Digest& dig = digests[proci];

            hostIDs[proci] = uniqDigests.find(dig);

            if (hostIDs[proci] < 0)
            {
                // First appearance of host. Encode as leader
                hostIDs[proci] = -(uniqDigests.size() + 1);
                uniqDigests.push_back(dig);
            }
        }
    }

    UPstream::broadcast
    (
        hostIDs.data_bytes(),
        hostIDs.size_bytes(),
        parentCommunicator,
        UPstream::masterNo()
    );

    return hostIDs;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const label nProcs, const bool haveThreads)
{
    parRun_ = (nProcs > 0);
    haveThreads_ = haveThreads;

    label comm = -1;
    labelRange singleProc(1);

    // Redo communicators that were created during static initialisation.
    // When parRun == true, redo with MPI components
    // When parRun == false, just redo in case of future changes

    if (!parRun_)
    {
        // These are already correct from the static initialisation,
        // but just in case of future changes

        // Using (world, self) ordering
        freeCommunicator(UPstream::commSelf());
        freeCommunicator(UPstream::commGlobal());

        // 0: COMM_WORLD : commWorld() / commGlobal()
        comm = allocateCommunicator(-1, singleProc, false);
        if (comm != UPstream::commGlobal())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-global:" << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }

        // 1: COMM_SELF
        comm = allocateCommunicator(-2, singleProc, false);
        if (comm != UPstream::commSelf())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-self:" << UPstream::commSelf()
                << Foam::exit(FatalError);
        }

        Pout.prefix().clear();
        Perr.prefix().clear();
    }
    else
    {
        // Redo communicators that were created during static initialisation
        // but this time with MPI components

        // Using (world, self) ordering
        freeCommunicator(UPstream::commSelf());
        freeCommunicator(UPstream::commGlobal());

        // 0: COMM_WORLD : commWorld() / commGlobal()
        comm = allocateCommunicator(-1, labelRange(nProcs), true);
        if (comm != UPstream::commGlobal())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-global:" << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }

        // 1: COMM_SELF
        // - Processor number wrt world communicator
        singleProc.start() = UPstream::myProcNo(UPstream::commGlobal());
        comm = allocateCommunicator(-2, singleProc, true);
        if (comm != UPstream::commSelf())
        {
            // Failed sanity check
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  expected comm-self:" << UPstream::commSelf()
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' + Foam::name(myProcNo(commGlobal())) + "] ";
        Perr.prefix() = Pout.prefix();
    }

    if (debug)
    {
        Perr<< "UPstream::setParRun :"
            << " nProcs:" << nProcs
            << " haveThreads:" << haveThreads
            << endl;
    }
}


Foam::label Foam::UPstream::getAvailableCommIndex(const label parentIndex)
{
    label index;
    if (!freeComms_.empty())
    {
        // LIFO pop
        index = freeComms_.back();
        freeComms_.pop_back();

        // Reset existing
        myProcNo_[index] = -1;
        parentComm_[index] = parentIndex;

        procIDs_[index].clear();
        linearCommunication_[index].clear();
        treeCommunication_[index].clear();
    }
    else
    {
        // Extend storage
        index = parentComm_.size();

        myProcNo_.push_back(-1);
        parentComm_.push_back(parentIndex);

        procIDs_.emplace_back();
        linearCommunication_.emplace_back();
        treeCommunication_.emplace_back();
    }

    return index;
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelRange& subRanks,
    const bool withComponents
)
{
    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Allocating communicator " << index << nl
            << "    parent : " << parentIndex << nl
            << "    procs  : " << subRanks << nl
            << endl;
    }

    // Initially treat as master,
    // overwritten by allocateCommunicatorComponents
    myProcNo_[index] = UPstream::masterNo();

    // The selected sub-ranks.
    // - transcribe from label to int
    // - already in incremental order
    auto& procIds = procIDs_[index];
    procIds.resize_nocopy(subRanks.size());

    label numSubRanks = 0;
    for (const label subRanki : subRanks)
    {
        procIds[numSubRanks] = subRanki;
        ++numSubRanks;
    }

    // Sizing and filling are demand-driven
    linearCommunication_[index].clear();
    treeCommunication_[index].clear();

    if (withComponents && parRun())
    {
        allocateCommunicatorComponents(parentIndex, index);

        // Could 'remember' locations of uninvolved ranks
        /// if (myProcNo_[index] < 0 && parentIndex >= 0)
        /// {
        ///     // As global rank
        ///     myProcNo_[index] = -(myProcNo_[worldComm]+1);
        ///
        /// OR:
        ///     // As parent rank number
        ///     if (myProcNo_[parentIndex] >= 0)
        ///     {
        ///         myProcNo_[index] = -(myProcNo_[parentIndex]+1);
        ///     }
        /// }
    }

    return index;
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelUList& subRanks,
    const bool withComponents
)
{
    const label index = getAvailableCommIndex(parentIndex);

    if (debug)
    {
        Perr<< "Allocating communicator " << index << nl
            << "    parent : " << parentIndex << nl
            << "    procs  : " << flatOutput(subRanks) << nl
            << endl;
    }

    // Initially treat as master,
    // overwritten by allocateCommunicatorComponents
    myProcNo_[index] = UPstream::masterNo();

    // The selected sub-ranks.
    // - transcribe from label to int. Treat negative values as 'ignore'
    // - enforce incremental order (so index is rank in next communicator)

    auto& procIds = procIDs_[index];
    procIds.resize_nocopy(subRanks.size());

    label numSubRanks = 0;
    bool monotonicOrder = true;
    for (const label subRanki : subRanks)
    {
        if (subRanki < 0)
        {
            continue;
        }
        if (monotonicOrder && numSubRanks)
        {
            monotonicOrder = (procIds[numSubRanks-1] < subRanki);
        }

        procIds[numSubRanks] = subRanki;
        ++numSubRanks;
    }

    if (!monotonicOrder)
    {
        auto last = procIds.begin() + numSubRanks;
        std::sort(procIds.begin(), last);
        last = std::unique(procIds.begin(), last);
        numSubRanks = label(last - procIds.begin());
    }

    procIds.resize(numSubRanks);

    // Sizing and filling are demand-driven
    linearCommunication_[index].clear();
    treeCommunication_[index].clear();

    if (withComponents && parRun())
    {
        allocateCommunicatorComponents(parentIndex, index);

        // Could 'remember' locations of uninvolved ranks
        /// if (myProcNo_[index] < 0 && parentIndex >= 0)
        /// {
        ///     // As global rank
        ///     myProcNo_[index] = -(myProcNo_[worldComm]+1);
        ///
        /// OR:
        ///     // As parent rank number
        ///     if (myProcNo_[parentIndex] >= 0)
        ///     {
        ///         myProcNo_[index] = -(myProcNo_[parentIndex]+1);
        ///     }
        /// }
    }

    return index;
}


Foam::label Foam::UPstream::allocateInterHostCommunicator
(
    const label parentCommunicator
)
{
    List<int> hostIDs = getHostGroupIds(parentCommunicator);

    DynamicList<label> subRanks(hostIDs.size());

    // From master to host-leader. Ranks between hosts.
    forAll(hostIDs, proci)
    {
        // Is host leader?
        if (hostIDs[proci] < 0)
        {
            subRanks.push_back(proci);
        }
    }

    return allocateCommunicator(parentCommunicator, subRanks);
}


Foam::label Foam::UPstream::allocateIntraHostCommunicator
(
    const label parentCommunicator
)
{
    List<int> hostIDs = getHostGroupIds(parentCommunicator);

    DynamicList<label> subRanks(hostIDs.size());

    // Intra-host ranks. Ranks within a host
    int myHostId = hostIDs[UPstream::myProcNo(parentCommunicator)];
    if (myHostId < 0) myHostId = -(myHostId + 1);  // Flip to generic id

    forAll(hostIDs, proci)
    {
        int id = hostIDs[proci];
        if (id < 0) id = -(id + 1);  // Flip to generic id

        if (id == myHostId)
        {
            subRanks.push_back(proci);
        }
    }

    return allocateCommunicator(parentCommunicator, subRanks);
}


bool Foam::UPstream::allocateHostCommunicatorPairs()
{
    // Use the world communicator (not global communicator)
    const label parentCommunicator = worldComm;

    // Skip if non-parallel
    if (!parRun())
    {
        return false;
    }

    if (interHostComm_ >= 0 || intraHostComm_ >= 0)
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Host communicator(s) already created!" << endl
            << Foam::exit(FatalError);
        return false;
    }

    interHostComm_ = getAvailableCommIndex(parentCommunicator);
    intraHostComm_ = getAvailableCommIndex(parentCommunicator);

    // Sorted order, purely cosmetic
    if (intraHostComm_ < interHostComm_)
    {
        std::swap(intraHostComm_, interHostComm_);
    }

    // Overwritten later
    myProcNo_[intraHostComm_] = UPstream::masterNo();
    myProcNo_[interHostComm_] = UPstream::masterNo();

    if (debug)
    {
        Perr<< "Allocating host communicators "
            << interHostComm_ << ", " << intraHostComm_ << nl
            << "    parent : " << parentCommunicator << nl
            << endl;
    }

    List<int> hostIDs = getHostGroupIds(parentCommunicator);

    DynamicList<int> subRanks(hostIDs.size());

    // From master to host-leader. Ranks between hosts.
    {
        subRanks.clear();
        forAll(hostIDs, proci)
        {
            // Is host leader?
            if (hostIDs[proci] < 0)
            {
                subRanks.push_back(proci);

                // Flip to generic host id
                hostIDs[proci] = -(hostIDs[proci] + 1);
            }
        }

        const label index = interHostComm_;

        // Direct copy (subRanks is also int)
        procIDs_[index] = subRanks;

        // Implicitly: withComponents = true
        if (parRun())  // Already checked...
        {
            allocateCommunicatorComponents(parentCommunicator, index);
        }

        // Sizing and filling are demand-driven
        linearCommunication_[index].clear();
        treeCommunication_[index].clear();
    }

    // Intra-host ranks. Ranks within a host
    {
        int myHostId = hostIDs[UPstream::myProcNo(parentCommunicator)];
        if (myHostId < 0) myHostId = -(myHostId + 1);  // Flip to generic id

        subRanks.clear();
        forAll(hostIDs, proci)
        {
            int id = hostIDs[proci];
            if (id < 0) id = -(id + 1);  // Flip to generic id

            if (id == myHostId)
            {
                subRanks.push_back(proci);
            }
        }

        const label index = intraHostComm_;

        // Direct copy (subRanks is also int)
        procIDs_[index] = subRanks;

        // Implicitly: withComponents = true
        if (parRun())  // Already checked...
        {
            allocateCommunicatorComponents(parentCommunicator, index);
        }

        // Sizing and filling are demand-driven
        linearCommunication_[index].clear();
        treeCommunication_[index].clear();
    }

    return true;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool withComponents
)
{
    // Filter out any placeholders
    if (communicator < 0)
    {
        return;
    }

    // Update demand-driven communicators
    if (interHostComm_ == communicator) interHostComm_ = -1;
    if (intraHostComm_ == communicator) intraHostComm_ = -1;

    if (debug)
    {
        Perr<< "Communicators : Freeing communicator " << communicator
            << " parent: " << parentComm_[communicator]
            << " myProcNo: " << myProcNo_[communicator]
            << endl;
    }

    if (withComponents && parRun())
    {
        freeCommunicatorComponents(communicator);
    }

    myProcNo_[communicator] = -1;
    parentComm_[communicator] = -1;
    //procIDs_[communicator].clear();
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    // LIFO push
    freeComms_.push_back(communicator);
}


int Foam::UPstream::baseProcNo(label comm, int procID)
{
    while (UPstream::parent(comm) >= 0 && procID >= 0)
    {
        const auto& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label comm, const int baseProcID)
{
    const auto& parentRanks = UPstream::procID(comm);
    label parentComm = UPstream::parent(comm);

    int procID = baseProcID;

    if (parentComm >= 0)
    {
        procID = UPstream::procNo(parentComm, baseProcID);
    }

    return parentRanks.find(procID);
}


Foam::label Foam::UPstream::procNo
(
    const label comm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return UPstream::procNo(comm, physProcID);
}


const Foam::List<Foam::UPstream::commsStruct>&
Foam::UPstream::linearCommunication(const label communicator)
{
    if (linearCommunication_[communicator].empty())
    {
        linearCommunication_[communicator] =
            List<commsStruct>(UPstream::nProcs(communicator));
    }

    return linearCommunication_[communicator];
}


const Foam::List<Foam::UPstream::commsStruct>&
Foam::UPstream::treeCommunication(const label communicator)
{
    if (treeCommunication_[communicator].empty())
    {
        treeCommunication_[communicator] =
            List<commsStruct>(UPstream::nProcs(communicator));
    }

    return treeCommunication_[communicator];
}


void Foam::UPstream::printCommTree(const label communicator)
{
    const auto& comms = UPstream::whichCommunication(communicator);

    if (UPstream::master(communicator))
    {
        commsStruct::printGraph(Info(), comms);
    }
}


Foam::label Foam::UPstream::commIntraHost()
{
    if (!parRun())
    {
        return worldComm;  // Don't know anything better to return
    }
    if (intraHostComm_ < 0)
    {
        allocateHostCommunicatorPairs();
    }
    return intraHostComm_;
}


Foam::label Foam::UPstream::commInterHost()
{
    if (!parRun())
    {
        return worldComm;  // Don't know anything better to return
    }
    if (interHostComm_ < 0)
    {
        allocateHostCommunicatorPairs();
    }
    return interHostComm_;
}


bool Foam::UPstream::hasHostComms()
{
    return (intraHostComm_ >= 0 || interHostComm_ >= 0);
}


void Foam::UPstream::clearHostComms()
{
    // Always with Pstream
    freeCommunicator(intraHostComm_, true);
    freeCommunicator(interHostComm_, true);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

int Foam::UPstream::msgType_(1);


Foam::wordList Foam::UPstream::allWorlds_(Foam::one{}, "");
Foam::labelList Foam::UPstream::worldIDs_(Foam::one{}, 0);


Foam::DynamicList<int> Foam::UPstream::myProcNo_(16);
Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(16);

Foam::DynamicList<Foam::label> Foam::UPstream::parentComm_(16);
Foam::DynamicList<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(16);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(16);


Foam::label Foam::UPstream::intraHostComm_(-1);
Foam::label Foam::UPstream::interHostComm_(-1);

Foam::label Foam::UPstream::worldComm(0);
Foam::label Foam::UPstream::warnComm(-1);


// Predefine world and self communicator slots.
// These are overwritten in parallel mode (by UPstream::setParRun())
const Foam::label nPredefinedComm = []()
{
    // 0: COMM_WORLD : commWorld() / commGlobal()
    (void) Foam::UPstream::allocateCommunicator(-1, Foam::labelRange(1), false);

    // 1: COMM_SELF
    (void) Foam::UPstream::allocateCommunicator(-2, Foam::labelRange(1), false);

    return Foam::UPstream::nComms();
}();


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 0)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

int Foam::UPstream::nProcsNonblockingExchange
(
    Foam::debug::optimisationSwitch("nbx.min", 0)
);
registerOptSwitch
(
    "nbx.min",
    int,
    Foam::UPstream::nProcsNonblockingExchange
);

int Foam::UPstream::tuning_NBX_
(
    Foam::debug::optimisationSwitch("nbx.tuning", 0)
);
registerOptSwitch
(
    "nbx.tuning",
    int,
    Foam::UPstream::tuning_NBX_
);


int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.get
    (
        "commsType",
        Foam::debug::optimisationSwitches()
    )
);


//! \cond file_scope
namespace Foam
{
    //- Registered reader for UPstream::defaultCommsType
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType =
                UPstream::commsTypeNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}
//! \endcond

int Foam::UPstream::maxCommsSize
(
    Foam::debug::optimisationSwitch("maxCommsSize", 0)
);
registerOptSwitch
(
    "maxCommsSize",
    int,
    Foam::UPstream::maxCommsSize
);


const int Foam::UPstream::mpiBufferSize
(
    Foam::debug::optimisationSwitch("mpiBufferSize", 0)
);


// ************************************************************************* //
