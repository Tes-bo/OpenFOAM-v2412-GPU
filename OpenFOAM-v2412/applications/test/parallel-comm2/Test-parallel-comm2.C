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

Application
    Test-parallel-comm2

Description
    Basic communicator tests

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pair.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "SHA1.H"

#include "openfoam_mpi.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rankInfo(const label comm)
{
    const int ranki = UPstream::myProcNo(comm);

    Pout<< "comm:" << comm
        << "(parent:" << UPstream::parent(comm) << ')'
        << " rank:" << ranki
        << "(sub:" << UPstream::is_subrank(comm)
        << ") nProcs:" << UPstream::nProcs(comm);
      // << " baseProcNo:" << UPstream::baseProcNo(comm, ranki);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Set UPstream::debug level");
    argList::addBoolOption("info", "information");
    argList::addBoolOption("print-tree", "Report tree(s) as graph");
    argList::addBoolOption("comm-split", "Test simple comm split");
    argList::addBoolOption("mpi-host-comm", "Test DIY host-comm split");
    argList::addBoolOption("host-comm", "Test Pstream host-comm");
    argList::addBoolOption("host-broadcast", "Test host-base broadcasts");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    #include "setRootCase.H"

    const bool optPrintTree = args.found("print-tree");

    Info<< nl
        << "parallel:" << UPstream::parRun()
        << " nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)."
        << " proc:" << UPstream::myProcNo() << nl;

    if (UPstream::parRun() && optPrintTree)
    {
        Info<< "comms: " << UPstream::whichCommunication() << endl;
        UPstream::printCommTree(UPstream::commWorld());
    }

    if (args.found("info"))
    {
        Info<< nl;

        // Process IDs within a given communicator
        for (label comm = 0; comm < UPstream::nComms(); ++comm)
        {
            Info<< "comm=" << comm
                << " procIDs: " << flatOutput(UPstream::procID(comm)) << endl;
            rankInfo(comm);
            Pout<< nl;
        }
        Pout<< endl;
    }

    bool generalTest = true;

    if (UPstream::parRun() && args.found("comm-split"))
    {
        generalTest = false;

        int world_nprocs = 0;
        int world_rank = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &world_nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        int host_nprocs = 0;
        int host_rank = -1;
        MPI_Comm commIntraHost;
        MPI_Comm_split_type
        (
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,  // OMPI_COMM_TYPE_NODE
            0, MPI_INFO_NULL, &commIntraHost
        );

        MPI_Comm_size(commIntraHost, &host_nprocs);
        MPI_Comm_rank(commIntraHost, &host_rank);

        int leader_nprocs = 0;
        int leader_rank = -1;
        MPI_Comm commInterHost;

        if (false)
        {
            // Easy enough to use MPI_Comm_split, but slightly annoying
            // that it returns MPI_COMM_NULL for unused ranks...
            MPI_Comm commInterHost;
            MPI_Comm_split
            (
                MPI_COMM_WORLD,
                (host_rank == 0) ? 0 : MPI_UNDEFINED,
                0, &commInterHost
            );

            if (commInterHost != MPI_COMM_NULL)
            {
                MPI_Comm_size(commInterHost, &leader_nprocs);
                MPI_Comm_rank(commInterHost, &leader_rank);
            }
        }
        else
        {
            boolList isHostLeader(world_nprocs, false);
            isHostLeader[world_rank] = (host_rank == 0);

            MPI_Allgather
            (
                // recv is also send
                MPI_IN_PLACE, 1, MPI_C_BOOL,
                isHostLeader.data(), 1, MPI_C_BOOL,
                MPI_COMM_WORLD
            );

            Pout<< "leaders: " << isHostLeader << endl;

            DynamicList<int> subRanks(isHostLeader.size());
            forAll(isHostLeader, proci)
            {
                if (isHostLeader[proci])
                {
                    subRanks.push_back(proci);
                }
            }
            // Starting from parent
            MPI_Group parent_group;
            MPI_Comm_group(MPI_COMM_WORLD, &parent_group);

            MPI_Group active_group;
            MPI_Group_incl
            (
                parent_group,
                subRanks.size(),
                subRanks.cdata(),
                &active_group
            );

            // Create new communicator for this group
            MPI_Comm_create_group
            (
                MPI_COMM_WORLD,
                active_group,
                UPstream::msgType(),
                &commInterHost
            );

            // Groups not needed after this...
            MPI_Group_free(&parent_group);
            MPI_Group_free(&active_group);

            MPI_Comm_size(commInterHost, &leader_nprocs);
            MPI_Comm_rank(commInterHost, &leader_rank);
        }

        Pout<< nl << "[MPI_Comm_split_type]" << nl
            << "Host rank " << host_rank << " / " << host_nprocs
            << " on " << hostName()
            << " inter-rank: " << leader_rank << " / " << leader_nprocs
            << " host leader:" << (leader_rank == 0)
            << " sub-rank:" << (leader_rank > 0)
            << nl;

        if (commInterHost != MPI_COMM_NULL)
        {
            MPI_Comm_free(&commInterHost);
        }
        if (commIntraHost != MPI_COMM_NULL)
        {
            MPI_Comm_free(&commIntraHost);
        }
    }

    if (UPstream::parRun() && args.found("mpi-host-comm"))
    {
        generalTest = false;

        // Host communicator, based on the current world communicator
        // Use hostname
        // Lowest rank per hostname is the IO rank

        label numprocs = UPstream::nProcs(UPstream::commGlobal());

        // Option 1: using hostnames
        // - pro: trivial coding
        // - con: unequal lengths, more allocations and 'hops'
        stringList hosts(numprocs);
        hosts[Pstream::myProcNo(UPstream::commGlobal())] = hostName();
        Pstream::gatherList(hosts, UPstream::msgType(), UPstream::commGlobal());


        // Option 2: using SHA1 of hostnames
        // - con: uglier coding (but only needed locally!)
        // - pro: fixed digest length enables direct MPI calls
        //        can avoid Pstream::gatherList() during setup...

        List<SHA1Digest> digests;
        if (UPstream::master(UPstream::commGlobal()))
        {
            digests.resize(numprocs);
        }

        {
            const SHA1Digest myDigest(SHA1(hostName()).digest());

            UPstream::mpiGather
            (
                myDigest.cdata_bytes(),     // Send
                digests.data_bytes(),       // Recv
                SHA1Digest::max_size(),     // Num send/recv per rank
                UPstream::commGlobal()
            );
        }


        labelList hostIDs(numprocs);
        DynamicList<label> subRanks(numprocs);

        Info<< "digests: " << digests << nl;

        // Compact numbering
        if (UPstream::master(UPstream::commGlobal()))
        {
            DynamicList<word> hostNames(numprocs);

            forAll(hosts, proci)
            {
                const word& host = hosts[proci];

                hostIDs[proci] = hostNames.find(host);

                if (hostIDs[proci] < 0)
                {
                    // First appearance of host (encode as leader)
                    hostIDs[proci] = -(hostNames.size() + 1);
                    hostNames.push_back(host);
                 }
            }
            hostIDs = -1;


            DynamicList<SHA1Digest> uniqDigests(numprocs);

            forAll(digests, proci)
            {
                const SHA1Digest& dig = digests[proci];

                hostIDs[proci] = uniqDigests.find(dig);

                if (hostIDs[proci] < 0)
                {
                    // First appearance of host (encode as leader)
                    hostIDs[proci] = -(uniqDigests.size() + 1);
                    uniqDigests.push_back(dig);
                }
            }
        }


        Info<< "hosts =  " << hosts << endl;
        Info<< "hostIDs =  " << hostIDs << endl;

        UPstream::broadcast
        (
            hostIDs.data_bytes(),
            hostIDs.size_bytes(),
            UPstream::commGlobal(),
            UPstream::masterNo()
        );

        // Ranks for world to inter-host communicator
        // - very straightforward

        #if 0
        subRanks.clear();
        forAll(hostIDs, proci)
        {
            // Is host leader?
            if (hostIDs[proci] < 0)
            {
                subRanks.push_back(proci);

                // Flip back to generic host id
                hostIDs[proci] = -(hostIDs[proci] + 1);
            }
        }

        // From world to hostMaster
        const label commInterHost =
            UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);
        #endif

        const label myWorldProci = UPstream::myProcNo(UPstream::commGlobal());

        label myHostId = hostIDs[myWorldProci];
        if (myHostId < 0) myHostId = -(myHostId + 1);  // Flip to generic id

        // Ranks for within a host
        subRanks.clear();
        forAll(hostIDs, proci)
        {
            label id = hostIDs[proci];
            if (id < 0) id = -(id + 1);  // Flip to generic id

            if (id == myHostId)
            {
                subRanks.push_back(proci);
            }
        }

        // The intra-host ranks
        const label commIntraHost =
            UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);


        // Test what if we have intra-host comm and we want host-master

        List<bool> isHostMaster(numprocs, false);
        if (UPstream::master(commIntraHost))
        {
            isHostMaster[myWorldProci] = true;
        }

        UPstream::mpiAllGather
        (
            isHostMaster.data_bytes(),
            sizeof(bool),
            UPstream::commGlobal()
        );

        // Ranks for world to hostMaster
        // - very straightforward
        subRanks.clear();
        forAll(isHostMaster, proci)
        {
            if (isHostMaster[proci])
            {
                subRanks.push_back(proci);
            }
        }

        // From world to hostMaster
        const label commInterHost =
            UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);


        Pout<< nl << "[manual split]" << nl
            << nl << "Host rank " << UPstream::myProcNo(commIntraHost)
            << " / " << UPstream::nProcs(commIntraHost)
            << " on " << hostName()
            << ", inter-rank: " << UPstream::myProcNo(commInterHost)
            << " / " << UPstream::nProcs(commInterHost)
            << " host leader:" << UPstream::master(commInterHost)
            << " sub-rank:" << UPstream::is_subrank(commInterHost)
            << nl;

        UPstream::freeCommunicator(commInterHost);
        UPstream::freeCommunicator(commIntraHost);
    }

    if (UPstream::parRun() && args.found("host-comm"))
    {
        generalTest = false;
        Info<< nl << "[pstream host-comm]" << nl << endl;

        const label commInterHost = UPstream::commInterHost();
        const label commIntraHost = UPstream::commIntraHost();

        Pout<< "Host rank " << UPstream::myProcNo(commIntraHost)
            << " / " << UPstream::nProcs(commIntraHost)
            << " on " << hostName()
            << ", inter-rank: " << UPstream::myProcNo(commInterHost)
            << " / " << UPstream::nProcs(commInterHost)
            << ", host leader:" << UPstream::master(commInterHost)
            << " sub-rank:" << UPstream::is_subrank(commInterHost)
            << endl;


        {
            Info<< "host-master: "
                << UPstream::whichCommunication(commInterHost) << endl;

            UPstream::printCommTree(commInterHost);
            UPstream::printCommTree(commIntraHost);
        }
    }

    if (UPstream::parRun() && args.found("host-broadcast"))
    {
        generalTest = false;
        Info<< nl << "[pstream host-broadcast]" << nl << endl;

        const label commInterHost = UPstream::commInterHost();
        const label commIntraHost = UPstream::commIntraHost();

        Pout<< "world rank: " << UPstream::myProcNo(UPstream::commWorld())
            << " host-leader rank: "
            << UPstream::myProcNo(UPstream::commInterHost())
            << " intra-host rank: "
            << UPstream::myProcNo(UPstream::commIntraHost())
            << endl;

        label value1(0), value2(0), value3(0);
        label hostIndex = UPstream::myProcNo(commInterHost);

        if (UPstream::master(commInterHost))
        {
            value1 = 100;
            value2 = 200;
        }
        if (UPstream::master(commIntraHost))
        {
            value3 = 300;
        }

        Pstream::broadcast(value1, commInterHost);
        Pstream::broadcast(value2, commIntraHost);
        Pstream::broadcast(hostIndex, commIntraHost);

        Pout<< "host: " << hostIndex
            << " broadcast 1: "
            << value1 << ' '
            << value2 << ' '
            << value3 << endl;

        // re-broadcast
        Pstream::broadcast(value1, commIntraHost);
        Pout<< "host: " << hostIndex
            << " broadcast 2: "
            << value1 << endl;


        label reduced1 = value1;
        label reduced2 = value1;

        reduce
        (
            reduced1,
            sumOp<label>(),
            UPstream::msgType(),
            commIntraHost
        );

        reduce
        (
            reduced2,
            sumOp<label>(),
            UPstream::msgType(),
            commInterHost
        );

        Pout<< "value1: (host) " << reduced1
            << " (leader) " << reduced2 << endl;

        // Pout<< "ranks: " << UPstream::nProcs(commInterHost) << endl;

        wordList strings;
        if (UPstream::is_rank(commInterHost))
        {
            strings.resize(UPstream::nProcs(commInterHost));
            strings[UPstream::myProcNo(commInterHost)] = name(pid());
        }

        // Some basic gather/scatter
        Pstream::allGatherList(strings, UPstream::msgType(), commInterHost);

        Pout<< "pids " << flatOutput(strings) << endl;

        Foam::reverse(strings);

        Pstream::broadcast(strings, commIntraHost);
        Pout<< "PIDS " << flatOutput(strings) << endl;
    }


    if (UPstream::parRun() && generalTest)
    {
        #if 1
        // With first ranks
        labelList subRanks =
            identity(UPstream::nProcs(UPstream::commWorld()) / 2);

        UPstream::communicator newComm;

        newComm.reset(UPstream::commWorld(), subRanks);
        label localRanki = UPstream::myProcNo(newComm);

        const int myProci = UPstream::myProcNo(UPstream::commWorld());

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " a:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif

        #if 1
        // With every other rank
        subRanks = identity(UPstream::nProcs(UPstream::commWorld()));

        for (label& val : subRanks)
        {
            if (val % 2) val = -1;
        }

        newComm.reset(UPstream::commWorld(), subRanks);
        localRanki = UPstream::myProcNo(newComm);

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " b:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
