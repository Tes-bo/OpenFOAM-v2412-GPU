/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-parallel-nbx2

Description
    Test for send/receive data

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addBoolOption("non-blocking", "Test with non-blocking receives");

    #include "setRootCase.H"

    const bool optNonBlocking = args.found("non-blocking");

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    Info<< "\nTesting with non-blocking receives: " << optNonBlocking << nl;


    const int tag = (UPstream::msgType() + 314159);
    const label comm = UPstream::worldComm;

    Random rnd(20*UPstream::myProcNo());

    // Looks a bit like a DIY PstreamBuffers...
    Map<DynamicList<char>> sendBufs;
    Map<DynamicList<char>> recvBufs;

    DynamicList<UPstream::Request> sendRequests(10);
    DynamicList<UPstream::Request> recvRequests(10);

    if (!Pstream::master())
    {
        // Send some random length to master

        const int toProci = UPstream::masterNo();

        label len = rnd.position<label>(10, 20);
        if (UPstream::myProcNo() && (UPstream::myProcNo() % 3) == 0) len = 0;

        scalarField fld(len, scalar(UPstream::myProcNo()));

        // Format for sending
        if (!fld.empty())
        {
            auto& buf = sendBufs(toProci);
            UOPstream os(buf);
            os << fld;
        }

        // Start nonblocking synchronous send to process dest

        if (sendBufs.found(toProci) && !sendBufs[toProci].empty())
        {
            Pout<< "send: [" << sendBufs[toProci].size() << " bytes] "
                << flatOutput(fld) << endl;

            // Has data to send
            UOPstream::write
            (
                sendRequests.emplace_back(),
                UPstream::masterNo(),
                sendBufs[toProci],
                tag,
                comm,
                UPstream::sendModes::sync
            );
        }
    }


    // Probe and receive

    UPstream::Request barrierReq;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        std::pair<int, int64_t> probed =
            UPstream::probeMessage
            (
                UPstream::commsTypes::nonBlocking,
                -1,  // ANY_SOURCE
                tag,
                comm
            );

        if (probed.second > 0)
        {
            // Message found and had size: receive it

            const label proci(probed.first);
            const label count(probed.second);

            if (optNonBlocking)
            {
                recvBufs(proci).resize_nocopy(count);

                // Non-blocking read - MPI_Irecv()
                UIPstream::read
                (
                    recvRequests.emplace_back(),
                    proci,
                    recvBufs[proci],
                    tag,
                    comm
                );
                // Pout<< "Done: "
                //     << UPstream::finishedRequests(recvRequests) << endl;
            }
            else
            {
                IPstream is
                (
                    UPstream::commsTypes::scheduled,  // ie, MPI_Recv()
                    proci,
                    count,  // bufSize
                    tag,
                    comm
                );

                scalarField fld(is);

                Info<< "from [" << probed.first
                    << "] : " << flatOutput(fld) << endl;
            }
        }

        if (barrier_active)
        {
            // Test barrier for completion
            if (UPstream::finishedRequest(barrierReq))
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            if (UPstream::finishedRequests(sendRequests))
            {
                UPstream::barrier(comm, &barrierReq);
                barrier_active = true;
            }
        }
    }

    Pout<< "pending receives: " << recvRequests.size() << endl;

    // Wait for receives to complete
    UPstream::waitRequests(recvRequests);

    // It could be we need this type of synchronization point
    // if the receives are non-blocking
    if (optNonBlocking)
    {
        UPstream::barrier(comm);
    }

    if (!recvBufs.empty())
    {
        Pout<< "Receives from: " << flatOutput(recvBufs.sortedToc()) << endl;

        forAllConstIters(recvBufs, iter)
        {
            Pout<< "proc:" << iter.key() << " len:" << iter.val().size() << nl;

            if (!iter.val().empty())
            {
                UIPstream is(iter.val());
                scalarField fld(is);

                Pout<< "recv:" << iter.key()
                    << " : " << flatOutput(fld) << nl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
