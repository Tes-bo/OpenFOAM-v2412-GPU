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
    Test-parallel-comm3b

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
#include "StringStream.H"
#include "Random.H"

using namespace Foam;


void printRequests(const UList<UPstream::Request>& requests)
{
    OStringStream buf;

    buf << "request: " << requests.size() << '(';

    for (const auto& req : requests)
    {
        if (req.good())
        {
            buf << " " << Foam::name(req.pointer());
        }
        else
        {
            buf << " null";
        }
    }

    buf << " )";
    Pout << buf.str().c_str() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    const int tag = (UPstream::msgType() + 314159);
    const label comm = UPstream::worldComm;

    Random rnd(20*UPstream::myProcNo());

    Map<DynamicList<char>> sendBufs;
    Map<DynamicList<char>> recvBufs;

    DynamicList<UPstream::Request> sendRequests(10);
    DynamicList<UPstream::Request> recvRequests(10);

    // Map request indices to procs
    Map<label> recvFromProc(20);

    if (!Pstream::master())
    {
        // Send some random length to master

        const int toProci = UPstream::masterNo();

        label len = rnd.position<label>(10, 20);
        if (UPstream::myProcNo() && (UPstream::myProcNo() % 3) == 0) len = 0;

        // Has data to send
        if (len)
        {
            auto& buf = sendBufs(toProci);
            buf.resize(len, 'x');

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

            recvBufs(proci).resize_nocopy(count);
            recvFromProc(recvRequests.size()) = proci;

            // Non-blocking read
            UIPstream::read
            (
                recvRequests.emplace_back(),
                proci,
                recvBufs[proci],
                tag,
                comm
            );
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

    if (recvRequests.empty())
    {
        Pout << "No receive requests" << endl;
    }
    else
    {
        printRequests(recvRequests);
    }


    // Either MPI_Waitall, or MPI_Waitany...

    label loop = 0;
    for (bool dispatched = recvRequests.empty(); !dispatched; /*nil*/)
    {
        label index = UPstream::waitAnyRequest(recvRequests);

        if (index < 0)
        {
            Pout<< "Waitany (loop:" << loop << ") : done" << endl;
            dispatched = true;
        }
        else
        {
            Pout<< "Waitany (loop:"
                << loop << ") "
                << index << " of " << recvRequests.size()
                << " from proc:" << recvFromProc.lookup(index, -1)
                << endl;

            printRequests(recvRequests);
        }
        ++loop;
    }

    // Not needed: all tested...
    // UPstream::waitRequest(recvRequests);

    UPstream::barrier(UPstream::worldComm);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
