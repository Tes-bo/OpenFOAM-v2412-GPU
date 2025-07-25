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
    Test-parallel-comm3a

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

#include "openfoam_mpi.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void printRequests(const UList<MPI_Request>& requests)
{
    OStringStream buf;

    buf << "request: " << requests.size() << '(';

    for (const auto& req : requests)
    {
        if (req == MPI_REQUEST_NULL)
        {
            buf << " null";
        }
        else
        {
            buf << " " << Foam::name(req);
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
    // const label comm = UPstream::worldComm;

    Random rnd(20*UPstream::myProcNo());

    Map<DynamicList<char>> sendBufs;
    Map<DynamicList<char>> recvBufs;

    DynamicList<MPI_Request> sendRequests(10);
    DynamicList<MPI_Request> recvRequests(10);


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

            MPI_Issend
            (
                buf.cdata_bytes(),
                buf.size_bytes(),
                MPI_BYTE,
                toProci,
                tag,
                MPI_COMM_WORLD,
               &sendRequests.emplace_back()
            );
        }
    }


    // Probe and receive

    MPI_Request barrierReq;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        int flag = 0;
        MPI_Message message;
        MPI_Status status;

        MPI_Improbe
        (
            MPI_ANY_SOURCE,
            tag,
            MPI_COMM_WORLD,
           &flag,
           &message,
           &status
        );

        if (flag)
        {
            // Message found, receive into dest buffer location
            const label fromProci = status.MPI_SOURCE;

            int count = 0;
            MPI_Get_count(&status, MPI_BYTE, &count);

            auto& buf = recvBufs(fromProci);
            buf.resize_nocopy(count);

            MPI_Imrecv
            (
                buf.data_bytes(),
                buf.size_bytes(),
                MPI_BYTE,
               &message,
               &recvRequests.emplace_back()
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            MPI_Test(&barrierReq, &flag, MPI_STATUS_IGNORE);

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
                sendRequests.size(), sendRequests.data(),
                &flag, MPI_STATUSES_IGNORE
            );

            if (flag)
            {
                MPI_Ibarrier(MPI_COMM_WORLD, &barrierReq);
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
        int index = 0;
        MPI_Waitany
        (
            recvRequests.size(),
            recvRequests.data(),
            &index,
            MPI_STATUS_IGNORE
        );

        if (index == MPI_UNDEFINED)
        {
            //Pout<< "Testany is " << (flag ? "done" : "waiting") << endl;
            Pout<< "Waitany (loop:" << loop << ") : done" << endl;
            dispatched = true;
        }
        else
        {
            Pout<< "Waitany (loop:"
                << loop << ") "
                << index << " of " << recvRequests.size() << endl;

            printRequests(recvRequests);
        }

        ++loop;
    }

    // Not needed: all tested...
    // MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
