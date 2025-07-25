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
    Test-parallel-waitSome

Description
    Test polling versus wait-all for processing receive data.
    Will not see much difference between -wait-all and -no-polling though
    since the master doesn't have enough other work.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"
#include "Switch.H"
#include "clockTime.H"

using namespace Foam;


// The 'classic' waiting receive, but also only waiting for recv request
template<class Type>
void waitingReceive
(
    const labelRange& recvRequests,
    const List<List<Type>>& recvBuffers,
    const bool waitAll = false
)
{
    clockTime waitTiming;

    if (waitAll)
    {
        // Wait for send and recv (assumes recv followed by send)
        UPstream::waitRequests(recvRequests.start(), -1);
    }
    else
    {
        // Wait for receives only
        UPstream::waitRequests(recvRequests.start(), recvRequests.size());
    }

    double waited = waitTiming.timeIncrement();
    if (waited > 1e-3)
    {
        Pout<< "waited: " << waited << " before processing" << endl;
    }

    forAll(recvBuffers, proci)
    {
        const auto& slice = recvBuffers[proci];

        if (!slice.empty())
        {
            // Process data from proci
            Pout<< "proc:" << proci
                << ' ' << flatOutput(slice) << nl;
        }
    }
}


// Polling receive
template<class Type>
void pollingReceive
(
    const labelRange& recvRequests,
    const UList<int>& recvProcs,
    const List<List<Type>>& recvBuffers
)
{
    clockTime waitTiming;

    DynamicList<int> indices(recvRequests.size());

    if (!recvRequests.empty()) Pout<< "..." << endl;

    for
    (
        label loop = 0;
        UPstream::waitSomeRequests
        (
            recvRequests.start(),
            recvRequests.size(),
            &indices
        );
        ++loop
    )
    {
        double waited = waitTiming.timeIncrement();
        if (waited <= 1e-3)
        {
            waited = 0;
        }
        Pout<< "loop:" << loop
            << " waited: " << waited
            << " before processing" << endl;

        for (const int idx : indices)
        {
            const int proci = recvProcs[idx];
            const auto& slice = recvBuffers[proci];

            // Process data from proci
            Pout<< "loop:" << loop << " polled:" << indices.size()
                << " proc:" << proci
                << ' ' << flatOutput(slice) << endl;
        }
        Pout<< "..." << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("timings etc");
    argList::addBoolOption("no-polling", "wait all instead of polling");
    argList::addBoolOption("wait-all", "wait all instead of polling");
    argList::addOption("sleep", "s", "change sleep (default: 5)");
    argList::noCheckProcessorDirectories();

    const label transferSize = 10;
    label sleepSeconds = 5;

    #include "setRootCase.H"

    args.readIfPresent("sleep", sleepSeconds);
    const bool waitAll = args.found("wait-all");
    const bool nonPolling = args.found("no-polling");

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    Info<< "Calling with sleep=" << sleepSeconds
        << ", polling=" << Switch::name(!nonPolling)
        << ", wait-all=" << Switch::name(waitAll) << nl;

    labelList sendBuffer;
    List<labelList> recvBuffers;

    if (UPstream::master())
    {
        recvBuffers.resize(UPstream::nProcs());
    }
    else
    {
        recvBuffers.resize(1);
    }

    clockTime timing;

    const label startOfRequests = UPstream::nRequests();

    // Setup receives
    labelRange recvRequests(UPstream::nRequests(), 0);
    DynamicList<int> recvProcs(UPstream::nProcs());

    if (UPstream::master())
    {
        for (const int proci : UPstream::subProcs())
        {
            // The rank corresponding to the request
            recvProcs.push_back(proci);
            auto& slice = recvBuffers[proci];
            slice.resize_nocopy(transferSize);

            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slice
            );
        }
    }
    else
    {
        const int proci = UPstream::masterNo();

        if ((UPstream::myProcNo() % 2) == 0)
        {
            recvProcs.push_back(proci);
            auto& slice = recvBuffers[proci];
            slice.resize_nocopy(transferSize);

            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slice
            );
        }
    }
    // OR: recvRequests.size() = (UPstream::nRequests() - recvRequests.start());
    recvRequests += recvProcs.size();


    labelList overallRecvRequests
    (
        UPstream::listGatherValues<label>(recvRequests.size())
    );

    Info<< "Number of recv requests: "
        << flatOutput(overallRecvRequests) << nl << nl;


    // Setup sends
    sendBuffer.resize_nocopy(transferSize);
    sendBuffer = UPstream::myProcNo();

    const auto startBufferSend = [&]() -> void
    {
        if (sleepSeconds > 0)
        {
            // Dispatch some immediately, others with a delay
            if ((UPstream::myProcNo() % 2) == 0)
            {
                sleep(sleepSeconds);
            }
            else if ((UPstream::myProcNo() % 3) == 0)
            {
                sleep(1.5*sleepSeconds);
            }
        }

        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            sendBuffer
        );
    };


    if (UPstream::master())
    {
        for (const int proci : UPstream::subProcs())
        {
            if ((UPstream::myProcNo() % 2) == 0)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBuffer
                );
            }
        }
    }
    else if (waitAll)
    {
        startBufferSend();
    }


    // Some skulduggery to get a differential in timings...

    const int nloops = (UPstream::master() ? 1 : 2);

    for (int loopi = 0; loopi < nloops; ++loopi)
    {
        if (waitAll || nonPolling)
        {
            waitingReceive(recvRequests, recvBuffers, waitAll);
        }
        else
        {
            pollingReceive(recvRequests, recvProcs, recvBuffers);
        }

        // Timing for processing all the receives
        if (args.verbose())
        {
            Pout<< "receive: " << timing.timeIncrement() << 's' << endl;
        }

        if (!UPstream::master() && loopi == 0 && !waitAll)
        {
            startBufferSend();
        }
    }

    if (args.verbose())
    {
        Pout<< "timing: " << timing.elapsedTime() << 's' << endl;
    }

    // Final
    UPstream::waitRequests(startOfRequests);

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
