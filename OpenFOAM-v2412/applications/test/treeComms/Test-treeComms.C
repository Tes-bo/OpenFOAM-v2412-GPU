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

Description
    Print/test tree communication patterns

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "Time.H"
#include "vector.H"
#include "Fstream.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

using namespace Foam;

void printConnection(Ostream& os, const label proci, const labelUList& below)
{
    for (const label connectProci : below)
    {
        os << indent << proci << " -- " << connectProci << nl;
    }
}


// The number of receives - as per gatherList (v2112)
void printRecvCount_gatherList
(
    const UList<UPstream::commsStruct>& comms,
    const label comm = UPstream::worldComm
)
{
    const label np = UPstream::nProcs(comm);
    labelList nMesg;
    labelList nRecv;

    if (UPstream::parRun() && np > 1 && UPstream::master(comm))
    {
        nMesg.resize(np, Zero);
        nRecv.resize(np, Zero);

        forAll(comms, proci)
        {
            nMesg[proci] += comms[proci].below().size();

            // Receive from my downstairs neighbours
            for (const label belowID : comms[proci].below())
            {
                const labelList& belowLeaves = comms[belowID].allBelow();
                nRecv[proci] += (belowLeaves.size() + 1);
            }
        }
    }

    label nTotalMesg = sum(nMesg);
    label nTotalRecv = sum(nRecv);

    Info<< "gatherList communication:" << nl
        << "messages = " << nTotalMesg << " -> "
        << flatOutput(nMesg) << nl
        << "gathers = " << nTotalRecv << " -> "
        << flatOutput(nRecv) << nl;
}


// The number of sends - as per scatterList (v2112)
void printSendCount_scatterList
(
    const UList<UPstream::commsStruct>& comms,
    const label comm = UPstream::worldComm
)
{
    const label np = UPstream::nProcs(comm);
    labelList nMesg;
    labelList nSend;

    if (UPstream::parRun() && np > 1 && UPstream::master(comm))
    {
        nMesg.resize(np, Zero);
        nSend.resize(np, Zero);

        forAll(comms, proci)
        {
            nMesg[proci] += comms[proci].below().size();

            // Send to my downstairs neighbours
            for (const label belowID : comms[proci].below())
            {
                const labelList& notBelowLeaves = comms[belowID].allNotBelow();
                nSend[proci] += notBelowLeaves.size();
            }
        }
    }

    label nTotalMesg = sum(nMesg);
    label nTotalSend = sum(nSend);

    Info<< "scatterList communication:" << nl
        << "messages = " << nTotalMesg << " -> "
        << flatOutput(nMesg) << nl
        << "scatters = " << nTotalSend << " -> "
        << flatOutput(nSend) << nl;
}


// Transmission widths (contiguous data)
void printWidths
(
    const UList<UPstream::commsStruct>& comms,
    const label comm = UPstream::worldComm
)
{
    const label np = UPstream::nProcs(comm);
    labelList maxBelow;
    labelList maxNotBelow;

    if (UPstream::parRun() && np > 1 && UPstream::master(comm))
    {
        maxBelow.resize(np, Zero);
        maxNotBelow.resize(np, Zero);

        forAll(comms, proci)
        {
            label& max0 = maxBelow[proci];
            label& max1 = maxNotBelow[proci];

            for (const label belowID : comms[proci].below())
            {
                // Receive from my downstairs neighbours
                const labelList& belowLeaves = comms[belowID].allBelow();

                // Send to my downstairs neighbours
                const labelList& notBelowLeaves = comms[belowID].allNotBelow();

                max0 = max(max0, belowLeaves.size());
                max1 = max(max1, notBelowLeaves.size());
            }
        }
    }

    Info<< "buffer width:" << nl
        << "gathers  -> " << flatOutput(maxBelow) << nl
        << "scatters -> " << flatOutput(maxNotBelow) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"
    #include "createTime.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Please run in parallel" << exit(FatalError);
    }
    const auto& comms = UPstream::treeCommunication();

    printRecvCount_gatherList(comms);
    printSendCount_scatterList(comms);
    printWidths(comms);

    // My communication order
    const UPstream::commsStruct& myComm = comms[UPstream::myProcNo()];

    // Info<< "allComms: " << comms << nl;

    if (UPstream::master())
    {
        OFstream os("treeComm.dot");

        os << "// tree communication" << nl << nl;
        os.beginBlock("graph");

        printConnection(os, 0, myComm.below());
        // Pout<< flatOutput(myComm.allBelow()) << nl;

        for (const int proci : UPstream::subProcs())
        {
            labelList below;
            IPstream::recv(below, proci);

            printConnection(os, proci, below);
        }

        os.endBlock();

        Info<< "Wrote processorTopology graph: "
            << runTime.relativePath(os.name()) << nl;

        Info<< nl
            << "Use dot (or other graphviz tools)" << nl;
    }
    else
    {
        OPstream::send(myComm.below(), UPstream::masterNo());
        // Pout<< flatOutput(myComm.allBelow()) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
