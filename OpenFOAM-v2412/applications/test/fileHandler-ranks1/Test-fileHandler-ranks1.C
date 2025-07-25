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
    Test-fileHandler-ranks1

Description
    Test IO ranks and ranks selection

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "fileOperation.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "SHA1.H"
#include "stringOps.H"

using namespace Foam;

// Parse space, comma, semicolon separated list of integers, floats etc...
template<class PrimitiveType>
static List<PrimitiveType> splitStringToList(const std::string& str)
{
    const auto items = stringOps::splitAny(str, " ,;");

    DynamicList<PrimitiveType> values(items.size());

    for (const auto& item : items)
    {
        const std::string s(item.str());

        PrimitiveType val;

        if (Foam::read(s, val))
        {
            values.push_back(val);
        }
        else
        {
            // Report errors? Could get noisy...
        }
    }

    return List<PrimitiveType>(std::move(values));
}


//- Construct by parsing string for scalar ranges
//  The individual items are space, comma or semicolon delimited.
static labelList parseIOranks
(
    const Foam::string& str,
    label numProcs = -1
)
{
    bool byHostName = false;

    labelList ranks;

    // Info<< "parsing: " << str << endl;

    if (!str.empty())
    {
        if (str.contains('('))
        {
            // Looks like a list - tokenise it

            ITstream is(str);
            if (!is.empty())
            {
                is >> ranks;
            }
        }
        else if (str == "host")
        {
            // Select by host
            byHostName = true;
        }
        else
        {
            // Manual parse
            ranks = splitStringToList<label>(str);
        }
    }

    if (ranks.size())
    {
        // Never trust user input.
        // Sort and eliminate any duplicates

        std::sort(ranks.begin(), ranks.end());

        auto last = std::unique(ranks.begin(), ranks.end());

        label newLen = label(last - ranks.begin());

        // Detect any values that are too large
        if (numProcs > 0)
        {
            auto iter = std::find_if
            (
                ranks.begin(),
                last,
                [=](label proci) { return proci >= numProcs; }
            );

            if (last != iter)
            {
                newLen = label(iter - ranks.begin());
            }
        }

        ranks.resize(newLen);
    }
    else if (byHostName)
    {
        // Use hostname
        // Lowest rank per hostname is the IO rank

        numProcs = UPstream::nProcs(UPstream::worldComm);

        List<SHA1Digest> digests;
        if (UPstream::master(UPstream::worldComm))
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
            SHA1Digest::max_size(),     // Num send/recv per rank
            UPstream::worldComm
        );

        if (UPstream::master(UPstream::worldComm))
        {
            DynamicList<label> dynRanks(numProcs);

            dynRanks.push_back(0);  // Always include master
            label previ = 0;

            for (label proci = 1; proci < digests.size(); ++proci)
            {
                if (digests[proci] != digests[previ])
                {
                    dynRanks.push_back(proci);
                    previ = proci;
                }
            }

            ranks.transfer(dynRanks);
        }

        Pstream::broadcast(ranks, UPstream::worldComm);
    }

    return ranks;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    argList::addOption("io-ranks", "list", "shadow for -ioRanks (testing)");
    argList::addOption("pick", "list", "limited subset of procs");

    labelList ioRanks;

    // Pre-check
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-io-ranks") == 0)
        {
            if (argi < argc-1)
            {
                ioRanks = parseIOranks(argv[argi+1]);
                break;
            }
        }
    }

    #include "setRootCase.H"

    bitSet useProc;

    if (args.found("pick"))
    {
        useProc = bitSet
        (
            UPstream::nProcs(),
            parseIOranks(args.get<string>("pick"))
        );
    }

    Info<< "procs: " << UPstream::nProcs() << endl;
    Info<< "io-ranks: " << flatOutput(ioRanks) << endl;
    Info<< "-ioRanks: "
        << flatOutput(fileOperation::getGlobalIORanks()) << endl;

    Info<< "pick: " << flatOutput(useProc.toc()) << endl;

    // labelList subRanks = fileOperation::getGlobalSubRanks(useProc);
    //
    // Pout<< "sub ranks: "
    //     << flatOutput(fileOperation::getGlobalSubRanks(useProc)) << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
