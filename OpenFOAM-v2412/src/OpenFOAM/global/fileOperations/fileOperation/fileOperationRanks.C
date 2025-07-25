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

\*---------------------------------------------------------------------------*/

#include "fileOperation.H"
#include "stringOps.H"
#include "ITstream.H"
#include "Pstream.H"
#include "SHA1.H"
#include "OSspecific.H"  // for hostName()
#include <cinttypes>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

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

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::labelRange Foam::fileOperation::subRanks(const labelUList& mainIOranks)
{
    // Fast path - no IO ranks.
    if (mainIOranks.empty())
    {
        return labelRange();
    }

    // The lowest numbered rank is the IO rank
    // - linear search for the enclosing range
    // - fallback is proc = 0, which silently adds master (0) into IO ranks

    label begProc = 0;
    label endProc = UPstream::nProcs(UPstream::worldComm);

    const label myProci = UPstream::myProcNo(UPstream::worldComm);

    forAllReverse(mainIOranks, i)
    {
        if (mainIOranks[i] <= myProci)
        {
            begProc = mainIOranks[i];
            if (i+1 < mainIOranks.size())
            {
                endProc = mainIOranks[i+1];
            }
            break;
        }
    }

    return labelRange(begProc, (endProc-begProc));
}


Foam::labelList Foam::fileOperation::getGlobalHostIORanks()
{
    const label numProcs = UPstream::nProcs(UPstream::worldComm);

    // Use hostname
    // Lowest rank per hostname is the IO rank

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

    labelList ranks;
    DynamicList<label> dynRanks;

    if (UPstream::master(UPstream::worldComm))
    {
        dynRanks.reserve(numProcs);

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
    return ranks;
}


Foam::labelList Foam::fileOperation::getGlobalIORanks
(
    // const bool useHost
)
{
    // bool byHostName = useHost;
    bool byHostName = false;

    DynamicList<label> dynRanks;

    Foam::string str(Foam::getEnv("FOAM_IORANKS"));

    if (!str.empty())
    {
        if (str.contains('('))
        {
            // Looks like a list - tokenise it
            ITstream is(str);
            if (!is.empty())
            {
                is >> dynRanks;
            }
        }
        else if (str == "host")
        {
            // Select by hostname
            byHostName = true;
        }
        else
        {
            // Manual parse
            dynRanks = splitStringToList<label>(str);
        }
    }

    if (dynRanks.size())
    {
        if (!dynRanks.contains(0))
        {
            // Could also add silently
            // dynRanks.push_back(0);
            FatalErrorInFunction
                << "Rank 0 (master) should be in the IO ranks. Currently:" << nl
                << "    " << flatOutput(dynRanks) << nl
                << exit(FatalError);
        }

        // Never trust user input.
        // Sort and eliminate any duplicates

        std::sort(dynRanks.begin(), dynRanks.end());

        if (dynRanks.front() < 0)
        {
            FatalErrorInFunction
                << "Cannot have negative ranks! Currently:" << nl
                << "    " << flatOutput(dynRanks) << nl
                << exit(FatalError);
        }

        labelList ranks;

        auto last = std::unique(dynRanks.begin(), dynRanks.end());

        if (last < dynRanks.end())
        {
            ranks = dynRanks.slice(0, (last - dynRanks.begin()));
        }
        else
        {
            ranks = dynRanks;
        }

        return ranks;
    }
    else if (byHostName)
    {
        return fileOperation::getGlobalHostIORanks();
    }

    return labelList();
}


bool Foam::fileOperation::isIOrank(const label proci) const
{
    return
    (
        UPstream::parRun()
      ? UPstream::master(comm_)
      : ioRanks_.empty()
      ? (proci == 0)                // No io-ranks, assume single communicator
      : ioRanks_.contains(proci)    // Found proci in IO rank
    );
}


void Foam::fileOperation::printRanks() const
{
    // Collect the names of the IO masters
    stringList hosts(UPstream::nProcs(UPstream::worldComm));
    if (UPstream::master(comm_))
    {
        hosts[UPstream::myProcNo(UPstream::worldComm)] = hostName();
    }
    Pstream::gatherList(hosts, UPstream::msgType(), UPstream::worldComm);


    DynamicList<label> offsetMaster;

    // Calculate the offsets/counts

    if (UPstream::master(UPstream::worldComm))
    {
        label nHostRanks = 0;
        forAll(hosts, ranki)
        {
            if (!hosts[ranki].empty())
            {
                ++nHostRanks;
            }
        }
        offsetMaster.reserve(nHostRanks+1);

        forAll(hosts, ranki)
        {
            if (!hosts[ranki].empty())
            {
                offsetMaster.push_back(ranki);
            }
        }

        // End of range is nProcs
        offsetMaster.push_back(hosts.size());
    }

    if (offsetMaster.size() > 2)
    {
        DetailInfo
            << "I/O on :" << nl << '(' << nl;
        for (label group = 1; group < offsetMaster.size(); ++group)
        {
            const label beg = offsetMaster[group-1];
            const label end = offsetMaster[group];

            DetailInfo
                << "    (" << hosts[beg].c_str() << ' '
                << (end-beg) << ')' << nl;
        }
        DetailInfo
            << ')' << nl;
    }
}


// ************************************************************************* //
