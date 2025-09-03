/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "mapDistributeBase.H"
#include "bitSet.H"
#include "commSchedule.H"
#include "labelPairHashes.H"
#include "globalIndex.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapDistributeBase, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::mapDistributeBase::hasFlipAddressing(const labelUList& map)
{
    for (const label val : map)
    {
        if (!val)
        {
            // Cannot be flipped addressing if it contains zero.
            return false;
        }
        else if (val < 0)
        {
            // Must be flipped addressing if it contains negatives.
            return true;
        }
    }

    return false;
}


bool Foam::mapDistributeBase::hasFlipAddressing(const labelListList& maps)
{
    for (const labelList& map : maps)
    {
        for (const label val : map)
        {
            if (!val)
            {
                // Cannot be flipped addressing if it contains zero.
                return false;
            }
            else if (val < 0)
            {
                // Must be flipped addressing if it contains negatives.
                return true;
            }
        }
    }

    return false;
}


Foam::label Foam::mapDistributeBase::getMappedSize
(
    const labelListList& maps,
    const bool hasFlip
)
{
    label maxIndex = -1;

    for (const labelList& map : maps)
    {
        for (label index : map)
        {
            if (hasFlip)
            {
                index = mag(index)-1;
            }

            maxIndex = max(maxIndex, index);
        }
    }

    return (maxIndex+1);
}


Foam::label Foam::mapDistributeBase::countUnmapped
(
    const labelUList& elements,
    const labelListList& maps,
    const bool hasFlip
)
{
    if (elements.empty())
    {
        return 0;
    }

    // Moderately efficient markup/search

    bitSet unvisited(elements);
    label nUnmapped = unvisited.count();

    if (hasFlip)
    {
        for (const labelList& map : maps)
        {
            for (label index : map)
            {
                index = mag(index)-1;

                if (unvisited.unset(index))
                {
                    --nUnmapped;
                    if (!nUnmapped) break;
                }
            }
        }
    }
    else
    {
        for (const labelList& map : maps)
        {
            for (label index : map)
            {
                if (unvisited.unset(index))
                {
                    --nUnmapped;
                    if (!nUnmapped) break;
                }
            }
        }
    }

    return nUnmapped;
}


void Foam::mapDistributeBase::checkReceivedSize
(
    const label proci,
    const label expectedSize,
    const label receivedSize
)
{
    if (receivedSize != expectedSize)
    {
        FatalErrorInFunction
            << "From processor " << proci
            << " : expected " << expectedSize
            << " but received " << receivedSize << " elements" << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::mapDistributeBase::schedule
(
    const labelListList& subMap,
    const labelListList& constructMap,
    const int tag,
    const label comm
)
{
    const label myRank = UPstream::myProcNo(comm);
    const label nProcs = UPstream::nProcs(comm);

    // Communications: send and receive processor
    List<labelPair> allComms;

    {
        labelPairHashSet commsSet(nProcs);

        // Find what communication is required
        forAll(subMap, proci)
        {
            if (proci != myRank)
            {
                if (subMap[proci].size())
                {
                    // I need to send to proci
                    commsSet.insert(labelPair(myRank, proci));
                }
                if (constructMap[proci].size())
                {
                    // I need to receive from proci
                    commsSet.insert(labelPair(proci, myRank));
                }
            }
        }
        allComms = commsSet.toc();
    }


    // Gather/reduce
    if (UPstream::master(comm))
    {
        // Receive and merge
        for (const int proci : UPstream::subProcs(comm))
        {
            List<labelPair> nbrData;
            IPstream::recv(nbrData, proci, tag, comm);

            for (const labelPair& connection : nbrData)
            {
                allComms.push_uniq(connection);
            }
        }
    }
    else
    {
        if (UPstream::parRun())
        {
            OPstream::send(allComms, UPstream::masterNo(), tag, comm);
        }
    }

    // Broadcast: send comms information to all
    Pstream::broadcast(allComms, comm);

    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            nProcs,
            allComms
        ).procSchedule()[myRank]
    );

    // Processors involved in my schedule
    return List<labelPair>(allComms, mySchedule);
}


const Foam::List<Foam::labelPair>& Foam::mapDistributeBase::schedule() const
{
    if (!schedulePtr_)
    {
        schedulePtr_.reset
        (
            new List<labelPair>
            (
                schedule(subMap_, constructMap_, UPstream::msgType(), comm_)
            )
        );
    }

    return *schedulePtr_;
}


const Foam::List<Foam::labelPair>& Foam::mapDistributeBase::whichSchedule
(
    const UPstream::commsTypes commsType
) const
{
    if (commsType == UPstream::commsTypes::scheduled)
    {
        return schedule();
    }

    return List<labelPair>::null();
}


void Foam::mapDistributeBase::printLayout(Ostream& os) const
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // Determine offsets of remote data.
    labelList minIndex(nProcs, labelMax);
    labelList maxIndex(nProcs, labelMin);
    forAll(constructMap_, proci)
    {
        const labelList& construct = constructMap_[proci];
        if (constructHasFlip_)
        {
            forAll(construct, i)
            {
                label index = mag(construct[i])-1;
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
        else
        {
            forAll(construct, i)
            {
                label index = construct[i];
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
    }

    label localSize(0);

    if (maxIndex[myRank] != labelMin)
    {
        localSize = maxIndex[myRank]+1;
    }

    os  << "Layout: (constructSize:" << constructSize_
        << " subHasFlip:" << subHasFlip_
        << " constructHasFlip:" << constructHasFlip_
        << ")" << nl
        << "local (processor " << myRank << "):" << nl
        << "    start : 0" << nl
        << "    size  : " << localSize << endl;

    label offset = localSize;
    forAll(minIndex, proci)
    {
        if (proci != myRank && !constructMap_[proci].empty())
        {
            label size(0);

            if (maxIndex[proci] != labelMin)
            {
                size = maxIndex[proci]-minIndex[proci]+1;
                if (minIndex[proci] != offset)
                {
                    FatalErrorInFunction
                        << "offset:" << offset
                        << " proci:" << proci
                        << " minIndex:" << minIndex[proci]
                        << abort(FatalError);
                }
            }

            os  << "processor " << proci << ':' << nl
                << "    start : " << offset << nl
                << "    size  : " << size << endl;

            offset += size;
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelUList& elements,
    List<Map<label>>& compactMap
) const
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(nProcs, Zero);

    for (const label globalIdx : elements)
    {
        if (globalIdx != -1 && !globalNumbering.isLocal(myRank, globalIdx))
        {
            label proci = globalNumbering.whichProcID(myRank, globalIdx);
            nNonLocal[proci]++;
        }
    }

    compactMap.resize_nocopy(nProcs);

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != myRank)
        {
            compactMap[proci].reserve(nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    for (const label globalIdx : elements)
    {
        if (globalIdx != -1 && !globalNumbering.isLocal(myRank, globalIdx))
        {
            label proci = globalNumbering.whichProcID(myRank, globalIdx);
            label index = globalNumbering.toLocal(proci, globalIdx);
            compactMap[proci].insert(index, compactMap[proci].size());
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelListList& cellCells,
    List<Map<label>>& compactMap
) const
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(nProcs, Zero);

    for (const labelList& cCells : cellCells)
    {
        for (const label globalIdx : cCells)
        {
            if (globalIdx != -1 && !globalNumbering.isLocal(myRank, globalIdx))
            {
                label proci = globalNumbering.whichProcID(myRank, globalIdx);
                nNonLocal[proci]++;
            }
        }
    }

    compactMap.resize_nocopy(nProcs);

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != myRank)
        {
            compactMap[proci].reserve(nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    for (const labelList& cCells : cellCells)
    {
        for (const label globalIdx : cCells)
        {
            if (globalIdx != -1 && !globalNumbering.isLocal(myRank, globalIdx))
            {
                label proci = globalNumbering.whichProcID(myRank, globalIdx);
                label index = globalNumbering.toLocal(proci, globalIdx);
                compactMap[proci].insert(index, compactMap[proci].size());
            }
        }
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(nProcs);
    compactStart[myRank] = 0;
    constructSize_ = globalNumbering.localSize(myRank);
    forAll(compactStart, proci)
    {
        if (proci != myRank)
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }


    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(nProcs);
    // Compact addressing for received data
    constructMap_.setSize(nProcs);
    forAll(compactMap, proci)
    {
        if (proci == myRank)
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize(myRank);
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIters(compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter.val();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter.val() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(nProcs);
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        comm_
    );

    // Renumber elements
    for (label& elem : elements)
    {
        elem = renumber(globalNumbering, comm_, compactMap, elem);
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(nProcs);
    compactStart[myRank] = 0;
    constructSize_ = globalNumbering.localSize(myRank);
    forAll(compactStart, proci)
    {
        if (proci != myRank)
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }


    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(nProcs);
    // Compact addressing for received data
    constructMap_.setSize(nProcs);
    forAll(compactMap, proci)
    {
        if (proci == myRank)
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize(myRank);
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIters(compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter.val();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter.val() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(nProcs);
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        comm_
    );

    // Renumber elements
    for (labelList& cCells : cellCells)
    {
        for (label& celli : cCells)
        {
            celli = renumber(globalNumbering, comm_, compactMap, celli);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributeBase::mapDistributeBase() noexcept
:
    mapDistributeBase(UPstream::worldComm)
{}


Foam::mapDistributeBase::mapDistributeBase(const label comm) noexcept
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase(const mapDistributeBase& map)
:
    constructSize_(map.constructSize_),
    subMap_(map.subMap_),
    constructMap_(map.constructMap_),
    subHasFlip_(map.subHasFlip_),
    constructHasFlip_(map.constructHasFlip_),
    comm_(map.comm_),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase(mapDistributeBase&& map)
:
    mapDistributeBase(map.comm())
{
    transfer(map);
}


Foam::mapDistributeBase::mapDistributeBase
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    constructSize_(constructSize),
    subMap_(std::move(subMap)),
    constructMap_(std::move(constructMap)),
    subHasFlip_(subHasFlip),
    constructHasFlip_(constructHasFlip),
    comm_(comm),
    schedulePtr_(nullptr)
{}


Foam::mapDistributeBase::mapDistributeBase
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    if (sendProcs.size() != recvProcs.size())
    {
        FatalErrorInFunction
            << "The send and receive data is not the same length. sendProcs:"
            << sendProcs.size() << " recvProcs:" << recvProcs.size()
            << abort(FatalError);
    }

    // Per processor the number of samples we have to send/receive.
    labelList nSend(nProcs, Zero);
    labelList nRecv(nProcs, Zero);

    forAll(sendProcs, sampleI)
    {
        const label sendProc = sendProcs[sampleI];
        const label recvProc = recvProcs[sampleI];

        // Note that also need to include local communication (both
        // RecvProc and sendProc on local processor)

        if (myRank == sendProc)
        {
            // I am the sender.
            nSend[recvProc]++;
        }
        if (myRank == recvProc)
        {
            // I am the receiver.
            nRecv[sendProc]++;
        }
    }

    subMap_.setSize(nProcs);
    constructMap_.setSize(nProcs);
    forAll(nSend, proci)
    {
        subMap_[proci].setSize(nSend[proci]);
        constructMap_[proci].setSize(nRecv[proci]);
    }
    nSend = 0;
    nRecv = 0;

    // Largest entry inside constructMap
    label maxRecvIndex = -1;

    forAll(sendProcs, sampleI)
    {
        const label sendProc = sendProcs[sampleI];
        const label recvProc = recvProcs[sampleI];

        if (myRank == sendProc)
        {
            // I am the sender. Store index I need to send.
            subMap_[recvProc][nSend[recvProc]++] = sampleI;
        }
        if (myRank == recvProc)
        {
            // I am the receiver.
            constructMap_[sendProc][nRecv[sendProc]++] = sampleI;
            maxRecvIndex = sampleI;
        }
    }

    constructSize_ = maxRecvIndex+1;
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    //forAll(compactMap, proci)
    //{
    //    if (proci != myRank)
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        const List<label> sorted(globalMap.sortedToc());
    //
    //        forAll(sorted, i)
    //        {
    //            globalMap(sorted[i]) = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    constructSize_(0),
    subMap_(),
    constructMap_(),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(comm),
    schedulePtr_(nullptr)
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    //forAll(compactMap, proci)
    //{
    //    if (proci != myRank)
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        const List<label> sorted(globalMap.sortedToc());
    //
    //        forAll(sorted, i)
    //        {
    //            globalMap(sorted[i]) = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    const layoutTypes constructLayout,
    labelListList&& subMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    constructSize_(0),
    subMap_(std::move(subMap)),
    constructMap_(),
    subHasFlip_(subHasFlip),
    constructHasFlip_(constructHasFlip),
    comm_(comm),
    schedulePtr_(nullptr)
{
    const label myRank = UPstream::myProcNo(comm_);
    const label nProcs = UPstream::nProcs(comm_);

    // Send over how many i need to receive.
    labelList recvSizes;
    Pstream::exchangeSizes(subMap_, recvSizes, comm_);

    constructSize_ = 0;
    constructMap_.resize(nProcs);

    // The order of receiving:

    if (constructLayout == layoutTypes::linear)
    {
        forAll(constructMap_, proci)
        {
            const label len = recvSizes[proci];

            constructMap_[proci] = identity(len, constructSize_);
            constructSize_ += len;
        }
    }
    else
    {
        // layoutTypes::localFirst

        // My data first
        {
            const label len = recvSizes[myRank];

            constructMap_[myRank] = identity(len, constructSize_);
            constructSize_ += len;
        }

        // What the other processors are sending to me
        forAll(constructMap_, proci)
        {
            if (proci != myRank)
            {
                const label len = recvSizes[proci];

                constructMap_[proci] = identity(len, constructSize_);
                constructSize_ += len;
            }
        }
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    labelListList&& subMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    mapDistributeBase
    (
        layoutTypes::localFirst,
        std::move(subMap),
        subHasFlip,
        constructHasFlip,
        comm
    )
{}


Foam::mapDistributeBase::mapDistributeBase
(
    const UPtrList<const mapDistributeBase>& maps,
    const labelList& localRanks,
    const label newComm,
    const labelListList& newToOldRanks,// from newComm to comm_
    labelList& startOfLocal,
    List<Map<label>>& compactMaps
)
:
    constructSize_(0),
    subHasFlip_(false),
    constructHasFlip_(false),
    comm_(-1),
    schedulePtr_(nullptr)
{
    if (maps.empty())
    {
        return;
    }

    comm_ = newComm;
    subHasFlip_ = maps[0].subHasFlip();
    constructHasFlip_ = maps[0].constructHasFlip();

    const label nNewRanks = newToOldRanks.size();
    const label myNewRank = UPstream::myProcNo(newComm);
    if (nNewRanks != UPstream::nProcs(newComm))
    {
        FatalErrorInFunction<< "nNewRanks:" << nNewRanks
            << " nProcs:" << UPstream::nProcs(newComm)
            << exit(FatalError);
    }

    if (localRanks.size() != maps.size())
    {
        FatalErrorInFunction
            << "Number of maps:" << maps.size()
            << " number of localRanks:" << localRanks.size()
            << exit(FatalError);
    }

    // Sanity checks
    const auto& map0 = maps[0];
    forAll(maps, mapi)
    {
        const auto& map = maps[mapi];

        if
        (
            (map.comm() != map0.comm())
         || (map.subHasFlip() != map0.subHasFlip())
         || (map.constructHasFlip() != map0.constructHasFlip())
        )
        {
            FatalErrorInFunction
                << "Maps should all be the same form"
                << " Map " << mapi
                << " has comm:" << map.comm()
                << " subHasFlip:" << map.subHasFlip()
                << " constructHasFlip:" << map.constructHasFlip()
                << " which is different from map 0"
                << exit(FatalError);
        }

        const label localRank = localRanks[mapi];
        const auto& constructOwn = maps[mapi].constructMap()[localRank];
        forAll(constructOwn, i)
        {
            if (constructOwn[i] != i)
            {
                FatalErrorInFunction
                    << "Maps constructMap not identity."
                    << " Map " << mapi
                    << " constructMap:" << flatOutput(constructOwn)
                    << exit(FatalError);
            }
        }
    }


    constructMap_.resize_nocopy(nNewRanks);
    subMap_.resize_nocopy(nNewRanks);


    // Store starts
    startOfLocal.setSize(maps.size()+1);
    compactMaps.resize_nocopy(maps.size());

    label constructi = 0;
    forAll(maps, mapi)
    {
        startOfLocal[mapi] = constructi;
        const label localRank = localRanks[mapi];
        const auto& map = maps[mapi].constructMap()[localRank];

        // Presize compaction array
        const label nRemote = maps[mapi].constructSize()-map.size();
        compactMaps[mapi].resize(2*nRemote);

        constructi += map.size();
    }
    startOfLocal.last() = constructi;


    // Determine start of constructed remote data. This is used to get the
    // local offset which can then be used to get the relative subMap location.
    labelListList startOfRemote(maps.size());
    forAll(maps, mapi)
    {
        const label nOldProcs = maps[mapi].constructMap().size();
        labelList& starts = startOfRemote[mapi];

        starts.setSize(nOldProcs, labelMax);
        forAll(maps[mapi].constructMap(), oldProci)
        {
            const labelList& map = maps[mapi].constructMap()[oldProci];
            forAll(map, i)
            {
                const label index
                (
                    constructHasFlip_
                  ? mag(map[i])-1
                  : map[i]
                );
                starts[oldProci] = min(starts[oldProci], index);
            }
        }
    }


    // Construct map
    // ~~~~~~~~~~~~~
    // - all localRanks:
    //      - data gets appended in map order
    //      - map is just an offset (startOfLocal)
    // - all previously remote ranks:
    //      - data is already present according to startOfLocal
    //      - map is old-to-new index
    // - all still remote ranks:
    //      - data gets appended in map order after the startOfLocal
    //      - map is old-to-new index


    // Append local (= myRank) data. TBD: assumes subMap and constructMap
    // are identity maps.
    {
        labelList& myConstruct = constructMap_[myNewRank];
        myConstruct.resize_nocopy(constructi);
        constructi = 0;
        forAll(maps, mapi)
        {
            const label localRank = localRanks[mapi];
            const auto& map = maps[mapi].constructMap()[localRank];
            const label offset = startOfLocal[mapi];

            forAll(map, i)
            {
                if (constructHasFlip_)
                {
                    forAll(map, i)
                    {
                        if (map[i] < 0)
                        {
                            myConstruct[constructi++] = map[i]-offset;
                        }
                        else
                        {
                            myConstruct[constructi++] = map[i]+offset;
                        }
                    }
                }
                else
                {
                    myConstruct[constructi++] = map[i]+offset;
                }
            }
        }
    }

    // Filter remote construct data
    {
        // Remote ranks that are now local
        //  - store new index for mapping stencils
        //  - no need to construct since already
        const auto& oldProcs = newToOldRanks[myNewRank];

        forAll(maps, mapi)
        {
            for (const label oldProci : oldProcs)
            {
                if (oldProci != localRanks[mapi])
                {
                    const auto& map = maps[mapi].constructMap()[oldProci];

                    if (!map.size())
                    {
                        continue;
                    }


                    // The slots come from a local map so we can look up the
                    // new location
                    const label sourceMapi = localRanks.find(oldProci);
                    const auto& subMap =
                        maps[sourceMapi].subMap()[localRanks[mapi]];

                    //Pout<< "From oldRank:" << oldProci
                    //    << " sending to masterRank:" << localRanks[mapi]
                    //    << " elements:" << flatOutput(subMap)
                    //    << nl
                    //    << "   received as elements:" << flatOutput(map)
                    //    << endl;

                    if (map.size() != subMap.size())
                    {
                        FatalErrorInFunction << "Problem:"
                            << "oldProci:" << oldProci
                            << " mapi:" << mapi
                            << " constructMap:" << map.size()
                            << " sourceMapi:" << sourceMapi
                            << " subMap:" << subMap.size()
                            << exit(FatalError);
                    }

                    const label offset = startOfLocal[sourceMapi];
                    // Construct map starts after the local data
                    const label nMapLocal = startOfRemote[mapi][oldProci];

                    auto& cptMap = compactMaps[mapi];
                    forAll(map, i)
                    {
                        // old slot position to new slot position
                        const label index
                        (
                            constructHasFlip_
                          ? mag(map[i])-1
                          : map[i]
                        );
                        const label newIndex = subMap[index-nMapLocal]+offset;

                        // Note: should always warn for duplicates? Or only if
                        // different?
                        if
                        (
                           !cptMap.insert(index, newIndex)
                         && cptMap[index] != newIndex
                        )
                        {
                            FatalErrorInFunction<< "Duplicate insertion"
                                << "From oldProc:" << oldProci
                                << " on map:" << mapi
                                << " at index:" << i
                                << " have construct slot:" << index
                                << " new index:" << newIndex
                                << " but already have entry:" << cptMap[index]
                                << " on for that slot"
                                << exit(FatalError);
                        }
                    }
                }
            }
        }


        // Remote ranks that are still remote
        //  - store new index for mapping stencils
        //  - append to construction

        // Either loop over all old ranks and filter out ones already handled
        // or loop over all new ranks and avoid myNewRank

        forAll(newToOldRanks, newProci)
        {
            if (newProci != myNewRank)
            {
                const auto& oldProcs = newToOldRanks[newProci];

                label allSize = 0;
                forAll(maps, mapi)
                {
                    for (const label oldProci : oldProcs)
                    {
                        allSize += maps[mapi].constructMap()[oldProci].size();
                    }
                }

                labelList& myConstruct = constructMap_[newProci];
                myConstruct.resize_nocopy(allSize);

                allSize = 0;
                forAll(maps, mapi)
                {
                    for (const label oldProci : oldProcs)
                    {
                        const auto& map = maps[mapi].constructMap()[oldProci];
                        // Construct map starts after the local data
                        const label nMapLocal = startOfRemote[mapi][oldProci];
                        SubList<label> slice(myConstruct, map.size(), allSize);

                        if (constructHasFlip_)
                        {
                            forAll(map, i)
                            {
                                if (map[i] < 0)
                                {
                                    slice[i] = map[i]+nMapLocal-constructi;
                                }
                                else
                                {
                                    slice[i] = map[i]-nMapLocal+constructi;
                                }
                            }

                            auto& cptMap = compactMaps[mapi];
                            forAll(map, i)
                            {
                                cptMap.insert(mag(map[i])-1,mag(slice[i])-1);
                            }
                        }
                        else
                        {
                            forAll(map, i)
                            {
                                slice[i] = map[i]-nMapLocal+constructi;
                                compactMaps[mapi].insert(map[i], slice[i]);
                            }
                        }
                        allSize += map.size();
                        constructi += map.size();
                    }
                }
            }
        }
    }


    // Sub (=send) map
    // ~~~~~~~~~~~~~~~
    // - all localRanks:
    //      - get appended in map order
    // - all previously remote ranks:
    //      - not needed. Stay empty
    // - all still remote ranks:
    //      - convert to new local index

    // Append local (= myRank) data
    {
        label allSize = 0;
        forAll(maps, mapi)
        {
            const label localRank = localRanks[mapi];
            allSize += maps[mapi].subMap()[localRank].size();
        }

        labelList& mySub = subMap_[myNewRank];
        mySub.resize_nocopy(allSize);
        allSize = 0;
        forAll(maps, mapi)
        {
            const label localRank = localRanks[mapi];
            const auto& map = maps[mapi].subMap()[localRank];
            SubList<label> slice(mySub, map.size(), allSize);

            if (subHasFlip_)
            {
                forAll(slice, i)
                {
                    if (map[i] < 0)
                    {
                        slice[i] = map[i]-startOfLocal[mapi];
                    }
                    else
                    {
                        slice[i] = map[i]+startOfLocal[mapi];
                    }
                }
            }
            else
            {
                forAll(slice, i)
                {
                    slice[i] = map[i]+startOfLocal[mapi];
                }
            }
            allSize += map.size();
        }
    }
    // Filter remote sub data
    forAll(newToOldRanks, newProci)
    {
        if (newProci != myNewRank)
        {
            const auto& oldProcs = newToOldRanks[newProci];

            label allSize = 0;
            forAll(maps, mapi)
            {
                for (const label oldProci : oldProcs)
                {
                    allSize += maps[mapi].subMap()[oldProci].size();
                }
            }

            labelList& mySub = subMap_[newProci];
            mySub.resize_nocopy(allSize);

            allSize = 0;
            for (const label oldProci : oldProcs)
            {
                forAll(maps, mapi)
                {
                    const auto& map = maps[mapi].subMap()[oldProci];
                    SubList<label> slice(mySub, map.size(), allSize);
                    if (subHasFlip_)
                    {
                        forAll(map, i)
                        {
                            if (map[i] < 0)
                            {
                                slice[i] = map[i]-startOfLocal[mapi];
                            }
                            else
                            {
                                slice[i] = map[i]+startOfLocal[mapi];
                            }
                        }
                    }
                    else
                    {
                        forAll(map, i)
                        {
                            slice[i] = map[i]+startOfLocal[mapi];
                        }
                    }
                    allSize += map.size();
                }
            }
        }
    }


    constructSize_ = constructi;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::mapDistributeBase::subMapSizes() const
{
    labelList sizes(subMap_.size());
    forAll(subMap_, i)
    {
        sizes[i] = subMap_[i].size();
    }
    return sizes;
}


Foam::labelList Foam::mapDistributeBase::constructMapSizes() const
{
    labelList sizes(constructMap_.size());
    forAll(constructMap_, i)
    {
        sizes[i] = constructMap_[i].size();
    }
    return sizes;
}


Foam::label Foam::mapDistributeBase::subMapTotalSize() const noexcept
{
    label total = 0;
    for (const auto& list : subMap_)
    {
        total += list.size();
    }
    return total;
}


Foam::label Foam::mapDistributeBase::constructMapTotalSize() const noexcept
{
    label total = 0;
    for (const auto& list : constructMap_)
    {
        total += list.size();
    }
    return total;
}


void Foam::mapDistributeBase::clear()
{
    constructSize_ = 0;
    subMap_.clear();
    constructMap_.clear();
    subHasFlip_ = false;
    constructHasFlip_ = false;
    // Leave comm_ intact
    schedulePtr_.reset(nullptr);
}


void Foam::mapDistributeBase::transfer(mapDistributeBase& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
        return;
    }

    constructSize_ = rhs.constructSize_;
    subMap_.transfer(rhs.subMap_);
    constructMap_.transfer(rhs.constructMap_);
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    comm_ = rhs.comm_;
    schedulePtr_.reset(nullptr);

    rhs.constructSize_ = 0;
    rhs.subHasFlip_ = false;
    rhs.constructHasFlip_ = false;
}


Foam::label Foam::mapDistributeBase::renumber
(
    const globalIndex& globalNumbering,
    const label comm,
    const List<Map<label>>& compactMap,
    const label globalI
)
{
    const label myRank = Pstream::myProcNo(comm);

    if (globalI == -1)
    {
        return globalI;
    }
    if (globalNumbering.isLocal(myRank, globalI))
    {
        return globalNumbering.toLocal(myRank, globalI);
    }
    else
    {
        label proci = globalNumbering.whichProcID(myRank, globalI);
        label index = globalNumbering.toLocal(proci, globalI);
        return compactMap[proci][index];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistributeBase::operator=(const mapDistributeBase& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    constructSize_ = rhs.constructSize_;
    subMap_ = rhs.subMap_;
    constructMap_ = rhs.constructMap_;
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    comm_ = rhs.comm_;
    schedulePtr_.reset(nullptr);
}


void Foam::mapDistributeBase::operator=(mapDistributeBase&& rhs)
{
    if (this != &rhs)
    {
        // Avoid self assignment
        transfer(rhs);
    }
}


// ************************************************************************* //
