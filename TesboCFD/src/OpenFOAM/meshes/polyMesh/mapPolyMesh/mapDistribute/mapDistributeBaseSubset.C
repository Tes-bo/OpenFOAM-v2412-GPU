/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
#include "ListOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Setup array of element masks to match maps sizes
//
// \param[out] masks Sized for each position in the maps
// \param maps  The element maps
static void blankElementMasks(List<bitSet>& masks, const labelListList& maps)
{
    // If base container not already sized
    if (masks.empty())
    {
        masks.resize(maps.size());
    }

    forAll(masks, proci)
    {
        masks[proci].reset();  // zero all bits
        masks[proci].resize(maps[proci].size());
    }
}


// Calculate the element mask correspondig to allowedElems in the maps
//
// \param allowedElems Permissible mapped elements (true/false)
// \param[out] masks   True/false for each position within the maps
// \param maps     The element maps
// \param hasFlip  Map has flip indexing
//
// \return the max index used.
static label calcElementMasks
(
    const bitSet& allowedElems,
    List<bitSet>& masks,   // [out] - often presized before calling
    const labelListList& maps,
    const bool hasFlip
)
{
    // Index after flipping
    const auto unflippedIndex =
    (
        hasFlip
      ? [](label idx) -> label { return mag(idx)-1; }
      : [](label idx) -> label { return idx; }
    );


    // If not already sized
    if (masks.empty())
    {
        masks.resize(maps.size());
    }

    label maxIndex = -1;

    forAll(masks, proci)
    {
        bitSet& mask = masks[proci];
        const labelList& map = maps[proci];

        mask.reset();  // zero all bits
        mask.resize(map.size());

        forAll(map, i)
        {
            // Element is used (or not)
            const label index = unflippedIndex(map[i]);

            if (allowedElems.test(index))
            {
                mask.set(i);
                maxIndex = max(maxIndex, index);
            }
        }
    }

    return maxIndex;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::mapDistributeBase::exchangeMasks
(
    const UList<bitSet>& sendMasks,
    UList<bitSet>& recvMasks,
    const int tag,
    const label comm
)
{
    // Require properly sized mask buffers.
    // The information *is* known from the maps, so always use that
    // instead having a needless all-to-all for the sizes.

    if (sendMasks.size() != recvMasks.size())
    {
        FatalErrorInFunction
            << "Mismatched mask sizes: "
            << sendMasks.size() << " != "
            << recvMasks.size() << nl
            << Foam::abort(FatalError);
    }

    const label myRank = UPstream::myProcNo(comm);

    if (UPstream::parRun())
    {
        #ifdef FULLDEBUG
        if (sendMasks.size() > UPstream::nProcs(comm))
        {
            FatalErrorInFunction
                << "Mask sizes (" << sendMasks.size()
                << ") are larger than number of procs:"
                << UPstream::nProcs(comm) << nl
                << Foam::abort(FatalError);
        }
        #endif

        const label startOfRequests = UPstream::nRequests();

        forAll(recvMasks, proci)
        {
            if (proci != myRank && recvMasks[proci].size())
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    recvMasks[proci].data_bytes(),
                    recvMasks[proci].size_bytes(),
                    tag,
                    comm
                );
            }
        }

        forAll(sendMasks, proci)
        {
            if (proci != myRank && sendMasks[proci].size())
            {
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendMasks[proci].cdata_bytes(),
                    sendMasks[proci].size_bytes(),
                    tag,
                    comm
                );
            }
        }

        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }

    // Receiving myself is just a copy
    recvMasks[myRank] = sendMasks[myRank];
}


void Foam::mapDistributeBase::unionCombineMasks
(
    UList<bitSet>& sendMasks,
    UList<bitSet>& recvMasks,
    const int tag,
    const label comm
)
{
    // Require properly sized mask buffers.
    // The information *is* known from the maps, so always use that
    // instead having a needless all-to-all for the sizes.

    if (sendMasks.size() != recvMasks.size())
    {
        FatalErrorInFunction
            << "Mismatched mask sizes: "
            << sendMasks.size() << " != "
            << recvMasks.size() << nl
            << Foam::abort(FatalError);
    }

    if (UPstream::parRun())
    {
        // Scratch buffers for union operations
        List<bitSet> scratch(recvMasks.size());

        // Size for receives
        forAll(scratch, proci)
        {
            scratch[proci].resize(recvMasks[proci].size());
        }

        // Exchange: from sendMasks -> scratch (intermediate receive)
        exchangeMasks(sendMasks, scratch, tag, comm);

        // Update recvMasks (as union)
        forAll(recvMasks, proci)
        {
            recvMasks[proci] &= scratch[proci];
        }

        // Size for sends
        forAll(scratch, proci)
        {
            scratch[proci].resize(sendMasks[proci].size());
        }

        // Exchange: from recvMasks -> scratch (intermediate send)
        exchangeMasks(recvMasks, scratch, tag, comm);

        // Final synchronization
        forAll(sendMasks, proci)
        {
            sendMasks[proci] &= scratch[proci];
        }
    }
    else
    {
        // Non-parallel: 'synchronize' myself
        const label myRank = UPstream::myProcNo(comm);

        recvMasks[myRank] &= sendMasks[myRank];
        sendMasks[myRank] = recvMasks[myRank];
    }

    // Done with parallel exchanges so can shrink the masks to
    // the min-size actually needed.

    for (auto& mask : sendMasks)
    {
        mask.resize_last();
    }

    for (auto& mask : recvMasks)
    {
        mask.resize_last();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::mapDistributeBase::renumberMap
(
    labelListList& mapElements,
    const labelUList& oldToNew,
    const bool hasFlip
)
{
    label maxIndex = -1;

    // Transcribe the map
    if (hasFlip)
    {
        for (labelList& map : mapElements)
        {
            for (label& val : map)
            {
                // Unflip indexed value
                const label index = oldToNew[mag(val)-1];

                if (index >= 0)   // Not certain this check is needed
                {
                    maxIndex = max(maxIndex, index);

                    // Retain flip information from original
                    val = (val < 0 ? (-index-1) : (index+1));
                }
            }
        }
    }
    else
    {
        for (labelList& map : mapElements)
        {
            for (label& val : map)
            {
                // Get indexed value (no flipping)

                const label index = oldToNew[val];

                if (index >= 0)   // Not certain this check is needed
                {
                    maxIndex = max(maxIndex, index);
                    val = index;
                }
            }
        }
    }

    return (maxIndex+1);
}


Foam::label Foam::mapDistributeBase::renumberMap
(
    labelList& map,
    const label localSize,
    const label offset,
    const Map<label>& cMap,
    const bool hasFlip
)
{
    label maxIndex = -1;

    // Transcribe the map
    if (hasFlip)
    {
        for (label& val : map)
        {
            // Unflip indexed value
            const label index = mag(val)-1;
            if (index < localSize)
            {
                // Local element
                if (val < 0)
                {
                    val -= offset;
                }
                else
                {
                    val += offset;
                }
            }
            else
            {
                // Remote element
                if (val < 0)
                {
                    val = -cMap[index]-1;
                }
                else
                {
                    val = cMap[index]+1;
                }
            }
            maxIndex = max(maxIndex, mag(val)-1);
        }
    }
    else
    {
        for (label& val : map)
        {
            // Get indexed value (no flipping)
            if (val < localSize)
            {
                val += offset;
            }
            else
            {
                val = cMap[val];
            }
            maxIndex = max(maxIndex, val);
        }
    }

    return (maxIndex+1);
}


Foam::label Foam::mapDistributeBase::renumberMap
(
    labelListList& mapElements,
    const label localSize,
    const label offset,
    const Map<label>& cMap,
    const bool hasFlip
)
{
    label maxIndex = -1;

    // Transcribe the map
    for (labelList& map : mapElements)
    {
        maxIndex = max
        (
            maxIndex,
            renumberMap
            (
                map,
                localSize,
                offset,
                cMap,
                hasFlip
            )
        );
    }

    return (maxIndex+1);
}


void Foam::mapDistributeBase::renumberVisitOrder
(
    const labelUList& origElements,
    labelList& oldToNew,
    labelListList& maps,
    const bool hasFlip
)
{
    // Both oldToNew and maps refer to compacted numbers in simple
    // ascending order, but we want to recover the original walk order.

    // CAUTION:
    // The following is ill-defined (ie, really bad idea) if the original
    // elements contained duplicates!

    // Inverse mapping:
    //   Original id -> compact id -> walked id

    labelList compactToWalkOrder(origElements.size(), -1);

    forAll(origElements, walkIndex)
    {
        const label origIndex = origElements[walkIndex];
        const label compactIndex = oldToNew[origIndex];

        if (compactIndex >= origElements.size())
        {
            FatalErrorInFunction
                << "Compact index: " << compactIndex
                << " is not represented in the original ("
                << origElements.size()
                << ") elements - indicates an addressing problem" << nl
                << Foam::abort(FatalError);
        }
        else if (compactIndex >= 0)
        {
            compactToWalkOrder[compactIndex] = walkIndex;
            oldToNew[origIndex] = walkIndex;
        }
    }

    renumberMap(maps, compactToWalkOrder, hasFlip);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapDistributeBase::calcCompactDataRequirements
(
    const bitSet& allowedLocalElems,
    const bitSet& allowedRemoteElems,
    List<bitSet>& sendMasks,     // [out]
    List<bitSet>& recvMasks,     // [out]
    const int tag
)
{
    sendMasks.resize_nocopy(UPstream::nProcs(comm_));
    recvMasks.resize_nocopy(UPstream::nProcs(comm_));

    // Determine local elements sent to which procs
    calcElementMasks
    (
        allowedLocalElems,
        sendMasks,
        subMap_,
        subHasFlip_
    );

    // Determine remote elements received from which procs
    calcElementMasks
    (
        allowedRemoteElems,
        recvMasks,
        constructMap_,
        constructHasFlip_
    );

    // Synchronize - combine as '&' union
    unionCombineMasks(sendMasks, recvMasks, tag, comm_);
}


void Foam::mapDistributeBase::calcCompactLocalDataRequirements
(
    const bitSet& allowedLocalElems,
    List<bitSet>& sendMasks,        // [out]
    List<bitSet>& recvMasks,        // [out]
    const int tag
)
{
    sendMasks.resize_nocopy(UPstream::nProcs(comm_));
    recvMasks.resize_nocopy(UPstream::nProcs(comm_));

    // Determine local elements sent to which procs
    calcElementMasks
    (
        allowedLocalElems,
        sendMasks,
        subMap_,
        subHasFlip_
    );

    blankElementMasks(recvMasks, constructMap_);

    // Exchange: from sendMasks -> recvMasks
    exchangeMasks(sendMasks, recvMasks, tag, comm_);
}


void Foam::mapDistributeBase::calcCompactRemoteDataRequirements
(
    const bitSet& allowedRemoteElems,
    List<bitSet>& sendMasks,        // [out]
    List<bitSet>& recvMasks,        // [out]
    const int tag
)
{
    sendMasks.resize_nocopy(UPstream::nProcs(comm_));
    recvMasks.resize_nocopy(UPstream::nProcs(comm_));

    // Determine remote elements received from which procs
    calcElementMasks
    (
        allowedRemoteElems,
        recvMasks,
        constructMap_,
        constructHasFlip_
    );

    blankElementMasks(sendMasks, subMap_);

    // Exchange: from recvMasks -> sendMasks
    exchangeMasks(recvMasks, sendMasks, tag, comm_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapDistributeBase::compactData
(
    const UList<bitSet>& sendMasks,
    const UList<bitSet>& recvMasks,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize  // (known) max sizing for subMap
)
{
    // Linear address (subMap) after any flipping
    const auto unflippedSendIndex =
    (
        subHasFlip_
      ? [](label idx) -> label { return mag(idx)-1; }
      : [](label idx) -> label { return idx; }
    );

    // Linear address (constructMap) after any flipping
    const auto unflippedRecvIndex =
    (
        constructHasFlip_
      ? [](label idx) -> label { return mag(idx)-1; }
      : [](label idx) -> label { return idx; }
    );


    // Compact renumbering enabled if oldToNew maps are notNull

    bitSet indexUsed;

    // The subMap old-to-new mapping
    if (notNull(oldToNewSub))
    {
        label subMapSize(localSize);
        if (subMapSize < 0)
        {
            subMapSize = getMappedSize(subMap_, subHasFlip_);
        }

        oldToNewSub.resize_nocopy(subMapSize);
        oldToNewSub = -1;

        indexUsed.reset();  // zero all bits
        indexUsed.resize(subMapSize);

        forAll(sendMasks, proci)
        {
            const bitSet& mask = sendMasks[proci];
            const auto& map = subMap_[proci];

            for (const label i : mask)
            {
                const label index = unflippedSendIndex(map[i]);

                indexUsed.set(index);
            }
        }

        label nCompact = 0;
        for (const label i : indexUsed)
        {
            oldToNewSub[i] = nCompact++;
        }
    }


    // The constructMap old-to-new mapping
    if (notNull(oldToNewConstruct))
    {
        oldToNewConstruct.resize_nocopy(constructSize_);
        oldToNewConstruct = -1;

        indexUsed.reset();  // zero all bits
        indexUsed.resize(constructSize_);

        forAll(recvMasks, proci)
        {
            const bitSet& mask = recvMasks[proci];
            const auto& map = constructMap_[proci];

            for (const label i : mask)
            {
                const label index = unflippedRecvIndex(map[i]);

                indexUsed.set(index);
            }
        }

        label nCompact = 0;
        for (const label i : indexUsed)
        {
            oldToNewConstruct[i] = nCompact++;
        }
    }


    // Compact out subMap entries referring to unused elements
    forAll(sendMasks, proci)
    {
        const bitSet& mask = sendMasks[proci];
        labelList& map = subMap_[proci];

        label nCompact = 0;

        for (const label i : mask)
        {
            // const label index = unflippedSendIndex(map[i]);
            // maxLocalIndex = max(maxLocalIndex, index);

            map[nCompact++] = map[i];
        }

        map.resize(nCompact);
    }


    // Compact out constructMap entries referring to unused elements

    label maxRemoteIndex = -1;

    forAll(recvMasks, proci)
    {
        const bitSet& mask = recvMasks[proci];
        labelList& map = constructMap_[proci];

        label nCompact = 0;

        for (const label i : mask)
        {
            const label index = unflippedRecvIndex(map[i]);
            maxRemoteIndex = max(maxRemoteIndex, index);

            map[nCompact++] = map[i];
        }

        map.resize(nCompact);
    }

    constructSize_ = maxRemoteIndex+1;


    // Do compact renumbering...

    if (notNull(oldToNewSub))
    {
        renumberMap(subMap_, oldToNewSub, subHasFlip_);
    }

    if (notNull(oldToNewConstruct))
    {
        constructSize_ =
            renumberMap(constructMap_, oldToNewConstruct, constructHasFlip_);
    }

    // Clear the schedule (note:not necessary if nothing changed)
    schedulePtr_.reset(nullptr);
}


void Foam::mapDistributeBase::compactData
(
    const labelUList& localElements,
    const labelUList& remoteElements,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize,
    const int tag
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactDataRequirements
    (
        bitSet(localElements),
        bitSet(remoteElements),
        sendMasks,
        recvMasks,
        tag
    );

    // Perform compaction and renumbering
    compactData
    (
        sendMasks,
        recvMasks,
        oldToNewSub,
        oldToNewConstruct,
        localSize
    );

    // Renumber according to visit order
    renumberVisitOrder
    (
        localElements,
        oldToNewSub,
        subMap_,
        subHasFlip_
    );

    // Renumber according to visit order
    renumberVisitOrder
    (
        remoteElements,
        oldToNewConstruct,
        constructMap_,
        constructHasFlip_
    );
}


void Foam::mapDistributeBase::compactLocalData
(
    const labelUList& localElements,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize,
    const int tag
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactLocalDataRequirements
    (
        // Retain items required on the local side
        bitSet(localElements),
        sendMasks,
        recvMasks,
        tag
    );

    // Perform compaction and renumbering
    compactData
    (
        sendMasks,
        recvMasks,
        oldToNewSub,
        oldToNewConstruct,
        localSize
    );

    // Renumber according to visit order
    renumberVisitOrder
    (
        localElements,
        oldToNewSub,
        subMap_,
        subHasFlip_
    );
}


void Foam::mapDistributeBase::compactRemoteData
(
    const labelUList& remoteElements,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize,
    const int tag
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactRemoteDataRequirements
    (
        // Retain items required on the remote side
        bitSet(remoteElements),
        sendMasks,
        recvMasks,
        tag
    );

    // Perform compaction and renumbering
    compactData
    (
        sendMasks,
        recvMasks,
        oldToNewSub,
        oldToNewConstruct,
        localSize
    );

    // Renumber according to visit order
    renumberVisitOrder
    (
        remoteElements,
        oldToNewConstruct,
        constructMap_,
        constructHasFlip_
    );
}


void Foam::mapDistributeBase::compactDataImpl
(
    const UList<bitSet>& sendMasks,
    const UList<bitSet>& recvMasks,
    const bool doRenumber
)
{
    if (doRenumber)
    {
        labelList oldToNewSub;
        labelList oldToNewConstruct;

        compactData
        (
            sendMasks,
            recvMasks,
            oldToNewSub,
            oldToNewConstruct,
            -1  // localSize: automatic
        );
    }
    else
    {
        // Call with placeholder values
        compactData
        (
            sendMasks,
            recvMasks,
            const_cast<labelList&>(labelList::null()),  // disabled
            const_cast<labelList&>(labelList::null()),  // disabled
            -1  // localSize: automatic
        );
    }
}


void Foam::mapDistributeBase::compactLocalData
(
    const bitSet& allowedLocalElems,
    const int tag,
    const bool doRenumber
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactLocalDataRequirements
    (
        allowedLocalElems,
        sendMasks,
        recvMasks,
        tag
    );

    compactDataImpl(sendMasks, recvMasks, doRenumber);
}


void Foam::mapDistributeBase::compactRemoteData
(
    const bitSet& allowedRemoteElems,
    const int tag,
    const bool doRenumber
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactRemoteDataRequirements
    (
        allowedRemoteElems,
        sendMasks,
        recvMasks,
        tag
    );

    compactDataImpl(sendMasks, recvMasks, doRenumber);
}


void Foam::mapDistributeBase::compactLocalData
(
    const bitSet& allowedLocalElems,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize,
    const int tag
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactLocalDataRequirements
    (
        allowedLocalElems,
        sendMasks,
        recvMasks,
        tag
    );

    compactData
    (
        sendMasks,
        recvMasks,
        oldToNewSub,
        oldToNewConstruct,
        localSize
    );
}


void Foam::mapDistributeBase::compactRemoteData
(
    const bitSet& allowedRemoteElems,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const label localSize,
    const int tag
)
{
    List<bitSet> sendMasks;
    List<bitSet> recvMasks;

    calcCompactRemoteDataRequirements
    (
        allowedRemoteElems,
        sendMasks,
        recvMasks,
        tag
    );

    compactData
    (
        sendMasks,
        recvMasks,
        oldToNewSub,
        oldToNewConstruct,
        localSize
    );
}


// * * * * * * * * * * * * * * * Housekeeping  * * * * * * * * * * * * * * * //

void Foam::mapDistributeBase::compact
(
    const boolList& remoteElemUsed,
    const int tag
)
{
    // Forward to bitSet version
    compactRemoteData(bitSet(remoteElemUsed), tag);
}


void Foam::mapDistributeBase::compact
(
    const boolList& remoteElemUsed,
    const label localSize,
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const int tag
)
{
    // Forward to bitSet version
    compactRemoteData
    (
        bitSet(remoteElemUsed),
        oldToNewSub,
        oldToNewConstruct,
        localSize,
        tag
    );
}


// ************************************************************************* //
