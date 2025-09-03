/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "globalIndex.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// Cannot use non-blocking for non-contiguous data.
// template<class Type>
// inline Foam::UPstream::commsTypes getCommsType
// (
//     const UPstream::commsTypes preferred
// )
// {
//     return
//     (
//         (
//             !is_contiguous<Type>::value
//          && UPstream::commsTypes::nonBlocking == preferred
//         )
//       ? UPstream::commsTypes::scheduled
//       : preferred
//     );
// }


template<class Addr>
Foam::labelList
Foam::globalIndex::calcOffsets
(
    const IndirectListBase<label, Addr>& counts,
    const bool checkOverflow
)
{
    labelList values;

    const label len = counts.size();

    if (len)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = counts[i];

            values[i] = start;
            start += count;

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, values[i], count);
            }
        }
        values[len] = start;
    }

    return values;
}


template<class SubListType>
Foam::labelList
Foam::globalIndex::calcListOffsets
(
    const List<SubListType>& lists,
    const bool checkOverflow
)
{
    labelList values;

    const label len = lists.size();

    if (len)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            const label count = lists[i].size();

            values[i] = start;
            start += count;

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, values[i], count);
            }
        }
        values[len] = start;
    }

    return values;
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gatherValues
(
    const label comm,
    const ProcIDsContainer& procIDs,
    const Type& localValue,
    List<Type>& allValues,
    const int tag,
    const UPstream::commsTypes preferredCommsType
)
{
    // low-level: no parRun guard

    // Cannot use non-blocking for non-contiguous data.
    const UPstream::commsTypes commsType =
    (
        (
            !is_contiguous<Type>::value
         && UPstream::commsTypes::nonBlocking == preferredCommsType
        )
      ? UPstream::commsTypes::scheduled
      : preferredCommsType
    );

    const label startOfRequests = UPstream::nRequests();

    const int masterProci = procIDs.size() ? procIDs[0] : 0;

    if (UPstream::myProcNo(comm) == masterProci)
    {
        allValues.resize_nocopy(procIDs.size());
        allValues[0] = localValue;

        for (label i = 1; i < procIDs.size(); ++i)
        {
            if (is_contiguous<Type>::value)
            {
                UIPstream::read
                (
                    commsType,
                    procIDs[i],
                    reinterpret_cast<char*>(&allValues[i]),
                    sizeof(Type),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(allValues[i], procIDs[i], tag, comm);
            }
        }
    }
    else
    {
        allValues.clear();  // safety: zero-size on non-master

        if (is_contiguous<Type>::value)
        {
            UOPstream::write
            (
                commsType,
                masterProci,
                reinterpret_cast<const char*>(&localValue),
                sizeof(Type),
                tag,
                comm
            );
        }
        else
        {
            OPstream::send(localValue, commsType, masterProci, tag, comm);
        }
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag,
    const UPstream::commsTypes preferredCommsType
)
{
    // low-level: no parRun guard

    // Cannot use non-blocking for non-contiguous data.
    const UPstream::commsTypes commsType =
    (
        (
            !is_contiguous<Type>::value
         && UPstream::commsTypes::nonBlocking == preferredCommsType
        )
      ? UPstream::commsTypes::scheduled
      : preferredCommsType
    );

    const label startOfRequests = UPstream::nRequests();

    const int masterProci = procIDs.size() ? procIDs[0] : 0;

    if (UPstream::myProcNo(comm) == masterProci)
    {
        allFld.resize_nocopy(off.back());  // == totalSize()

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.

        SubList<Type>(allFld, off[1]-off[0], off[0]) =
           SubList<Type>(fld, off[1]-off[0]);

        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

            if (procSlot.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UIPstream::read
                (
                    commsType,
                    procIDs[i],
                    procSlot.data_bytes(),
                    procSlot.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(procSlot, procIDs[i], tag, comm);
            }
        }
    }
    else
    {
        if (fld.empty())
        {
            // Nothing to do
        }
        else if (is_contiguous<Type>::value)
        {
            UOPstream::write
            (
                commsType,
                masterProci,
                fld.cdata_bytes(),
                fld.size_bytes(),
                tag,
                comm
            );
        }
        else
        {
            OPstream::send(fld, commsType, masterProci, tag, comm);
        }
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


template<class ProcIDsContainer, class Type, class Addr>
void Foam::globalIndex::gather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const IndirectListBase<Type, Addr>& fld,
    List<Type>& allFld,
    const int tag,
    const UPstream::commsTypes preferredCommsType
)
{
    // low-level: no parRun guard

    if (is_contiguous<Type>::value)
    {
        // Flatten list (locally) so that we can benefit from using direct
        // read/write of contiguous data

        gather
        (
            off,
            comm,
            procIDs,
            List<Type>(fld),
            allFld,
            tag,
            preferredCommsType
        );
        return;
    }

    // Cannot use non-blocking for non-contiguous data.
    const UPstream::commsTypes commsType =
    (
        (
            !is_contiguous<Type>::value
         && UPstream::commsTypes::nonBlocking == preferredCommsType
        )
      ? UPstream::commsTypes::scheduled
      : preferredCommsType
    );

    const label startOfRequests = UPstream::nRequests();

    const int masterProci = procIDs.size() ? procIDs[0] : 0;

    if (UPstream::myProcNo(comm) == masterProci)
    {
        allFld.resize_nocopy(off.back());  // == totalSize()

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied

        SubList<Type> localSlot(allFld, off[1]-off[0], off[0]);
        if (!localSlot.empty())
        {
            localSlot = fld;
        }

        // Already verified commsType != nonBlocking
        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

            if (procSlot.empty())
            {
                // Nothing to do
            }
            else
            {
                IPstream::recv(procSlot, procIDs[i], tag, comm);
            }
        }
    }
    else
    {
        if (fld.empty())
        {
            // Nothing to do
        }
        else
        {
            OPstream::send(fld, commsType, masterProci, tag, comm);
        }
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalIndex::gather
(
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    {
        globalIndex::gather
        (
            offsets_,  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType
        );
        if (!UPstream::master(comm))
        {
            allData.clear();  // safety: zero-size on non-master
        }
    }
}


template<class Type, class Addr>
void Foam::globalIndex::gather
(
    const IndirectListBase<Type, Addr>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    {
        globalIndex::gather
        (
            offsets_,  // needed on master only
            comm,
            UPstream::allProcs(comm),   // All communicator ranks
            sendData,
            allData,
            tag,
            commsType
        );
        if (!UPstream::master(comm))
        {
            allData.clear();  // safety: zero-size on non-master
        }
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::gather
(
    const UList<Type>& sendData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    OutputContainer allData;
    gather(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type, class Addr, class OutputContainer>
OutputContainer Foam::globalIndex::gather
(
    const IndirectListBase<Type, Addr>& sendData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    OutputContainer allData;
    gather(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type>
void Foam::globalIndex::gatherInplace
(
    List<Type>& fld,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        gather(fld, allData, tag, commsType, comm);

        if (UPstream::master(comm))
        {
            fld.transfer(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type, class OutputContainer>
void Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    OutputContainer& allData,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
) const
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    // MPI_Gatherv requires contiguous data, but a byte-wise transfer can
    // quickly exceed the 'int' limits used for MPI sizes/offsets.
    // Thus gather label/scalar components when possible to increase the
    // effective size limit.
    //
    // Note: cannot rely on pTraits (cmptType, nComponents) since this method
    // needs to compile (and work) even with things like strings etc.

    // Single char ad hoc "enum":
    // - b(yte):  gather bytes
    // - f(loat): gather scalars components
    // - i(nt):   gather label components
    // - 0:       gather with Pstream read/write etc.

    List<int> recvCounts;
    List<int> recvOffsets;

    char dataMode(0);
    int nCmpts(0);

    if (is_contiguous<Type>::value)
    {
        if (is_contiguous_scalar<Type>::value)
        {
            dataMode = 'f';
            nCmpts = static_cast<int>(sizeof(Type)/sizeof(scalar));
        }
        else if (is_contiguous_label<Type>::value)
        {
            dataMode = 'i';
            nCmpts = static_cast<int>(sizeof(Type)/sizeof(label));
        }
        else
        {
            dataMode = 'b';
            nCmpts = static_cast<int>(sizeof(Type));
        }

        // Offsets must fit into int
        if (UPstream::master(comm))
        {
            const globalIndex& globalAddr = *this;

            if (globalAddr.totalSize() > (INT_MAX/nCmpts))
            {
                // Offsets do not fit into int - revert to manual.
                dataMode = 0;
            }
            else
            {
                // Must be same as Pstream::nProcs(comm), at least on master!
                const label nproc = globalAddr.nProcs();

                allData.resize_nocopy(globalAddr.totalSize());

                recvCounts.resize(nproc);
                recvOffsets.resize(nproc+1);

                for (label proci = 0; proci < nproc; ++proci)
                {
                    recvCounts[proci] = globalAddr.localSize(proci)*nCmpts;
                    recvOffsets[proci] = globalAddr.localStart(proci)*nCmpts;
                }
                recvOffsets[nproc] = globalAddr.totalSize()*nCmpts;

                // Assign local data directly

                recvCounts[0] = 0;  // ie, ignore for MPI_Gatherv
                SubList<Type>(allData, globalAddr.range(0)) =
                    SubList<Type>(sendData, globalAddr.range(0));
            }
        }

        // Consistent information for everyone
        UPstream::broadcast(&dataMode, 1, comm);
    }

    // Dispatch
    switch (dataMode)
    {
        case 'b':   // Byte-wise
        {
            UPstream::gather
            (
                sendData.cdata_bytes(),
                sendData.size_bytes(),
                allData.data_bytes(),
                recvCounts,
                recvOffsets,
                comm
            );
            break;
        }
        case 'f':   // Float (scalar) components
        {
            typedef scalar cmptType;

            UPstream::gather
            (
                reinterpret_cast<const cmptType*>(sendData.cdata()),
                (sendData.size()*nCmpts),
                reinterpret_cast<cmptType*>(allData.data()),
                recvCounts,
                recvOffsets,
                comm
            );
            break;
        }
        case 'i':   // Int (label) components
        {
            typedef label cmptType;

            UPstream::gather
            (
                reinterpret_cast<const cmptType*>(sendData.cdata()),
                (sendData.size()*nCmpts),
                reinterpret_cast<cmptType*>(allData.data()),
                recvCounts,
                recvOffsets,
                comm
            );
            break;
        }
        default:    // Regular (manual) gathering
        {
            globalIndex::gather
            (
                offsets_,  // needed on master only
                comm,
                UPstream::allProcs(comm),   // All communicator ranks
                sendData,
                allData,
                tag,
                commsType
            );
            break;
        }
    }

    if (!UPstream::master(comm))
    {
        allData.clear();  // safety: zero-size on non-master
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
) const
{
    OutputContainer allData;
    mpiGather(sendData, allData, comm, commsType, tag);
    return allData;
}


template<class Type>
void Foam::globalIndex::mpiGatherInplace
(
    List<Type>& fld,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
) const
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        mpiGather(fld, allData, comm, commsType, tag);

        if (UPstream::master(comm))
        {
            fld.transfer(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type, class OutputContainer>
void Foam::globalIndex::mpiGatherOp
(
    const UList<Type>& sendData,
    OutputContainer& allData,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .mpiGather(sendData, allData, comm, commsType, tag);
    }
    else
    {
        // Serial: direct copy
        allData = sendData;
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::mpiGatherOp
(
    const UList<Type>& sendData,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
)
{
    OutputContainer allData;
    mpiGatherOp(sendData, allData, comm, commsType, tag);
    return allData;
}


template<class Type>
void Foam::globalIndex::mpiGatherInplaceOp
(
    List<Type>& fld,
    const label comm,

    const UPstream::commsTypes commsType,
    const int tag
)
{
    if (UPstream::parRun())
    {
        List<Type> allData;
        mpiGatherOp(fld, allData, comm, commsType, tag);

        if (UPstream::master(comm))
        {
            fld.transfer(allData);
        }
        else
        {
            fld.clear();  // zero-size on non-master
        }
    }
    // Serial: (no-op)
}


template<class Type>
void Foam::globalIndex::gatherOp
(
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .gather(sendData, allData, tag, commsType, comm);
    }
    else
    {
        // Serial: direct copy
        allData = sendData;
    }
}


template<class Type, class Addr>
void Foam::globalIndex::gatherOp
(
    const IndirectListBase<Type, Addr>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, sendData.size(), comm)
            .gather(sendData, allData, tag, commsType, comm);
    }
    else
    {
        // Serial: direct copy
        allData = List<Type>(sendData);
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::gatherOp
(
    const UList<Type>& sendData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
)
{
    OutputContainer allData;
    gatherOp(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type, class Addr, class OutputContainer>
OutputContainer Foam::globalIndex::gatherOp
(
    const IndirectListBase<Type, Addr>& sendData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
)
{
    OutputContainer allData;
    gatherOp(sendData, allData, tag, commsType, comm);
    return allData;
}


template<class Type>
void Foam::globalIndex::gatherInplaceOp
(
    List<Type>& fld,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
)
{
    if (UPstream::parRun())
    {
        // Gather sizes - only needed on master
        globalIndex(globalIndex::gatherOnly{}, fld.size(), comm)
            .gather(fld, tag, commsType, comm);
    }
    // Serial: (no-op)
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::scatter
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& allFld,
    UList<Type>& fld,
    const int tag,
    const UPstream::commsTypes preferredCommsType
)
{
    // low-level: no parRun guard

    // Cannot use non-blocking for non-contiguous data.
    const UPstream::commsTypes commsType =
    (
        (
            !is_contiguous<Type>::value
         && UPstream::commsTypes::nonBlocking == preferredCommsType
        )
      ? UPstream::commsTypes::scheduled
      : preferredCommsType
    );

    const label startOfRequests = UPstream::nRequests();

    const int masterProci = procIDs.size() ? procIDs[0] : 0;

    if (UPstream::myProcNo(comm) == masterProci)
    {
        for (label i = 1; i < procIDs.size(); ++i)
        {
            const SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

            if (procSlot.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UOPstream::write
                (
                    commsType,
                    procIDs[i],
                    procSlot.cdata_bytes(),
                    procSlot.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream::send(procSlot, commsType, procIDs[i], tag, comm);
            }
        }

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.

        SubList<Type>(fld, off[1]-off[0]) =
            SubList<Type>(allFld, off[1]-off[0], off[0]);
    }
    else
    {
        // Note: we are receiving into UList, so sizes MUST match or we
        // have a problem. Can therefore reasonably assume that a zero-sized
        // send matches a zero-sized receive, and we can skip that.

        if (fld.empty())
        {
            // Nothing to do
        }
        else if (is_contiguous<Type>::value)
        {
            UIPstream::read
            (
                commsType,
                masterProci,
                fld.data_bytes(),
                fld.size_bytes(),
                tag,
                comm
            );
        }
        else
        {
            IPstream::recv(fld, masterProci, tag, comm);
        }
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


template<class Type>
void Foam::globalIndex::scatter
(
    const UList<Type>& allData,
    UList<Type>& localData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        scatter
        (
            offsets_,  // needed on master only
            comm,
            UPstream::allProcs(comm),  // All communicator ranks
            allData,
            localData,
            tag,
            commsType
        );
    }
    else
    {
        // Serial: direct copy
        // - fails miserably if incorrectly dimensioned!
        localData.deepCopy(allData);
    }
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::scatter
(
    const UList<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm
) const
{
    if (UPstream::parRun())
    {
        // The globalIndex might be correct on master only,
        // so scatter local sizes to ensure consistency

        const label count
        (
            UPstream::listScatterValues<label>(this->localSizes(), comm)
        );

        OutputContainer localData(count);
        this->scatter(allData, localData, tag, commsType, comm);

        return localData;
    }
    else
    {
        // Serial: direct copy
        return OutputContainer(allData);
    }
}


template<class Type, class CombineOp>
void Foam::globalIndex::get
(
    List<Type>& allFld,
    const labelUList& globalIds,
    const CombineOp& cop,
    const label comm,
    const int tag
) const
{
    allFld.resize_nocopy(globalIds.size());

    if (globalIds.size())
    {
        // Sort according to processor
        labelList order;
        DynamicList<label> validBins(Pstream::nProcs());

        CompactListList<label> bins
        (
            bin(offsets(), globalIds, order, validBins)
        );

        // Send local indices to individual processors as local index
        PstreamBuffers sendBufs(comm, tag);

        for (const auto proci : validBins)
        {
            labelList localIDs(bins[proci]);

            for (label& val : localIDs)
            {
                val = toLocal(proci, val);
            }

            UOPstream os(proci, sendBufs);
            os << localIDs;
        }
        sendBufs.finishedSends();


        PstreamBuffers returnBufs(comm, tag);

        for (const int proci : sendBufs.allProcs())
        {
            if (sendBufs.recvDataCount(proci))
            {
                UIPstream is(proci, sendBufs);
                labelList localIDs(is);

                // Collect entries
                List<Type> fld(localIDs.size());
                cop(fld, localIDs);

                UOPstream os(proci, returnBufs);
                os << fld;
            }
        }
        returnBufs.finishedSends();

        // Slot back
        for (const auto proci : validBins)
        {
            label start = bins.offsets()[proci];
            const SubList<label> es
            (
                order,
                bins.offsets()[proci+1]-start,  // start
                start
            );
            UIPstream is(proci, returnBufs);
            List<Type> fld(is);

            UIndirectList<Type>(allFld, es) = fld;
        }
    }
}


// ************************************************************************* //
