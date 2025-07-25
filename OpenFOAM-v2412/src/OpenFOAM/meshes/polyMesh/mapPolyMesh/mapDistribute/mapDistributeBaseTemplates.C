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

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp, class NegateOp>
void Foam::mapDistributeBase::flipAndCombine
(
    List<T>& lhs,
    const UList<T>& rhs,

    const labelUList& map,
    const bool hasFlip,
    const CombineOp& cop,
    const NegateOp& negOp
)
{
    const label len = map.size();

    if (hasFlip)
    {
        for (label i = 0; i < len; ++i)
        {
            const label index = map[i];

            if (index > 0)
            {
                cop(lhs[index-1], rhs[i]);
            }
            else if (index < 0)
            {
                cop(lhs[-index-1], negOp(rhs[i]));
            }
            else
            {
                FatalErrorInFunction
                    << "Illegal flip index '0' at " << i << '/' << map.size()
                    << " for list:" << rhs.size() << nl
                    << exit(FatalError);
            }
        }
    }
    else
    {
        for (label i = 0; i < len; ++i)
        {
            cop(lhs[map[i]], rhs[i]);
        }
    }
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::accessAndFlip
(
    List<T>& output,
    const UList<T>& values,
    const labelUList& map,
    const bool hasFlip,
    const NegateOp& negOp
)
{
    const label len = map.size();

    // FULLDEBUG: if (output.size() < len)  FatalError ...;

    if (hasFlip)
    {
        for (label i = 0; i < len; ++i)
        {
            const label index = map[i];

            if (index > 0)
            {
                output[i] = values[index-1];
            }
            else if (index < 0)
            {
                output[i] = negOp(values[-index-1]);
            }
            else
            {
                FatalErrorInFunction
                    << "Illegal flip index '0' at " << i << '/' << map.size()
                    << " for list:" << values.size() << nl
                    << exit(FatalError);
            }
        }
    }
    else
    {
        // Like indirect list
        for (label i = 0; i < len; ++i)
        {
            output[i] = values[map[i]];
        }
    }
}


template<class T, class NegateOp>
Foam::List<T> Foam::mapDistributeBase::accessAndFlip
(
    const UList<T>& values,
    const labelUList& map,
    const bool hasFlip,
    const NegateOp& negOp
)
{
    List<T> output(map.size());
    accessAndFlip(output, values, map, hasFlip, negOp);
    return output;
}


template<class T, class negateOp>
void Foam::mapDistributeBase::send
(
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    const UList<T>& field,

    labelRange& sendRequests,
    PtrList<List<T>>& sendFields,

    labelRange& recvRequests,
    PtrList<List<T>>& recvFields,

    const negateOp& negOp,
    const int tag,
    const label comm
)
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Only contiguous is currently supported"
            << exit(FatalError);
    }

    const auto myRank = UPstream::myProcNo(comm);
    const auto nProcs = UPstream::nProcs(comm);


    // Set up receives from neighbours
    recvRequests.start() = UPstream::nRequests();
    recvRequests.size() = 0;

    recvFields.resize(nProcs);

    for (const int proci : UPstream::allProcs(comm))
    {
        const labelList& map = constructMap[proci];

        if (proci == myRank)
        {
            // No communication for myself - but recvFields may be used
        }
        else if (map.empty())
        {
            // No receive necessary
            (void) recvFields.release(proci);
        }
        else
        {
            List<T>& subField = recvFields.try_emplace(proci);
            subField.resize_nocopy(map.size());

            UIPstream::read
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                subField.data_bytes(),
                subField.size_bytes(),
                tag,
                comm
            );
        }
    }

    // Finished setting up the receives
    recvRequests.size() = (UPstream::nRequests() - recvRequests.start());


    // Set up sends to neighbours
    sendRequests.start() = UPstream::nRequests();
    sendRequests.size() = 0;

    sendFields.resize(nProcs);

    for (const int proci : UPstream::allProcs(comm))
    {
        const labelList& map = subMap[proci];

        if (proci == myRank)
        {
            // No communication - sendFields not needed
            (void) sendFields.release(proci);
        }
        else if (map.empty())
        {
            // No send necessary
            (void) sendFields.release(proci);
        }
        else
        {
            List<T>& subField = sendFields.try_emplace(proci);
            subField.resize_nocopy(map.size());

            accessAndFlip(subField, field, map, subHasFlip, negOp);

            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                subField.cdata_bytes(),
                subField.size_bytes(),
                tag,
                comm
            );
        }
    }

    // Finished setting up the sends
    sendRequests.size() = (UPstream::nRequests() - sendRequests.start());


    // Set up 'send' to myself - copy directly into recvFields
    {
        const labelList& map = subMap[myRank];

        if (map.empty())
        {
            // Nothing to send/recv
            (void) recvFields.release(myRank);
        }
        else
        {
            List<T>& subField = recvFields.try_emplace(myRank);
            subField.resize_nocopy(map.size());

            accessAndFlip(subField, field, map, subHasFlip, negOp);
        }
    }
}


template<class T>
void Foam::mapDistributeBase::send
(
    const UList<T>& field,
    labelRange& sendRequests,
    PtrList<List<T>>& sendFields,
    labelRange& recvRequests,
    PtrList<List<T>>& recvFields,
    const int tag
) const
{
    send
    (
        subMap_,
        subHasFlip_,
        constructMap_,
        constructHasFlip_,
        field,
        sendRequests, sendFields,
        recvRequests, recvFields,
        flipOp(),
        tag,
        comm_
    );
}


template<class T, class CombineOp, class negateOp>
void Foam::mapDistributeBase::receive
(
    const label constructSize,
    const labelListList& constructMap,
    const bool constructHasFlip,
    const labelRange& requests,
    const UPtrList<List<T>>& recvFields,
    List<T>& field,
    const CombineOp& cop,
    const negateOp& negOp,
    const int tag,
    const label comm
)
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Only contiguous is currently supported"
            << exit(FatalError);
    }

    const auto myRank = UPstream::myProcNo(comm);
    const auto nProcs = UPstream::nProcs(comm);


    // Receiving from which procs - according to map information

    DynamicList<int> recvProcs(nProcs);
    for (const int proci : UPstream::allProcs(comm))
    {
        const labelList& map = constructMap[proci];

        if (proci != myRank && map.size())
        {
            recvProcs.push_back(proci);

            const auto* subFieldPtr = recvFields.get(proci);
            if (subFieldPtr)
            {
                checkReceivedSize(proci, map.size(), subFieldPtr->size());
            }
            else
            {
                FatalErrorInFunction
                    << "From processor " << proci
                    << " : unallocated receive field."
                    << " Expected size " << map.size()
                    << " on comm " << comm
                    << " with procs " << UPstream::nProcs(comm) << nl
                    << exit(FatalError);
            }
        }
    }


    // Combining bits - can reuse field storage
    field.resize_nocopy(constructSize);


    // Received sub field from myself : recvFields[myRank]
    if (recvFields.test(myRank))
    {
        const labelList& map = constructMap[myRank];
        const List<T>& subField = recvFields[myRank];

        // Unlikely to need a size check
        // checkReceivedSize(myRank, map.size(), subField.size());

        flipAndCombine
        (
            field,
            subField,
            map,
            constructHasFlip,
            cop,
            negOp
        );
    }


    // NB: do NOT use polling and dispatch here.
    // There is no certainty if the listed requests are still pending or
    // have already been waited on before calling this method.

    // Wait for (receive) requests, but the range may also include
    // other requests depending on what the caller provided

    UPstream::waitRequests(requests.start(), requests.size());


    // Process received fields

    {
        for (const int proci : recvProcs)
        {
            const labelList& map = constructMap[proci];
            const List<T>& subField = recvFields[proci];

            // Already checked the sizes previously
            // checkReceivedSize(proci, map.size(), subField.size());

            flipAndCombine
            (
                field,
                subField,
                map,
                constructHasFlip,
                cop,
                negOp
            );
        }
    }
}


template<class T>
void Foam::mapDistributeBase::receive
(
    const labelRange& requests,
    const UPtrList<List<T>>& recvFields,
    List<T>& field,
    const int tag
) const
{
    receive
    (
        constructSize_,
        constructMap_,
        constructHasFlip_,
        requests,
        recvFields,
        field,
        eqOp<T>(),
        flipOp(),
        tag,
        comm_
    );
}


template<class T, class CombineOp, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const T& nullValue,
    const CombineOp& cop,
    const NegateOp& negOp,
    const int tag,
    const label comm
)
{
    const auto myRank = UPstream::myProcNo(comm);
    const auto nProcs = UPstream::nProcs(comm);

    if (!UPstream::parRun())
    {
        // Do only me to me.

        List<T> subField
        (
            accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
        );

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[myRank];

        // Combining bits - can now reuse field storage
        field.resize_nocopy(constructSize);
        field = nullValue;

        flipAndCombine
        (
            field,
            subField,
            map,
            constructHasFlip,
            cop,
            negOp
        );

        return;
    }

    if (commsType == UPstream::commsTypes::buffered)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (const int proci : UPstream::allProcs(comm))
        {
            const labelList& map = subMap[proci];

            if (proci != myRank && map.size())
            {
                List<T> subField
                (
                    accessAndFlip(field, map, subHasFlip, negOp)
                );

                // buffered send
                OPstream os(commsType, proci, 0, tag, comm);
                os  << subField;
            }
        }

        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            // Combining bits - can now reuse field storage
            field.resize_nocopy(constructSize);
            field = nullValue;

            flipAndCombine
            (
                field,
                subField,
                map,
                constructHasFlip,
                cop,
                negOp
            );
        }

        // Receive and process sub-field from neighbours
        for (const int proci : UPstream::allProcs(comm))
        {
            const labelList& map = constructMap[proci];

            if (proci != myRank && map.size())
            {
                List<T> subField;
                IPstream::recv(subField, proci, tag, comm);

                checkReceivedSize(proci, map.size(), subField.size());

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    cop,
                    negOp
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField;
        newField.resize_nocopy(constructSize);
        newField = nullValue;

        // First handle self
        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            flipAndCombine
            (
                newField,
                subField,
                map,
                constructHasFlip,
                cop,
                negOp
            );
        }


        // Schedule will already have pruned 0-sized comms
        for (const labelPair& twoProcs : schedule)
        {
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            if (twoProcs.first() == myRank)
            {
                // I am send first, receive next
                const label nbrProc = twoProcs.second();

                {
                    const labelList& map = subMap[nbrProc];

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    OPstream::send(subField, nbrProc, tag, comm);
                }
                {
                    List<T> subField;
                    IPstream::recv(subField, nbrProc, tag, comm);

                    const labelList& map = constructMap[nbrProc];

                    checkReceivedSize(nbrProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        newField,
                        subField,
                        map,
                        constructHasFlip,
                        cop,
                        negOp
                    );
                }
            }
            else
            {
                // I am receive first, send next
                const label nbrProc = twoProcs.first();

                {
                    List<T> subField;
                    IPstream::recv(subField, nbrProc, tag, comm);

                    const labelList& map = constructMap[nbrProc];

                    checkReceivedSize(nbrProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        newField,
                        subField,
                        map,
                        constructHasFlip,
                        cop,
                        negOp
                    );
                }
                {
                    const labelList& map = subMap[nbrProc];

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    OPstream::send(subField, nbrProc, tag, comm);
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        const label startOfRequests = UPstream::nRequests();

        if (!is_contiguous<T>::value)
        {
            PstreamBuffers pBufs(comm, tag);

            // Stream data into buffer
            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = subMap[proci];

                if (proci != myRank && map.size())
                {
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    UOPstream os(proci, pBufs);
                    os  << subField;
                }
            }

            // Initiate receiving - do yet not block
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                List<T> subField
                (
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
                );

                // Combining bits - can now reuse field storage
                field.resize_nocopy(constructSize);
                field = nullValue;

                // Receive sub field from myself
                const labelList& map = constructMap[myRank];

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    cop,
                    negOp
                );
            }

            // Wait for receive requests (and the send requests too)
            UPstream::waitRequests(startOfRequests);

            // Receive and process neighbour fields
            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = constructMap[proci];

                if (proci != myRank && map.size())
                {
                    UIPstream is(proci, pBufs);
                    List<T> subField(is);

                    checkReceivedSize(proci, map.size(), subField.size());

                    flipAndCombine
                    (
                        field,
                        subField,
                        map,
                        constructHasFlip,
                        cop,
                        negOp
                    );
                }
            }
        }
        else
        {
            // Set up receives from neighbours

            List<List<T>> recvFields(nProcs);
            DynamicList<int> recvProcs(nProcs);

            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = constructMap[proci];

                if (proci != myRank && map.size())
                {
                    recvProcs.push_back(proci);
                    List<T>& subField = recvFields[proci];
                    subField.resize_nocopy(map.size());

                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        subField.data_bytes(),
                        subField.size_bytes(),
                        tag,
                        comm
                    );
                }
            }


            // Set up sends to neighbours

            List<List<T>> sendFields(nProcs);

            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = subMap[proci];

                if (proci != myRank && map.size())
                {
                    List<T>& subField = sendFields[proci];
                    subField.resize_nocopy(map.size());

                    accessAndFlip(subField, field, map, subHasFlip, negOp);

                    UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        subField.cdata_bytes(),
                        subField.size_bytes(),
                        tag,
                        comm
                    );
                }
            }

            // Set up 'send' to myself - copy directly into recvFields
            {
                const labelList& map = subMap[myRank];
                List<T>& subField = recvFields[myRank];
                subField.resize_nocopy(map.size());

                accessAndFlip(subField, field, map, subHasFlip, negOp);
            }


            // Combining bits - can now reuse field storage
            field.resize_nocopy(constructSize);
            field = nullValue;

            // Receive sub field from myself : recvFields[myRank]
            {
                const labelList& map = constructMap[myRank];
                const List<T>& subField = recvFields[myRank];

                // Probably don't need a size check
                // checkReceivedSize(myRank, map.size(), subField.size());

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    cop,
                    negOp
                );
            }


            // Poll for completed receive requests and dispatch
            DynamicList<int> indices(recvProcs.size());
            while
            (
                UPstream::waitSomeRequests
                (
                    startOfRequests,
                    recvProcs.size(),
                   &indices
                )
            )
            {
                for (const int idx : indices)
                {
                    const int proci = recvProcs[idx];
                    const labelList& map = constructMap[proci];
                    const List<T>& subField = recvFields[proci];

                    // No size check - was dimensioned above
                    // checkReceivedSize(proci, map.size(), subField.size());

                    flipAndCombine
                    (
                        field,
                        subField,
                        map,
                        constructHasFlip,
                        cop,
                        negOp
                    );
                }
            }

            // Wait for any remaining requests
            UPstream::waitRequests(startOfRequests);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const NegateOp& negOp,
    const int tag,
    const label comm
)
{
    const auto myRank = UPstream::myProcNo(comm);
    const auto nProcs = UPstream::nProcs(comm);

    if (!UPstream::parRun())
    {
        // Do only me to me.

        List<T> subField
        (
            accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
        );

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[myRank];

        // Combining bits - can now reuse field storage
        field.resize_nocopy(constructSize);

        flipAndCombine
        (
            field,
            subField,
            map,
            constructHasFlip,
            eqOp<T>(),
            negOp
        );

        return;
    }

    if (commsType == UPstream::commsTypes::buffered)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (const int proci : UPstream::allProcs(comm))
        {
            const labelList& map = subMap[proci];

            if (proci != myRank && map.size())
            {
                List<T> subField
                (
                    accessAndFlip(field, map, subHasFlip, negOp)
                );

                // buffered send
                OPstream os(commsType, proci, 0, tag, comm);
                os  << subField;
            }
        }

        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            // Combining bits - can now reuse field storage
            field.resize_nocopy(constructSize);

            flipAndCombine
            (
                field,
                subField,
                map,
                constructHasFlip,
                eqOp<T>(),
                negOp
            );
        }

        // Receive and process sub-field from neighbours
        for (const int proci : UPstream::allProcs(comm))
        {
            const labelList& map = constructMap[proci];

            if (proci != myRank && map.size())
            {
                List<T> subField;
                IPstream::recv(subField, proci, tag, comm);

                checkReceivedSize(proci, map.size(), subField.size());

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    eqOp<T>(),
                    negOp
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField;
        newField.resize_nocopy(constructSize);

        // First handle self
        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            flipAndCombine
            (
                newField,
                subField,
                map,
                constructHasFlip,
                eqOp<T>(),
                negOp
            );
        }

        // Schedule will already have pruned 0-sized comms
        for (const labelPair& twoProcs : schedule)
        {
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            if (twoProcs.first() == myRank)
            {
                // I am send first, receive next
                const label nbrProc = twoProcs.second();

                {
                    const labelList& map = subMap[nbrProc];
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    OPstream::send(subField, nbrProc, tag, comm);
                }
                {
                    List<T> subField;
                    IPstream::recv(subField, nbrProc, tag, comm);

                    const labelList& map = constructMap[nbrProc];

                    checkReceivedSize(nbrProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        newField,
                        subField,
                        map,
                        constructHasFlip,
                        eqOp<T>(),
                        negOp
                    );
                }
            }
            else
            {
                // I am receive first, send next
                const label nbrProc = twoProcs.first();

                {
                    List<T> subField;
                    IPstream::recv(subField, nbrProc, tag, comm);

                    const labelList& map = constructMap[nbrProc];

                    checkReceivedSize(nbrProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        newField,
                        subField,
                        map,
                        constructHasFlip,
                        eqOp<T>(),
                        negOp
                    );
                }
                {
                    const labelList& map = subMap[nbrProc];
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    OPstream::send(subField, nbrProc, tag, comm);
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        const label startOfRequests = UPstream::nRequests();

        if (!is_contiguous<T>::value)
        {
            PstreamBuffers pBufs(comm, tag);

            // Stream data into buffer
            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = subMap[proci];

                if (proci != myRank && map.size())
                {
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    UOPstream os(proci, pBufs);
                    os  << subField;
                }
            }

            // Initiate receiving - do yet not block
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                List<T> subField
                (
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
                );

                // Combining bits - can now reuse field storage
                field.resize_nocopy(constructSize);

                // Receive sub field from myself
                const labelList& map = constructMap[myRank];

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    eqOp<T>(),
                    negOp
                );
            }

            // Wait for receive requests (and the send requests too)
            UPstream::waitRequests(startOfRequests);

            // Receive and process neighbour fields
            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = constructMap[proci];

                if (proci != myRank && map.size())
                {
                    UIPstream is(proci, pBufs);
                    List<T> subField(is);

                    checkReceivedSize(proci, map.size(), subField.size());

                    flipAndCombine
                    (
                        field,
                        subField,
                        map,
                        constructHasFlip,
                        eqOp<T>(),
                        negOp
                    );
                }
            }
        }
        else
        {
            // Set up receives from neighbours

            List<List<T>> recvFields(nProcs);
            DynamicList<int> recvProcs(nProcs);

            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = constructMap[proci];

                if (proci != myRank && map.size())
                {
                    recvProcs.push_back(proci);
                    List<T>& subField = recvFields[proci];
                    subField.resize_nocopy(map.size());

                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        subField.data_bytes(),
                        subField.size_bytes(),
                        tag,
                        comm
                    );
                }
            }


            // Set up sends to neighbours

            List<List<T>> sendFields(nProcs);

            for (const int proci : UPstream::allProcs(comm))
            {
                const labelList& map = subMap[proci];

                if (proci != myRank && map.size())
                {
                    List<T>& subField = sendFields[proci];
                    subField.resize_nocopy(map.size());

                    accessAndFlip(subField, field, map, subHasFlip, negOp);

                    UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        subField.cdata_bytes(),
                        subField.size_bytes(),
                        tag,
                        comm
                    );
                }
            }

            // Set up 'send' to myself - copy directly into recvFields
            {
                const labelList& map = subMap[myRank];
                List<T>& subField = recvFields[myRank];
                subField.resize_nocopy(map.size());

                accessAndFlip(subField, field, map, subHasFlip, negOp);
            }


            // Combining bits - can now reuse field storage
            field.resize_nocopy(constructSize);


            // Receive sub field from myself : recvFields[myRank]
            {
                const labelList& map = constructMap[myRank];
                const List<T>& subField = recvFields[myRank];

                // Probably don't need a size check
                // checkReceivedSize(myRank, map.size(), subField.size());

                flipAndCombine
                (
                    field,
                    subField,
                    map,
                    constructHasFlip,
                    eqOp<T>(),
                    negOp
                );
            }


            // Poll for completed receive requests and dispatch
            DynamicList<int> indices(recvProcs.size());
            while
            (
                UPstream::waitSomeRequests
                (
                    startOfRequests,
                    recvProcs.size(),
                   &indices
                )
            )
            {
                for (const int idx : indices)
                {
                    const int proci = recvProcs[idx];
                    const labelList& map = constructMap[proci];
                    const List<T>& subField = recvFields[proci];

                    // No size check - was dimensioned above
                    // checkReceivedSize(proci, map.size(), subField.size());

                    flipAndCombine
                    (
                        field,
                        subField,
                        map,
                        constructHasFlip,
                        eqOp<T>(),
                        negOp
                    );
                }
            }

            // Wait for any remaining requests
            UPstream::waitRequests(startOfRequests);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


template<class T>
void Foam::mapDistributeBase::send
(
    PstreamBuffers& pBufs,
    const List<T>& field
) const
{
    // Stream data into buffer
    for (const int proci : UPstream::allProcs(comm_))
    {
        const labelList& map = subMap_[proci];

        if (map.size())
        {
            List<T> subField
            (
                accessAndFlip(field, map, subHasFlip_, flipOp())
            );

            UOPstream os(proci, pBufs);
            os  << subField;
        }
    }

    // Start sending and receiving but do not block.
    pBufs.finishedSends(false);
}


template<class T>
void Foam::mapDistributeBase::receive
(
    PstreamBuffers& pBufs,
    List<T>& field
) const
{
    // Consume
    field.resize_nocopy(constructSize_);

    for (const int proci : UPstream::allProcs(comm_))
    {
        const labelList& map = constructMap_[proci];

        if (map.size())
        {
            UIPstream is(proci, pBufs);
            List<T> subField(is);

            checkReceivedSize(proci, map.size(), subField.size());

            flipAndCombine
            (
                field,
                subField,
                map,
                constructHasFlip_,
                eqOp<T>(),
                flipOp()
            );
        }
    }
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize_,
        subMap_,
        subHasFlip_,
        constructMap_,
        constructHasFlip_,
        values,
        negOp,
        tag,
        comm_
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    const T& nullValue,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize_,
        subMap_,
        subHasFlip_,
        constructMap_,
        constructHasFlip_,
        values,
        nullValue,
        eqOp<T>(),
        negOp,
        tag,
        comm_
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        UPstream::defaultCommsType, values, negOp, tag
    );
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    List<T>& values,
    const int tag
) const
{
    distribute(commsType, values, flipOp(), tag);
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    const UPstream::commsTypes commsType,
    DynamicList<T>& values,
    const int tag
) const
{
    values.shrink();

    List<T>& list = static_cast<List<T>&>(values);

    distribute(commsType, list, tag);

    values.setCapacity(list.size());
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    List<T>& values,
    const int tag
) const
{
    distribute(UPstream::defaultCommsType, values, tag);
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    DynamicList<T>& values,
    const int tag
) const
{
    distribute(UPstream::defaultCommsType, values, tag);
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const UPstream::commsTypes commsType,
    const label constructSize,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute<T, flipOp>
    (
        commsType,
        constructSize,
        values,
        flipOp(),
        tag
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::reverseDistribute
(
    const UPstream::commsTypes commsType,
    const label constructSize,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize,
        constructMap_,
        constructHasFlip_,
        subMap_,
        subHasFlip_,
        values,
        negOp,
        tag,
        comm_
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const UPstream::commsTypes commsType,
    const label constructSize,
    const T& nullValue,
    List<T>& values,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize,
        constructMap_,
        constructHasFlip_,
        subMap_,
        subHasFlip_,
        values,

        nullValue,
        eqOp<T>(),
        flipOp(),

        tag,
        comm_
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute
    (
        UPstream::defaultCommsType,
        constructSize,
        values,
        tag
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    const T& nullValue,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute
    (
        UPstream::defaultCommsType,
        constructSize,
        nullValue,
        values,
        tag
    );
}


// ************************************************************************* //
