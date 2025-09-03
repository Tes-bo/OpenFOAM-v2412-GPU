/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::Pstream::gatherList
(
    const UList<UPstream::commsStruct>& comms,
    UList<T>& values,
    const int tag,
    const label comm
)
{
    if (!comms.empty() && UPstream::is_parallel(comm))
    {
        const label myProci = UPstream::myProcNo(comm);
        const label numProc = UPstream::nProcs(comm);

        if (values.size() < numProc)
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << numProc << nl
                << Foam::abort(FatalError);
        }

        // My communication order
        const auto& myComm = comms[myProci];

        // Receive from my downstairs neighbours
        for (const label belowID : myComm.below())
        {
            const labelList& belowLeaves = comms[belowID].allBelow();

            if (is_contiguous<T>::value)
            {
                List<T> received(belowLeaves.size() + 1);

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    received.data_bytes(),
                    received.size_bytes(),
                    tag,
                    comm
                );

                values[belowID] = received[0];

                forAll(belowLeaves, leafI)
                {
                    values[belowLeaves[leafI]] = received[leafI + 1];
                }
            }
            else
            {
                IPstream fromBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,  // bufsize
                    tag,
                    comm
                );
                fromBelow >> values[belowID];

                if (debug & 2)
                {
                    Perr<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                for (const label leafID : belowLeaves)
                {
                    fromBelow >> values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() >= 0)
        {
            const labelList& belowLeaves = myComm.allBelow();

            if (debug & 2)
            {
                Perr<< " sending to " << myComm.above()
                    << " data from me:" << myProci
                    << " data:" << values[myProci] << endl;
            }

            if (is_contiguous<T>::value)
            {
                List<T> sending(belowLeaves.size() + 1);
                sending[0] = values[myProci];

                forAll(belowLeaves, leafI)
                {
                    sending[leafI + 1] = values[belowLeaves[leafI]];
                }

                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    sending.cdata_bytes(),
                    sending.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,  // bufsize
                    tag,
                    comm
                );
                toAbove << values[myProci];

                for (const label leafID : belowLeaves)
                {
                    if (debug & 2)
                    {
                        Perr<< " sending to "
                            << myComm.above() << " data from:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                    toAbove << values[leafID];
                }
            }
        }
    }
}


template<class T>
void Foam::Pstream::scatterList
(
    const UList<UPstream::commsStruct>& comms,
    UList<T>& values,
    const int tag,
    const label comm
)
{
    // Apart from the additional size check, the only difference
    // between scatterList() and using broadcast(List<T>&) or a regular
    // scatter(List<T>&) is that processor-local data is skipped.

    if (!comms.empty() && UPstream::is_parallel(comm))
    {
        const label myProci = UPstream::myProcNo(comm);
        const label numProc = UPstream::nProcs(comm);

        if (values.size() < numProc)
        {
            FatalErrorInFunction
                << "List of values:" << values.size()
                << " < numProcs:" << numProc << nl
                << Foam::abort(FatalError);
        }

        // My communication order
        const auto& myComm = comms[myProci];

        // Receive from up
        if (myComm.above() >= 0)
        {
            const labelList& notBelowLeaves = myComm.allNotBelow();

            if (is_contiguous<T>::value)
            {
                List<T> received(notBelowLeaves.size());

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    received.data_bytes(),
                    received.size_bytes(),
                    tag,
                    comm
                );

                forAll(notBelowLeaves, leafI)
                {
                    values[notBelowLeaves[leafI]] = received[leafI];
                }
            }
            else
            {
                IPstream fromAbove
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    0,  // bufsize
                    tag,
                    comm
                );

                for (const label leafID : notBelowLeaves)
                {
                    fromAbove >> values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " received through "
                            << myComm.above() << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            const label belowID = myComm.below()[belowI];
            const labelList& notBelowLeaves = comms[belowID].allNotBelow();

            if (is_contiguous<T>::value)
            {
                List<T> sending(notBelowLeaves.size());

                forAll(notBelowLeaves, leafI)
                {
                    sending[leafI] = values[notBelowLeaves[leafI]];
                }

                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    sending.cdata_bytes(),
                    sending.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,  // bufsize
                    tag,
                    comm
                );

                // Send data destined for all other processors below belowID
                for (const label leafID : notBelowLeaves)
                {
                    toBelow << values[leafID];

                    if (debug & 2)
                    {
                        Perr<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << values[leafID] << endl;
                    }
                }
            }
        }
    }
}


template<class T>
void Foam::Pstream::gatherList
(
    UList<T>& values,
    const int tag,
    const label comm
)
{
    Pstream::gatherList
    (
        UPstream::whichCommunication(comm),
        values,
        tag,
        comm
    );
}


// Unused - slate for removal? (MAY-2023)
template<class T>
void Foam::Pstream::scatterList
(
    UList<T>& values,
    const int tag,
    const label comm
)
{
    Pstream::scatterList
    (
        UPstream::whichCommunication(comm),
        values,
        tag,
        comm
    );
}


template<class T>
void Foam::Pstream::allGatherList
(
    UList<T>& values,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        if (is_contiguous<T>::value)
        {
            if (values.size() < UPstream::nProcs(comm))
            {
                FatalErrorInFunction
                    << "List of values is too small:" << values.size()
                    << " vs numProcs:" << UPstream::nProcs(comm) << nl
                    << Foam::abort(FatalError);
            }

            UPstream::mpiAllGather(values.data_bytes(), sizeof(T), comm);
            return;
        }

        const auto& comms = UPstream::whichCommunication(comm);

        Pstream::gatherList(comms, values, tag, comm);
        Pstream::scatterList(comms, values, tag, comm);
    }
}


// ************************************************************************* //
