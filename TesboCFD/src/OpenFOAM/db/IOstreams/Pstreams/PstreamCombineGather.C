/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    Variant of gather.
    Normal gather uses:
    - default construct and read (>>) from Istream
    - binary operator and assignment operator to combine values

    combineGather uses:
    - construct from Istream
    - modify operator which modifies its lhs

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::Pstream::combineGather
(
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // Communication order
        const auto& comms = UPstream::whichCommunication(comm);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        for (const label belowID : myComm.below())
        {
            if (is_contiguous<T>::value)
            {
                T received;

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&received),
                    sizeof(T),
                    tag,
                    comm
                );

                if (debug & 2)
                {
                    Perr<< " received from "
                        << belowID << " data:" << received << endl;
                }

                cop(value, received);
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
                T received(fromBelow);

                if (debug & 2)
                {
                    Perr<< " received from "
                        << belowID << " data:" << received << endl;
                }

                cop(value, received);
            }
        }

        // Send up value
        if (myComm.above() >= 0)
        {
            if (debug & 2)
            {
                Perr<< " sending to " << myComm.above()
                    << " data:" << value << endl;
            }

            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&value),
                    sizeof(T),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream::send(value, myComm.above(), tag, comm);
            }
        }
    }
}


template<class T, class CombineOp>
void Foam::Pstream::combineReduce
(
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        Pstream::combineGather(value, cop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::Pstream::listCombineGather
(
    UList<T>& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // Communication order
        const auto& comms = UPstream::whichCommunication(comm);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        for (const label belowID : myComm.below())
        {
            if (is_contiguous<T>::value)
            {
                List<T> received(values.size());

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    received.data_bytes(),
                    received.size_bytes(),
                    tag,
                    comm
                );

                if (debug & 2)
                {
                    Perr<< " received from "
                        << belowID << " data:" << received << endl;
                }

                forAll(values, i)
                {
                    cop(values[i], received[i]);
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
                List<T> received(fromBelow);

                if (debug & 2)
                {
                    Perr<< " received from "
                        << belowID << " data:" << received << endl;
                }

                forAll(values, i)
                {
                    cop(values[i], received[i]);
                }
            }
        }

        // Send up values
        if (myComm.above() >= 0)
        {
            if (debug & 2)
            {
                Perr<< " sending to " << myComm.above()
                    << " data:" << values << endl;
            }

            if (is_contiguous<T>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    values.cdata_bytes(),
                    values.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream::send(values, myComm.above(), tag, comm);
            }
        }
    }
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineReduce
(
    List<T>& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        Pstream::listCombineGather(values, cop, tag, comm);
        Pstream::broadcast(values, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Container, class CombineOp>
void Foam::Pstream::mapCombineGather
(
    Container& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // Communication order
        const auto& comms = UPstream::whichCommunication(comm);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        for (const label belowID : myComm.below())
        {
            // Map/HashTable: non-contiguous

            IPstream fromBelow
            (
                UPstream::commsTypes::scheduled,
                belowID,
                0,  // bufsize
                tag,
                comm
            );
            Container received(fromBelow);

            if (debug & 2)
            {
                Perr<< " received from "
                    << belowID << " data:" << received << endl;
            }

            for
            (
                auto recvIter = received.cbegin();
                recvIter != received.cend();
                ++recvIter
            )
            {
                auto masterIter = values.find(recvIter.key());

                if (masterIter.good())
                {
                    // Combine with existing
                    cop(masterIter.val(), recvIter.val());
                }
                else
                {
                    // Insert new key/value
                    values.insert(recvIter.key(), recvIter.val());
                }
            }
        }

        // Send up values
        if (myComm.above() >= 0)
        {
            if (debug & 2)
            {
                Perr<< " sending to " << myComm.above()
                    << " data:" << values << endl;
            }

            OPstream::send(values, myComm.above(), tag, comm);
        }
    }
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineReduce
(
    Container& values,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    if (UPstream::is_parallel(comm))
    {
        Pstream::mapCombineGather(values, cop, tag, comm);
        Pstream::broadcast(values, comm);
    }
}


// ************************************************************************* //
