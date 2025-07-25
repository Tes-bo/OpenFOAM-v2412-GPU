/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Class
    Foam::Pstream

Description
    Inter-processor communications stream.

SourceFiles
    Pstream.C
    PstreamBroadcast.C
    PstreamGather.C
    PstreamCombineGather.C
    PstreamGatherList.C
    PstreamExchangeConsensus.C
    PstreamExchange.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Pstream_H
#define Foam_Pstream_H

#include "UPstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Pstream Declaration
\*---------------------------------------------------------------------------*/

class Pstream
:
    public UPstream
{
protected:

    // Protected Data

        //- Allocated transfer buffer (can be used for send or receive)
        DynamicList<char> transferBuf_;


public:

    //- Declare name of the class and its debug switch
    ClassName("Pstream");


    // Constructors

        //- Construct for given communication type, with optional buffer size
        explicit Pstream
        (
            const UPstream::commsTypes commsType,
            const label bufSize = 0
        )
        :
            UPstream(commsType)
        {
            if (bufSize > 0)
            {
                transferBuf_.setCapacity(bufSize + 2*sizeof(scalar) + 1);
            }
        }


    // Static Functions

    // Broadcast

        //- Broadcast buffer content to all processes in communicator.
        using UPstream::broadcast;

        //- Broadcast content (contiguous or non-contiguous) to all
        //- communicator ranks. Does nothing in \b non-parallel.
        template<class Type>
        static void broadcast
        (
            Type& value,
            const label comm = UPstream::worldComm
        );

        //- Broadcast multiple items to all communicator ranks.
        //- Does nothing in \b non-parallel.
        template<class Type, class... Args>
        static void broadcasts(const label comm, Type& arg1, Args&&... args);

        //- Broadcast list content (contiguous or non-contiguous) to all
        //- communicator ranks. Does nothing in \b non-parallel.
        //  For contiguous list data, this avoids serialization overhead,
        //  but at the expense of an additional broadcast call.
        template<class ListType>
        static void broadcastList
        (
            ListType& list,
            const label comm = UPstream::worldComm
        );


    // Gather

        //- Gather (reduce) data, applying \c bop to combine \c value
        //- from different processors. The basis for Foam::reduce().
        //  Uses linear/tree communication (with parallel guard).
        template<class T, class BinaryOp>
        static void gather
        (
            T& value,
            const BinaryOp& bop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Gather individual values into list locations.
        //  On master list length == nProcs, otherwise zero length.
        //  \n
        //  For \b non-parallel :
        //  the returned list length is 1 with localValue.
        template<class T>
        static List<T> listGatherValues
        (
            const T& localValue,
            const label comm = UPstream::worldComm,
            //! Only used for non-contiguous types
            const int tag = UPstream::msgType()
        );

        //- Scatter individual values from list locations.
        //  On master input list length == nProcs, ignored on other procs.
        //  \n
        //  For \b non-parallel :
        //  returns the first list element (or default initialized).
        template<class T>
        static T listScatterValues
        (
            const UList<T>& allValues,
            const label comm = UPstream::worldComm,
            //! Only used for non-contiguous types
            const int tag = UPstream::msgType()
        );


    // Gather/combine data
    // Inplace combine values from processors.
    // (Uses construct from Istream instead of \c << operator)

        //- Gather data, applying \c cop to inplace combine \c value
        //- from different processors.
        //  Uses linear/tree communication (with parallel guard).
        template<class T, class CombineOp>
        static void combineGather
        (
            //! [in,out]
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Reduce inplace (cf. MPI Allreduce)
        //- applying \c cop to inplace combine \c value
        //- from different processors.
        //- After completion all processors have the same data.
        //  Uses linear/tree communication.
        //  Wraps combineGather/broadcast (may change in the future).
        template<class T, class CombineOp>
        static void combineReduce
        (
            //! [in,out]
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Same as Pstream::combineReduce
        template<class T, class CombineOp>
        static void combineAllGather
        (
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            Pstream::combineReduce(value, cop, tag, comm);
        }


    // Combine variants working on whole List at a time.

        //- Combines List elements.
        //  Uses linear/tree communication (with parallel guard).
        template<class T, class CombineOp>
        static void listCombineGather
        (
            //! [in,out]
            UList<T>& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Combines List elements.
        //- After completion all processors have the same data.
        //  Uses linear/tree communication (with parallel guard).
        template<class T, class CombineOp>
        static void listCombineReduce
        (
            //! [in,out] - List (not UList) due to broadcast()
            List<T>& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Same as Pstream::listCombineReduce
        template<class T, class CombineOp>
        static void listCombineAllGather
        (
            //! [in,out] - List (not UList) due to broadcast()
            List<T>& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            Pstream::listCombineReduce(values, cop, tag, comm);
        }


    // Combine variants working on whole map at a time.
    // Container needs iterators, find() and insert methods defined.

        //- Combine Map elements.
        //  Uses linear/tree communication (with parallel guard).
        template<class Container, class CombineOp>
        static void mapCombineGather
        (
            Container& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Reduce inplace (cf. MPI Allreduce)
        //- applying \c cop to inplace combine map \c values
        //- from different processors.
        //- After completion all processors have the same data.
        //  Uses the specified communication schedule.
        //  Wraps mapCombineGather/broadcast (may change in the future).
        //- After completion all processors have the same data.
        template<class Container, class CombineOp>
        static void mapCombineReduce
        (
            Container& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Same as Pstream::mapCombineReduce
        template<class Container, class CombineOp>
        static void mapCombineAllGather
        (
            Container& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        )
        {
            Pstream::mapCombineReduce(values, cop, tag, comm);
        }


    // Gather/scatter keeping the individual processor data separate.
    // The values is a List of size UPstream::nProcs() where
    // values[UPstream::myProcNo()] is the data for the current processor.

        //- Gather data, but keep individual values separate.
        //- Uses the specified communication schedule.
        template<class T>
        static void gatherList
        (
            const UList<commsStruct>& comms,
            //! [in,out]
            UList<T>& values,
            const int tag,
            const label comm
        );

        //- Gather data, but keep individual values separate.
        //- Uses linear/tree communication.
        template<class T>
        static void gatherList
        (
            //! [in,out]
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Gather data, but keep individual values separate.
        //- Uses MPI_Allgather or manual linear/tree communication.
        //  After completion all processors have the same data.
        //  Wraps gatherList/scatterList (may change in the future).
        template<class T>
        static void allGatherList
        (
            //! [in,out]
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Scatter

        //- Inverse of gatherList.
        //- Uses the specified communication schedule.
        template<class T>
        static void scatterList
        (
            const UList<commsStruct>& comms,
            UList<T>& values,
            const int tag,
            const label comm
        );

        //- Inverse of gatherList.
        //- Uses linear/tree communication.
        template<class T>
        static void scatterList
        (
            UList<T>& values,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Exchange

        //- Helper: exchange sizes of sendBufs for specified send/recv ranks
        template<class Container>
        static void exchangeSizes
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const Container& sendBufs,
            labelList& sizes,
            const label tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Helper: exchange sizes of sendBufs for specified neighbour ranks
        template<class Container>
        static void exchangeSizes
        (
            const labelUList& neighProcs,
            const Container& sendBufs,
            labelList& sizes,
            const label tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Helper: exchange sizes of sendBufs.
        //- The sendBufs is the data per processor (in the communicator).
        //  Returns sizes of sendBufs on the sending processor.
        //  \n
        //  For \b non-parallel : copy sizes from sendBufs directly.
        template<class Container>
        static void exchangeSizes
        (
            const Container& sendBufs,
            labelList& recvSizes,
            const label comm = UPstream::worldComm
        );

        //- Exchange the \b non-zero sizes of sendBufs entries (sparse map)
        //- with other ranks in the communicator
        //- using non-blocking consensus exchange.
        //
        //  Since the recvData map always cleared before receipt and sizes
        //  of zero are never transmitted, a simple check
        //  of its keys is sufficient to determine connectivity.
        //
        //  For \b non-parallel : copy size of rank (if it exists and non-empty)
        //  from sendBufs to recvSizes.
        //
        //  \note The message tag is adjusted internally to improve uniqueness
        template<class Container>
        static void exchangeSizes
        (
            const Map<Container>& sendBufs,
            Map<label>& recvSizes,
            const label tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Helper: exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs using predetermined receive
        //- sizing.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const UList<Container>& sendBufs,
            const labelUList& recvSizes,
            List<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //  Data provided and received as container.
        //
        // No internal guards or resizing.
        template<class Container, class Type>
        static void exchange
        (
            const Map<Container>& sendBufs,
            const Map<label>& recvSizes,  //!< Num of recv elements (not bytes)
            Map<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //- Determines sizes to receive.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const UList<Container>& sendBufs,
            List<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );

        //- Exchange \em contiguous data.
        //- Sends sendBufs, receives into recvBufs.
        //- Determines sizes to receive.
        //  If wait=true will wait for all transfers to finish.
        template<class Container, class Type>
        static void exchange
        (
            const Map<Container>& sendBufs,
            Map<Container>& recvBufs,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool wait = true  //!< Wait for requests to complete
        );


    // Non-blocking exchange

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData, receives into recvData.
        //
        //  Each entry of the recvBufs list is cleared before receipt.
        //  For \b non-parallel : copy own rank from sendBufs to recvBufs.
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static void exchangeConsensus
        (
            const UList<Container>& sendBufs,
            List<Container>& recvBufs,
            const int tag,
            const label comm,
            const bool wait = true  //!< (ignored)
        );

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData, receives into recvData.
        //
        //  Each \em entry of the recvBufs map is cleared before receipt,
        //  but the map itself if not cleared. This allows the map to preserve
        //  allocated space (eg DynamicList entries) between calls.
        //
        //  For \b non-parallel : copy own rank (if it exists and non-empty)
        //  from sendBufs to recvBufs.
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static void exchangeConsensus
        (
            const Map<Container>& sendBufs,
            Map<Container>& recvBufs,
            const int tag,
            const label comm,
            const bool wait = true  //!< (ignored)
        );

        //- Exchange \em contiguous data using non-blocking consensus (NBX)
        //- Sends sendData returns receive information.
        //
        //  For \b non-parallel : copy own rank (if it exists and non-empty)
        //
        //  \note The message tag should be chosen to be a unique value
        //  since the implementation uses probing with ANY_SOURCE !!
        template<class Container, class Type>
        static Map<Container> exchangeConsensus
        (
            const Map<Container>& sendBufs,
            const int tag,
            const label comm,
            const bool wait = true  //!< (ignored)
        );


    // Housekeeping

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void scatter
        (
            T& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const label comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void combineScatter
        (
            T& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const label comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class T>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void listCombineScatter
        (
            List<T>& value,
            const int tag = UPstream::msgType(),  //!< ignored
            const label comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(value, comm);
        }

        //- \deprecated(2024-01) Broadcast data
        template<class Container>
        FOAM_DEPRECATED_FOR(2024-01, "Pstream::broadcast()")
        static void mapCombineScatter
        (
            Container& values,
            const int tag = UPstream::msgType(),  //!< ignored
            const label comm = UPstream::worldComm
        )
        {
            Pstream::broadcast(values, comm);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PstreamBroadcast.C"
    #include "PstreamGather.C"
    #include "PstreamCombineGather.C"
    #include "PstreamGatherList.C"
    #include "PstreamExchange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
