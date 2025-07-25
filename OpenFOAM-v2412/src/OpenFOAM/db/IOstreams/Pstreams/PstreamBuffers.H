/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::PstreamBuffers

Description
    Buffers for inter-processor communications streams (UOPstream, UIPstream).

    Use UOPstream to stream data into buffers, call finishedSends() to
    notify that data is in buffers and then use IUPstream to get data out
    of received buffers. Works with both buffered and non-blocking. Does
    not make much sense with scheduled since there you would not need these
    explicit buffers.

    Example usage:
    \code
        PstreamBuffers pBufs;

        for (const int proci : UPstream::allProcs())
        {
            if (proci != UPstream::myProcNo())
            {
                someObject vals;

                UOPstream send(proci, pBufs);
                send << vals;
            }
        }

        pBufs.finishedSends();   // no-op for buffered

        for (const int proci : UPstream::allProcs())
        {
            if (proci != UPstream::myProcNo())
            {
                UIPstream recv(proci, pBufs);
                someObject vals(recv);
            }
        }
    \endcode

    There are special versions of finishedSends() for
    restricted neighbour communication as well as for special
    one-to-all and all-to-one communication patterns.
    For example,
    \code
        PstreamBuffers pBufs;

        if (UPstream::master())
        {
            someObject vals;
            for (const int proci : UPstream::subProcs())
            {
                UOPstream send(proci, pBufs);
                send << vals;
            }
        }

        pBufs.finishedScatters();

        if (!UPstream::master())
        {
            UIPstream recv(UPstream::masterNo(), pBufs);
            someObject vals(recv);
        }
    \endcode

    Additionally there are some situations that use speculative sends
    that may not actually be required. In this case, it is possible to
    mark all sends as initially \em unregistered and subsequently
    mark the "real" sends as \em registered.

    For example,
    \code
        PstreamBuffers pBufs;

        pBufs.initRegisterSend();

        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);
            if (ppp)
            {
                const label nbrProci = ppp->neighbProcNo();

                // Gather some patch information...
                UOPstream toNbr(nbrProci, pBufs);
                toNbr << patchInfo;

                // The send is needed if patchInfo is non-empty
                pBufs.registerSend(nbrProci, !patchInfo.empty());
            }
        }

        // optional:  pBufs.clearUnregistered();

        pBufs.finishedSends();

        ...
    \endcode

SourceFiles
    PstreamBuffers.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_PstreamBuffers_H
#define Foam_PstreamBuffers_H

#include "DynamicList.H"
#include "UPstream.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                       Class PstreamBuffers Declaration
\*---------------------------------------------------------------------------*/

class PstreamBuffers
{
    // Private Data Types

        //- Private enumeration for handling PEX stage 1 (sizing) modes
        enum class modeOption : unsigned char
        {
            DEFAULT,    //!< Use statically configured algorithm
            GATHER,     //!< Use all-to-one (gather) of sizes
            SCATTER,    //!< Use one-to-all (scatter) of sizes
            ALL_TO_ALL, //!< Use allToAll to obtain sizes
            NBX_PEX     //!< Use consensus exchange (NBX) to obtain sizes
        };


    // Private Data

        //- Track if sends are complete
        bool finishedSendsCalled_;

        //- Permit clear of individual receive buffer by external access
        bool allowClearRecv_;

        //- Buffer format (ascii | binary)
        const IOstreamOption::streamFormat format_;

        //- Communications type of this stream
        const UPstream::commsTypes commsType_;

        //- The transfer message type
        const int tag_;

        //- Communicator
        const label comm_;

        //- Number of ranks associated with PstreamBuffers (at construction)
        const label nProcs_;


    // Buffer storage

        //- Send buffers. Size is nProcs()
        List<DynamicList<char>> sendBuffers_;

        //- Receive buffers. Size is nProcs()
        List<DynamicList<char>> recvBuffers_;

        //- Current read positions within recvBuffers_. Size is nProcs()
        //  This list is also misused for registerSend() bookkeeping
        labelList recvPositions_;


    // Private Member Functions

        //- Change status of finished sends called
        inline void setFinished(bool on) noexcept;

        //- Clear 'unregistered' send buffers, tag as being send-ready
        inline void initFinalExchange();

        //- Mark all sends as having been done.
        //  This will start receives (non-blocking comms).
        void finalExchange
        (
            enum modeOption mode,
            const bool wait,
            labelList& recvSizes
        );

        //- Mark sends as done.
        //  Only exchange sizes using the neighbour ranks
        //  (non-blocking comms).
        void finalExchange
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const bool wait,
            labelList& recvSizes
        );


    // Friendship Access

        //- Access a send buffer for given proc (in range 0-nProcs)
        DynamicList<char>& accessSendBuffer(const label proci);

        //- Access a recv buffer for given proc (in range 0-nProcs).
        DynamicList<char>& accessRecvBuffer(const label proci);

        //- Access the recv position within recv buffer for given proc
        //- (in range 0-nProcs).
        label& accessRecvPosition(const label proci);

        friend class UOPstreamBase;  // accessSendBuffer()
        friend class UIPstreamBase;  // accessRecvBuffer(), accessRecvPosition()


public:

    // Declare name of the class and its debug switch
    ClassName("PstreamBuffers");


    // Static Data

        //- Preferred exchange algorithm (may change or be removed in future)
        static int algorithm;


    // Constructors

        //- Construct given communication type (default: nonBlocking), message
        //- tag, communicator (default: worldComm), IO format (default: binary)
        explicit PstreamBuffers
        (
            UPstream::commsTypes commsType = UPstream::commsTypes::nonBlocking,
            int tag = UPstream::msgType(),
            label communicator = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given communicator, communication type
        //- (default: nonBlocking), message tag, IO format (default: binary)
        explicit PstreamBuffers
        (
            label communicator,
            UPstream::commsTypes commsType = UPstream::commsTypes::nonBlocking,
            int tag = UPstream::msgType(),
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            PstreamBuffers(commsType, tag, communicator, fmt)
        {}

        //- Construct given communicator, message tag, communication type
        //- (default: nonBlocking), IO format (default: binary)
        PstreamBuffers
        (
            label communicator,
            int tag,
            UPstream::commsTypes commsType = UPstream::commsTypes::nonBlocking,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            PstreamBuffers(commsType, tag, communicator, fmt)
        {}


    //- Destructor - checks that all data have been consumed
    ~PstreamBuffers();


    // Member Functions

    // Attributes

        //- The associated buffer format (ascii | binary)
        IOstreamOption::streamFormat format() const noexcept
        {
            return format_;
        }

        //- The communications type of the stream
        UPstream::commsTypes commsType() const noexcept
        {
            return commsType_;
        }

        //- The transfer message tag
        int tag() const noexcept
        {
            return tag_;
        }

        //- The communicator index
        label comm() const noexcept
        {
            return comm_;
        }

        //- Number of ranks associated with PstreamBuffers
        label nProcs() const noexcept
        {
            return nProcs_;
        }


    // Sizing

        //- Range of ranks indices associated with PstreamBuffers
        UPstream::rangeType allProcs() const noexcept
        {
            // Proc 0 -> nProcs (int value)
            return UPstream::rangeType(static_cast<int>(nProcs_));
        }

        //- Range of sub-processes indices associated with PstreamBuffers
        UPstream::rangeType subProcs() const noexcept
        {
            // Proc 1 -> nProcs (int value)
            return UPstream::rangeType(1, static_cast<int>(nProcs_-1));
        }


    // Queries

        //- True if finishedSends() or finishedNeighbourSends() has been called
        bool finished() const noexcept;

        //- Is clearStorage of individual receive buffer by external hooks
        //- allowed? (default: true)
        bool allowClearRecv() const noexcept;

        //- True if any (local) send buffers have data
        bool hasSendData() const;

        //- True if any (local) recv buffers have unconsumed data.
        //- Must call finishedSends() or other finished.. method first!
        bool hasRecvData() const;

        //- Number of send bytes for the specified processor.
        label sendDataCount(const label proci) const;

        //- Number of unconsumed receive bytes for the specified processor.
        //- Must call finishedSends() or other finished.. method first!
        label recvDataCount(const label proci) const;

        //- Number of unconsumed receive bytes for all processors.
        //- Must call finishedSends() or other finished.. method first!
        labelList recvDataCounts() const;

        //- Maximum receive size from any rocessor rank.
        //- Must call finishedSends() or other finished.. method first!
        label maxRecvCount() const;

        //- Maximum receive size, excluding current processor rank
        //- Must call finishedSends() or other finished.. method first!
        label maxNonLocalRecvCount() const;

        //- Maximum receive size, excluding the specified processor rank
        //- Must call finishedSends() or other finished.. method first!
        label maxNonLocalRecvCount(const label excludeProci) const;

        //- Number of unconsumed receive bytes for the specified processor.
        //- Must call finishedSends() or other finished.. method first!
        //  The method is only useful in limited situations, such as when
        //  PstreamBuffers has been used to fill contiguous data
        //  (eg, using OPstream::write).
        const UList<char> peekRecvData(const label proci) const;


    // Edit

        //- Clear all send/recv buffers and reset states.
        //  Does not remove the buffer storage.
        void clear();

        //- Clear all send buffers (does not remove buffer storage)
        void clearSends();

        //- Clear all recv buffer and positions (does not remove buffer storage)
        void clearRecvs();

        //- Clear an individual send buffer (eg, data not required)
        void clearSend(const label proci);

        //- Clear an individual receive buffer (eg, data not required)
        //  Does not remove the buffer storage.
        void clearRecv(const label proci);

        //- Clear storage for all send/recv buffers and reset states.
        void clearStorage();

        //- Change allowClearRecv, return previous value
        bool allowClearRecv(bool on) noexcept;


    // Registered Sending

        //- Initialise registerSend() bookkeeping by mark all send buffers
        //- as 'unregistered'
        //  Usually called immediately after construction or clear().
        void initRegisterSend();

        //- Toggle an individual send buffer as 'registered'.
        //- The setting is sticky (does not turn off)
        void registerSend(const label proci, const bool toggleOn = true);

        //- Clear any 'unregistered' send buffers.
        void clearUnregistered();


    // Regular Functions

        //- Mark the send phase as being finished.
        //
        //  Non-blocking mode: populates receive buffers using all-to-all
        //  or NBX (depending on tuning parameters).
        //  \param wait wait for requests to complete (in non-blocking mode)
        void finishedSends(const bool wait = true);

        //- Mark the send phase as being finished.
        //
        //  Non-blocking mode: populates receive buffers using NBX.
        //  \param wait wait for requests to complete (in non-blocking mode)
        void finishedSendsNBX(const bool wait = true);

        //- Mark the send phase as being finished.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers using all-to-all
        //  or NBX (depending on tuning parameters).
        //  \warning currently only valid for non-blocking comms.
        void finishedSends
        (
            //! [out] the sizes (bytes) received
            labelList& recvSizes,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );

        //- Mark the send phase as being finished.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers using NBX.
        //  \warning currently only valid for non-blocking comms.
        void finishedSendsNBX
        (
            //! [out] the sizes (bytes) received
            labelList& recvSizes,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );


    // Functions with restricted neighbours

        //- Mark the send phase as being finished, with communication
        //- being limited to a known subset of send/recv ranks.
        //
        //  Non-blocking mode: populates receive buffers.
        //
        //  \warning currently only valid for non-blocking comms.
        //  \note Same as finishedSends with identical sendProcs/recvProcs
        void finishedNeighbourSends
        (
            //! ranks used for sends/recvs
            const labelUList& neighProcs,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );

        //- Mark the send phase as being finished, with communication
        //- being limited to a known subset of send/recv ranks.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: it will populate receive buffers.
        //
        //  \warning currently only valid for non-blocking mode.
        void finishedNeighbourSends
        (
            //! ranks used for sends/recvs
            const labelUList& neighProcs,
            //! [out] the sizes (bytes) received
            labelList& recvSizes,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );

        //- A caching version that uses a limited send/recv connectivity.
        //
        //  Non-blocking mode: populates receive buffers.
        //  \return True if the send/recv connectivity changed
        //
        //  \warning currently only valid for non-blocking comms.
        bool finishedSends
        (
            //! inter-rank connections (on/off) for sending ranks
            bitSet& sendConnections,
            //! ranks used for sends
            DynamicList<label>& sendProcs,
            //! ranks used for recvs
            DynamicList<label>& recvProcs,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );


    // Gather/scatter modes

        //- Mark all sends to master as done.
        //
        //  Non-blocking mode: populates receive buffers.
        //  Can use recvDataCount, maxRecvCount etc to recover sizes received.
        //
        //  \param wait wait for requests to complete (in non-blocking mode)
        //
        //  \warning currently only valid for non-blocking comms.
        void finishedGathers(const bool wait = true);

        //- Mark all sends to master as done.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers (all-to-one).
        //  \warning currently only valid for non-blocking comms.
        void finishedGathers
        (
            //! [out] the sizes (bytes) received
            labelList& recvSizes,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );

        //- Mark all sends to sub-procs as done.
        //
        //  Non-blocking mode: populates receive buffers.
        //  Can use recvDataCount, maxRecvCount etc to recover sizes received.
        //
        //  \param wait wait for requests to complete (in non-blocking mode)
        //
        //  \warning currently only valid for non-blocking comms.
        void finishedScatters(const bool wait = true);

        //- Mark all sends to sub-procs as done.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers (all-to-one).
        //  \warning currently only valid for non-blocking comms.
        void finishedScatters
        (
            //! [out] the sizes (bytes) received
            labelList& recvSizes,
            //! wait for requests to complete (in non-blocking mode)
            const bool wait = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
