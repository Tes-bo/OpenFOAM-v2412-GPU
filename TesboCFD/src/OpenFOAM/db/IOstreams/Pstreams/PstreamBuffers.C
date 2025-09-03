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

\*---------------------------------------------------------------------------*/

#include "PstreamBuffers.H"
#include "bitSet.H"
#include "debug.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::PstreamBuffers::algorithm
(
    // Name may change in the future (JUN-2023)
    Foam::debug::optimisationSwitch("pbufs.tuning", 0)
);
registerOptSwitch
(
    "pbufs.tuning",
    int,
    Foam::PstreamBuffers::algorithm
);

namespace Foam
{
    defineTypeNameAndDebug(PstreamBuffers, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::PstreamBuffers::setFinished(bool on) noexcept
{
    finishedSendsCalled_ = on;
}


inline void Foam::PstreamBuffers::initFinalExchange()
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    setFinished(true);

    clearUnregistered();
}


void Foam::PstreamBuffers::finalExchange
(
    enum modeOption mode,
    const bool wait,
    labelList& recvSizes
)
{
    initFinalExchange();

    // Pre-flight checks
    switch (mode)
    {
        case modeOption::DEFAULT :
        {
            // Choose (ALL_TO_ALL | NBX_PEX) from static settings
            mode =
            (
                (algorithm <= 0)
              ? modeOption::ALL_TO_ALL
              : modeOption::NBX_PEX
            );
            break;
        }

        case modeOption::GATHER :
        {
            // gather mode (all-to-one) : master [0] <- everyone
            // - only send to master [0]
            // note: master [0] is also allowed to 'send' to itself

            for (label proci = 1; proci < sendBuffers_.size(); ++proci)
            {
                sendBuffers_[proci].clear();
            }
            break;
        }

        case modeOption::SCATTER :
        {
            // scatter mode (one-to-all) : master [0] -> everyone

            if (!UPstream::master(comm_))
            {
                // Non-master: has no sends
                clearSends();
            }
            break;
        }

        default :
            break;
    }


    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        // PEX algorithm with different flavours of exchanging sizes
        // PEX stage 1: exchange sizes

        labelList sendSizes;  // Not used by gather/scatter

        switch (mode)
        {
            case modeOption::GATHER :
            {
                // gather mode (all-to-one): master [0] <- everyone
                // - presumed that MPI_Gather will be the most efficient

                recvSizes =
                    UPstream::listGatherValues(sendBuffers_[0].size(), comm_);

                if (!UPstream::master(comm_))
                {
                    recvSizes.resize_nocopy(nProcs_);
                    recvSizes = Zero;
                }

                break;
            }

            case modeOption::SCATTER :
            {
                // scatter mode (one-to-all): master [0] -> everyone
                // - presumed that MPI_Scatter will be the most efficient

                recvSizes.resize_nocopy(nProcs_);

                if (UPstream::master(comm_))
                {
                    forAll(sendBuffers_, proci)
                    {
                        recvSizes[proci] = sendBuffers_[proci].size();
                    }
                }

                const label myRecv
                (
                    UPstream::listScatterValues(recvSizes, comm_)
                );

                recvSizes = Zero;
                recvSizes[0] = myRecv;

                break;
            }

            case modeOption::NBX_PEX :
            {
                // Assemble the send sizes (cf. Pstream::exchangeSizes)
                sendSizes.resize_nocopy(nProcs_);
                forAll(sendBuffers_, proci)
                {
                    sendSizes[proci] = sendBuffers_[proci].size();
                }
                recvSizes.resize_nocopy(nProcs_);

                // Exchange sizes (non-blocking consensus)
                UPstream::allToAllConsensus
                (
                    sendSizes,
                    recvSizes,
                    (tag_ + 314159),  // some unique tag?
                    comm_
                );
                break;
            }

            case modeOption::DEFAULT :
            case modeOption::ALL_TO_ALL :
            {
                // Assemble the send sizes (cf. Pstream::exchangeSizes)
                sendSizes.resize_nocopy(nProcs_);
                forAll(sendBuffers_, proci)
                {
                    sendSizes[proci] = sendBuffers_[proci].size();
                }
                recvSizes.resize_nocopy(nProcs_);

                // Exchange sizes (all-to-all)
                UPstream::allToAll(sendSizes, recvSizes, comm_);
                break;
            }
        }


        // PEX stage 2: point-to-point data exchange
        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


void Foam::PstreamBuffers::finalExchange
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const bool wait,
    labelList& recvSizes
)
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;

    initFinalExchange();

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        // Preparation. Temporarily abuse recvSizes as logic to clear
        // send buffers that are not in the neighbourhood connection
        {
            recvSizes.resize_nocopy(nProcs_);
            recvSizes = 0;

            // Preserve self-send, even if not described by neighbourhood
            recvSizes[UPstream::myProcNo(comm_)] = 1;

            for (const label proci : sendProcs)
            {
                recvSizes[proci] = 1;  // Connected
            }

            for (label proci = 0; proci < nProcs_; ++proci)
            {
                if (!recvSizes[proci])  // Not connected
                {
                    sendBuffers_[proci].clear();
                }
            }
        }

        // PEX stage 1: exchange sizes (limited neighbourhood)
        Pstream::exchangeSizes
        (
            sendProcs,
            recvProcs,
            sendBuffers_,
            recvSizes,
            tag_,
            comm_
        );

        // PEX stage 2: point-to-point data exchange
        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    UPstream::commsTypes commsType,
    int tag,
    label communicator,
    IOstreamOption::streamFormat fmt
)
:
    finishedSendsCalled_(false),
    allowClearRecv_(true),
    format_(fmt),
    commsType_(commsType),
    tag_(tag),
    comm_(communicator),
    nProcs_(UPstream::nProcs(comm_)),
    sendBuffers_(nProcs_),
    recvBuffers_(nProcs_),
    recvPositions_(nProcs_, Zero)
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;

    // Check that all data has been consumed.
    forAll(recvBuffers_, proci)
    {
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos >= 0 && pos < len)
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " Only consumed " << pos << " of " << len << " bytes" << nl
                << " comm " << comm_ << " tag " << tag_ << nl
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::DynamicList<char>& Foam::PstreamBuffers::accessSendBuffer
(
    const label proci
)
{
    return sendBuffers_[proci];
}


Foam::DynamicList<char>& Foam::PstreamBuffers::accessRecvBuffer
(
    const label proci
)
{
    return recvBuffers_[proci];
}


Foam::label& Foam::PstreamBuffers::accessRecvPosition(const label proci)
{
    return recvPositions_[proci];
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::clearSends()
{
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clear();
    }
}


void Foam::PstreamBuffers::clearRecvs()
{
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clear();
    }
    recvPositions_ = Zero;
}


void Foam::PstreamBuffers::clear()
{
    clearSends();
    clearRecvs();
    setFinished(false);
}


void Foam::PstreamBuffers::clearUnregistered()
{
    for (label proci = 0; proci < nProcs_; ++proci)
    {
        if (recvPositions_[proci] < 0)
        {
            recvPositions_[proci] = 0;
            sendBuffers_[proci].clear();
        }
    }
}


void Foam::PstreamBuffers::clearSend(const label proci)
{
    sendBuffers_[proci].clear();
    if (recvPositions_[proci] < 0)
    {
        // Reset the unregistered flag
        recvPositions_[proci] = 0;
    }
}


void Foam::PstreamBuffers::clearRecv(const label proci)
{
    recvBuffers_[proci].clear();
    recvPositions_[proci] = 0;
}


void Foam::PstreamBuffers::clearStorage()
{
    // Could also clear out entire sendBuffers_, recvBuffers_ and reallocate.
    // Not sure if it makes much difference
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clearStorage();
    }
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clearStorage();
    }
    recvPositions_ = Zero;

    setFinished(false);
}


void Foam::PstreamBuffers::initRegisterSend()
{
    if (!finished())
    {
        for (label proci = 0; proci < nProcs_; ++proci)
        {
            sendBuffers_[proci].clear();
            // Mark send buffer as 'unregistered'
            recvPositions_[proci] = -1;
        }
    }
}


void Foam::PstreamBuffers::registerSend(const label proci, const bool toggleOn)
{
    // Clear the 'unregistered' flag
    if (toggleOn && recvPositions_[proci] < 0)
    {
        recvPositions_[proci] = 0;
    }
}


bool Foam::PstreamBuffers::hasSendData() const
{
    for (const DynamicList<char>& buf : sendBuffers_)
    {
        if (!buf.empty())
        {
            return true;
        }
    }
    return false;
}


bool Foam::PstreamBuffers::hasRecvData() const
{
    if (finished())
    {
        forAll(recvBuffers_, proci)
        {
            if (recvPositions_[proci] < recvBuffers_[proci].size())
            {
                return true;
            }
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return false;
}


Foam::label Foam::PstreamBuffers::sendDataCount(const label proci) const
{
    return sendBuffers_[proci].size();
}


Foam::label Foam::PstreamBuffers::recvDataCount(const label proci) const
{
    if (finished())
    {
        const label len(recvBuffers_[proci].size() - recvPositions_[proci]);

        if (len > 0)
        {
            return len;
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return 0;
}


Foam::labelList Foam::PstreamBuffers::recvDataCounts() const
{
    labelList counts(nProcs_, Zero);

    if (finished())
    {
        forAll(recvBuffers_, proci)
        {
            const label len(recvBuffers_[proci].size() - recvPositions_[proci]);

            if (len > 0)
            {
                counts[proci] = len;
            }
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return counts;
}


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount
(
    const label excludeProci
) const
{
    label maxLen = 0;

    if (finished())
    {
        forAll(recvBuffers_, proci)
        {
            if (excludeProci != proci)
            {
                label len(recvBuffers_[proci].size() - recvPositions_[proci]);
                maxLen = Foam::max(maxLen, len);
            }
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return maxLen;
}


Foam::label Foam::PstreamBuffers::maxRecvCount() const
{
    // Use out-of-range proci to avoid excluding any processor
    return maxNonLocalRecvCount(-1);
}


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount() const
{
    return maxNonLocalRecvCount(UPstream::myProcNo(comm_));
}


const Foam::UList<char>
Foam::PstreamBuffers::peekRecvData(const label proci) const
{
    if (finished())
    {
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos < len)
        {
            return UList<char>
            (
                const_cast<char*>(recvBuffers_[proci].cbegin(pos)),
                (len - pos)
            );
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return UList<char>();
}


void Foam::PstreamBuffers::finishedSends(const bool wait)
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;

    labelList recvSizes;
    finalExchange(modeOption::DEFAULT, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedSendsNBX(const bool wait)
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;

    labelList recvSizes;
    finalExchange(modeOption::NBX_PEX, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedSends
(
    labelList& recvSizes,
    const bool wait
)
{
    DebugPoutInFunction
        << "tag:" << tag_
        << " comm:" << comm_
        << " nProcs:" << nProcs_
        << endl;

    // Resize for copying back
    recvSizes.resize_nocopy(sendBuffers_.size());

    finalExchange(modeOption::DEFAULT, wait, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(neighProcs, neighProcs, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    const bool wait
)
{
    labelList recvSizes;
    finalExchange(neighProcs, neighProcs, wait, recvSizes);
}


bool Foam::PstreamBuffers::finishedSends
(
    bitSet& sendConnections,
    DynamicList<label>& sendProcs,
    DynamicList<label>& recvProcs,
    const bool wait
)
{
    bool changed = (sendConnections.size() != nProcs());

    if (changed)
    {
        sendConnections.resize(nProcs());
    }

    // Update send connections
    forAll(sendBuffers_, proci)
    {
        if (sendConnections.set(proci, !sendBuffers_[proci].empty()))
        {
            // The state changed
            changed = true;
        }
    }

    UPstream::reduceOr(changed, comm_);

    if (changed)
    {
        // Update send/recv topology

        labelList recvSizes;
        finishedSends(recvSizes, wait);  // modeOption::DEFAULT (eg all-to-all)

        // The send ranks
        sendProcs.clear();
        forAll(sendBuffers_, proci)
        {
            if (!sendBuffers_[proci].empty())
            {
                sendProcs.push_back(proci);
            }
        }

        // The recv ranks
        recvProcs.clear();
        forAll(recvSizes, proci)
        {
            if (recvSizes[proci] > 0)
            {
                recvProcs.push_back(proci);
            }
        }
    }
    else
    {
        // Use existing send/recv ranks
        labelList recvSizes;
        finalExchange(sendProcs, recvProcs, wait, recvSizes);
    }

    return changed;
}


void Foam::PstreamBuffers::finishedGathers(const bool wait)
{
    labelList recvSizes;
    finalExchange(modeOption::GATHER, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedScatters(const bool wait)
{
    labelList recvSizes;
    finalExchange(modeOption::SCATTER, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedGathers
(
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(modeOption::GATHER, wait, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::finishedScatters
(
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(modeOption::SCATTER, wait, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Controls

bool Foam::PstreamBuffers::finished() const noexcept
{
    return finishedSendsCalled_;
}


bool Foam::PstreamBuffers::allowClearRecv() const noexcept
{
    return allowClearRecv_;
}


bool Foam::PstreamBuffers::allowClearRecv(bool on) noexcept
{
    bool old(allowClearRecv_);
    allowClearRecv_ = on;
    return old;
}


// ************************************************************************* //
