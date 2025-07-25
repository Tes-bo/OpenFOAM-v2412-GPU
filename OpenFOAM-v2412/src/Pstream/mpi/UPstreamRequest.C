/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "UPstreamWrapping.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Request::Request() noexcept
:
    UPstream::Request(MPI_REQUEST_NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Request::good() const noexcept
{
    return MPI_REQUEST_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Request::reset() noexcept
{
    *this = UPstream::Request(MPI_REQUEST_NULL);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// Foam::UPstream::Request
// Foam::UPstream::Request::lookup(const label req)
// {
//     if (req < 0 || req >= PstreamGlobals::outstandingRequests_.size())
//     {
//         WarningInFunction
//             << "Illegal request " << req << nl
//             << "Should be within range [0,"
//             << PstreamGlobals::outstandingRequests_.size()
//             << ')' << endl;
//
//         return UPstream::Communicator(MPI_REQUEST_NULL);
//     }
//
//     return UPstream::Request(PstreamGlobals::outstandingRequests_[req]);
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UPstream::nRequests() noexcept
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label n)
{
    if (n >= 0 && n < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.resize(n);
    }
}


void Foam::UPstream::addRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    // Transcribe as a MPI_Request
    PstreamGlobals::outstandingRequests_.push_back
    (
        PstreamUtils::Cast::to_mpi(req)
    );

    // Invalidate parameter
    req = UPstream::Request(MPI_REQUEST_NULL);
}


void Foam::UPstream::cancelRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return;
    }

    {
        auto& request = PstreamGlobals::outstandingRequests_[i];
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);  //<- Sets to MPI_REQUEST_NULL
        }
    }
}


void Foam::UPstream::cancelRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::cancelRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::removeRequests(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    const labelRange range(pos, count);

    for (const label i : range)
    {
        auto& request = PstreamGlobals::outstandingRequests_[i];
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);  //<- Sets to MPI_REQUEST_NULL
        }
    }

    // Remove from list of outstanding requests and move down
    PstreamGlobals::outstandingRequests_.remove(range);
}


void Foam::UPstream::freeRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            // if (cancel)
            // {
            //     MPI_Cancel(&request);
            // }
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::freeRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            // if (cancel)
            // {
            //     MPI_Cancel(&request);
            // }
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::waitRequests(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);
    bool trim = true;  // Can trim the trailing part of the list

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
        trim = false;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequests : starting wait for "
            << count << " requests starting at " << pos << endl;
    }

    profilingPstream::beginTiming();

    if (count == 1)
    {
        // On success: sets request to MPI_REQUEST_NULL
        if (MPI_Wait(waitRequests, MPI_STATUS_IGNORE))
        {
            FatalErrorInFunction
                << "MPI_Wait returned with error"
                << Foam::abort(FatalError);
        }
    }
    else if (count > 1)
    {
        // On success: sets each request to MPI_REQUEST_NULL
        if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error"
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addWaitTime();

    if (trim)
    {
        // Trim the length of outstanding requests
        PstreamGlobals::outstandingRequests_.resize(pos);
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        if (MPI_REQUEST_NULL != request)  // Apply some prefiltering
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // No active request handles
        return;
    }

    profilingPstream::beginTiming();

    // On success: sets each request to MPI_REQUEST_NULL
    if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitall returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    // Everything handled, reset all to MPI_REQUEST_NULL
    requests = UPstream::Request(MPI_REQUEST_NULL);
}


bool Foam::UPstream::waitAnyRequest(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return false;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitAnyRequest : starting wait for any of "
            << count << " requests starting at " << pos << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    int index = MPI_UNDEFINED;
    if (MPI_Waitany(count, waitRequests, &index, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitany returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (index == MPI_UNDEFINED)
    {
        // No active request handles
        return false;
    }

    return true;
}


bool Foam::UPstream::waitSomeRequests
(
    const label pos,
    label len,
    DynamicList<int>* indices
)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        if (indices) indices->clear();
        return false;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream:waitSomeRequest : starting wait for some of "
            << count << " requests starting at " << pos << endl;
    }


    // Local temporary storage, or return via calling parameter
    List<int> tmpIndices;
    if (indices)
    {
        indices->resize_nocopy(count);
    }
    else
    {
        tmpIndices.resize(count);
    }

    profilingPstream::beginTiming();

    // On success: sets non-blocking requests to MPI_REQUEST_NULL
    int outcount = 0;
    if
    (
        MPI_Waitsome
        (
            count,
            waitRequests,
           &outcount,
            (indices ? indices->data() : tmpIndices.data()),
            MPI_STATUSES_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Waitsome returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (outcount == MPI_UNDEFINED || outcount < 1)
    {
        // No active request handles
        if (indices) indices->clear();
        return false;
    }

    if (indices)
    {
        indices->resize(outcount);
    }

    return true;
}


bool Foam::UPstream::waitSomeRequests
(
    UList<UPstream::Request>& requests,
    DynamicList<int>* indices
)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        if (indices) indices->clear();
        return false;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        waitRequests[count] = PstreamUtils::Cast::to_mpi(req);
        ++count;
    }

    // Local temporary storage, or return via calling parameter
    List<int> tmpIndices;
    if (indices)
    {
        indices->resize_nocopy(count);
    }
    else
    {
        tmpIndices.resize(count);
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream:waitSomeRequest : starting wait for some of "
            << requests.size() << " requests" << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets non-blocking requests to MPI_REQUEST_NULL
    int outcount = 0;
    if
    (
        MPI_Waitsome
        (
            count,
            waitRequests,
           &outcount,
            (indices ? indices->data() : tmpIndices.data()),
            MPI_STATUSES_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Waitsome returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (outcount == MPI_UNDEFINED || outcount < 1)
    {
        // No active request handles
        if (indices) indices->clear();

        // Everything handled or inactive, reset all to MPI_REQUEST_NULL
        requests = UPstream::Request(MPI_REQUEST_NULL);
        return false;
    }

    if (indices)
    {
        indices->resize(outcount);
    }

    // Transcribe MPI_Request back into UPstream::Request
    // - do in reverse order - see note in finishedRequests()
    {
        for (label i = requests.size()-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }
    }

    return true;
}


Foam::label Foam::UPstream::waitAnyRequest(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return -1;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    // Transcribe UPstream::Request into MPI_Request
    // - do not change locations within the list since these are relevant
    //   for the return index.
    for (auto& req : requests)
    {
        waitRequests[count] = PstreamUtils::Cast::to_mpi(req);
        ++count;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    int index = MPI_UNDEFINED;
    if (MPI_Waitany(count, waitRequests, &index, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitany returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (index == MPI_UNDEFINED)
    {
        index = -1;  // No outstanding requests
    }

    // Transcribe MPI_Request back into UPstream::Request
    // - do in reverse order - see note in finishedRequests()
    {
        for (label i = count-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }

        // Trailing portion
        for (label i = count; i < requests.size(); ++i)
        {
            requests[i] = UPstream::Request(MPI_REQUEST_NULL);
        }
    }

    return index;
}


// FUTURE?
//
/// void Foam::UPstream::waitRequests
/// (
///     UPstream::Request& req0,
///     UPstream::Request& req1
/// )
/// {
///     // No-op for non-parallel
///     if (!UPstream::parRun())
///     {
///         return;
///     }
///
///     int count = 0;
///     MPI_Request waitRequests[2];
///
///     waitRequests[count] = PstreamUtils::Cast::to_mpi(req0);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         ++count;
///     }
///
///     waitRequests[count] = PstreamUtils::Cast::to_mpi(req1);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         ++count;
///     }
///
///     // Flag in advance as being handled
///     req0 = UPstream::Request(MPI_REQUEST_NULL);
///     req1 = UPstream::Request(MPI_REQUEST_NULL);
///
///     if (!count)
///     {
///         return;
///     }
///
///     profilingPstream::beginTiming();
///
///     // On success: sets each request to MPI_REQUEST_NULL
///     if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
///     {
///         FatalErrorInFunction
///             << "MPI_Waitall returned with error"
///             << Foam::abort(FatalError);
///     }
///
///     profilingPstream::addWaitTime();
/// }


void Foam::UPstream::waitRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return;
    }

    auto& request = PstreamGlobals::outstandingRequests_[i];

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return;
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequest : starting wait for request:"
            << i << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    if (MPI_Wait(&request, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequest : finished wait for request:"
            << i << endl;
    }
}


void Foam::UPstream::waitRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    MPI_Request request = PstreamUtils::Cast::to_mpi(req);

    // No-op for null request
    if (MPI_REQUEST_NULL == request)
    {
        return;
    }

    profilingPstream::beginTiming();

    if (MPI_Wait(&request, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
}


bool Foam::UPstream::finishedRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return true;
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::finishedRequest : check request:"
            << i << endl;
    }

    auto& request = PstreamGlobals::outstandingRequests_[i];

    // Fast-path (no-op) for null request
    if (MPI_REQUEST_NULL == request)
    {
        return true;
    }

    // On success: sets request to MPI_REQUEST_NULL
    int flag = 0;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);

    return flag != 0;
}


bool Foam::UPstream::finishedRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return true;
    }

    MPI_Request request = PstreamUtils::Cast::to_mpi(req);

    // Fast-path (no-op) for null request
    if (MPI_REQUEST_NULL == request)
    {
        return true;
    }

    int flag = 0;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);

    if (flag)
    {
        // Success: now inactive
        req = UPstream::Request(MPI_REQUEST_NULL);
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequests(const label pos, label len)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return true;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    if (UPstream::debug)
    {
        Perr<< "UPstream::finishedRequests : check " << count
            << " requests starting at " << pos << endl;
    }

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    int flag = 1;

    if (count == 1)
    {
        // Fast-path (no-op) for single null request
        if (MPI_REQUEST_NULL == *waitRequests)
        {
            return true;
        }

        // On success: sets request to MPI_REQUEST_NULL
        MPI_Test(waitRequests, &flag, MPI_STATUS_IGNORE);
    }
    else if (count > 1)
    {
        // On success: sets each request to MPI_REQUEST_NULL
        // On failure: no request is modified
        MPI_Testall(count, waitRequests, &flag, MPI_STATUSES_IGNORE);
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return true;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        if (MPI_REQUEST_NULL != request)  // Apply some prefiltering
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // No active handles
        return true;
    }

    // On success: sets each request to MPI_REQUEST_NULL
    // On failure: no request is modified
    int flag = 0;
    MPI_Testall(count, waitRequests, &flag, MPI_STATUSES_IGNORE);

    if (flag)
    {
        // Success: reset all requests to MPI_REQUEST_NULL
        requests = UPstream::Request(MPI_REQUEST_NULL);
    }
    else
    {
        // Not all done. Recover wrapped representation but in reverse order
        // since sizeof(MPI_Request) can be smaller than
        // sizeof(UPstream::Request::value_type)
        // eg, mpich has MPI_Request as 'int'
        //
        // This is uglier that we'd like, but much better than allocating
        // and freeing a scratch buffer each time we query things.

        for (label i = count-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }

        // Trailing portion
        for (label i = count; i < requests.size(); ++i)
        {
            requests[i] = UPstream::Request(MPI_REQUEST_NULL);
        }
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequestPair(label& req0, label& req1)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        req0 = -1;
        req1 = -1;
        return true;
    }

    bool anyActive = false;
    MPI_Request waitRequests[2];

    // No-op for out-of-range (eg, placeholder indices)

    if (req0 >= 0 && req0 < PstreamGlobals::outstandingRequests_.size())
    {
        waitRequests[0] = PstreamGlobals::outstandingRequests_[req0];
    }
    else
    {
        waitRequests[0] = MPI_REQUEST_NULL;
    }

    if (req1 >= 0 && req1 < PstreamGlobals::outstandingRequests_.size())
    {
        waitRequests[1] = PstreamGlobals::outstandingRequests_[req1];
    }
    else
    {
        waitRequests[1] = MPI_REQUEST_NULL;
    }

    if (MPI_REQUEST_NULL != waitRequests[0])  // An active handle
    {
        anyActive = true;
    }
    else
    {
        req0 = -1;
    }

    if (MPI_REQUEST_NULL != waitRequests[1])  // An active handle
    {
        anyActive = true;
    }
    else
    {
        req1 = -1;
    }

    if (!anyActive)
    {
        // No active handles
        return true;
    }

    profilingPstream::beginTiming();

    // On success: sets each request to MPI_REQUEST_NULL
    int indices[2];
    int outcount = 0;
    if
    (
        MPI_Testsome
        (
            2,
            waitRequests,
           &outcount,
            indices,
            MPI_STATUSES_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Testsome returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (outcount == MPI_UNDEFINED)
    {
        // No active request handles.
        // Slight pedantic, but copy back requests in case they were altered

        if (req0 >= 0)
        {
            PstreamGlobals::outstandingRequests_[req0] = waitRequests[0];
        }

        if (req1 >= 0)
        {
            PstreamGlobals::outstandingRequests_[req1] = waitRequests[1];
        }

        // Flag indices as 'done'
        req0 = -1;
        req1 = -1;
        return true;
    }

    // Copy back requests to their 'stack' locations
    for (int i = 0; i < outcount; ++i)
    {
        const int idx = indices[i];

        if (idx == 0)
        {
            if (req0 >= 0)
            {
                PstreamGlobals::outstandingRequests_[req0] = waitRequests[0];
                req0 = -1;
            }
        }
        if (idx == 1)
        {
            if (req1 >= 0)
            {
                PstreamGlobals::outstandingRequests_[req1] = waitRequests[1];
                req1 = -1;
            }
        }
    }

    return (outcount > 0);
}


void Foam::UPstream::waitRequestPair(label& req0, label& req1)
{
    // No-op for non-parallel. Flag indices as 'done'
    if (!UPstream::parRun())
    {
        req0 = -1;
        req1 = -1;
        return;
    }

    int count = 0;
    MPI_Request waitRequests[2];

    // No-op for out-of-range (eg, placeholder indices)
    // Prefilter inactive handles

    if (req0 >= 0 && req0 < PstreamGlobals::outstandingRequests_.size())
    {
        waitRequests[count] = PstreamGlobals::outstandingRequests_[req0];
        PstreamGlobals::outstandingRequests_[req0] = MPI_REQUEST_NULL;

        if (MPI_REQUEST_NULL != waitRequests[count])  // An active handle
        {
            ++count;
        }
    }

    if (req1 >= 0 && req1 < PstreamGlobals::outstandingRequests_.size())
    {
        waitRequests[count] = PstreamGlobals::outstandingRequests_[req1];
        PstreamGlobals::outstandingRequests_[req1] = MPI_REQUEST_NULL;

        if (MPI_REQUEST_NULL != waitRequests[count])  // An active handle
        {
            ++count;
        }
    }

    // Flag in advance as being handled
    req0 = -1;
    req1 = -1;

    if (!count)
    {
        // No active handles
        return;
    }

    profilingPstream::beginTiming();

    // On success: sets each request to MPI_REQUEST_NULL
    if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitall returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();
}


// ************************************************************************* //
