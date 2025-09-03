/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for MPI.

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamGlobals_H
#define Foam_PstreamGlobals_H

#include "DynamicList.H"
#include "UPstream.H"  // For UPstream::Request
#include "openfoam_mpi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamGlobals
{

// Track if MPI_Comm_free is needed for communicator index in MPICommunicators_
extern DynamicList<bool> pendingMPIFree_;

// Current communicators, which may be allocated or predefined
// (eg, MPI_COMM_SELF, MPI_COMM_WORLD)
extern DynamicList<MPI_Comm> MPICommunicators_;

//- Outstanding non-blocking operations.
extern DynamicList<MPI_Request> outstandingRequests_;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Fatal if comm is outside the allocated range
void checkCommunicator(const label comm, const label toProcNo);


//- Reset UPstream::Request to null and/or the index of the outstanding
//- request to -1.
//  Does not affect the stack of outstanding requests.
inline void reset_request
(
    UPstream::Request* requestPtr,
    label* requestIdx = nullptr
)
{
    if (requestPtr) *requestPtr = UPstream::Request(MPI_REQUEST_NULL);
    if (requestIdx) *requestIdx = -1;
}


//- Transcribe MPI_Request to UPstream::Request
//- (does not affect the stack of outstanding requests)
//- or else push onto list of outstanding requests
//- and (optionally) record its location
inline void push_request
(
    MPI_Request request,
    UPstream::Request* requestPtr = nullptr,
    label* requestIdx = nullptr
)
{
    if (requestPtr)
    {
        // Transcribe as UPstream::Request
        *requestPtr = UPstream::Request(request);

        // Not on stack of outstanding requests
        if (requestIdx) *requestIdx = -1;
    }
    else
    {
        if (requestIdx)
        {
            // Its index into outstanding requests
            *requestIdx = PstreamGlobals::outstandingRequests_.size();
        }

        PstreamGlobals::outstandingRequests_.push_back(request);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamGlobals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
