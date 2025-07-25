/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::broadcast
(
    char* buf,
    const std::streamsize bufSize,
    const label comm,
    const int rootProcNo
)
{
    if (!UPstream::is_parallel(comm))
    {
        // Nothing to do - ignore
        return true;
    }

    //Needed?  PstreamGlobals::checkCommunicator(comm, rootProcNo);

    if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
    {
        Perr<< "UPstream::broadcast : root:" << rootProcNo
            << " comm:" << comm
            << " size:" << label(bufSize)
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Perr);
    }
    else if (UPstream::debug)
    {
        Perr<< "UPstream::broadcast : root:" << rootProcNo
            << " comm:" << comm
            << " size:" << label(bufSize)
            << Foam::endl;
    }

    profilingPstream::beginTiming();

    const int returnCode = MPI_Bcast
    (
        buf,
        bufSize,
        MPI_BYTE,
        rootProcNo,
        PstreamGlobals::MPICommunicators_[comm]
    );

    profilingPstream::addBroadcastTime();

    return (returnCode == MPI_SUCCESS);
}


// ************************************************************************* //
