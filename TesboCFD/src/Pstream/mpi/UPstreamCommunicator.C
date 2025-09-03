/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Communicator::Communicator() noexcept
:
    UPstream::Communicator(MPI_COMM_NULL)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::UPstream::Communicator
Foam::UPstream::Communicator::lookup(const label comm)
{
    if (comm < 0 || comm >= PstreamGlobals::MPICommunicators_.size())
    {
        WarningInFunction
            << "Illegal communicator " << comm << nl
            << "Should be within range [0,"
            << PstreamGlobals::MPICommunicators_.size()
            << ')' << endl;

        return UPstream::Communicator(MPI_COMM_NULL);
    }

    return UPstream::Communicator(PstreamGlobals::MPICommunicators_[comm]);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Communicator::good() const noexcept
{
    return MPI_COMM_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Communicator::reset() noexcept
{
    *this = UPstream::Communicator(MPI_COMM_NULL);
}


// ************************************************************************* //
