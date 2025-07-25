/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

Application
    Test-parallel-external-init

Description
    Simulate starting MPI outside of OpenFOAM

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Pstream.H"
#include <iostream>

#include "openfoam_mpi.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool startMPI()
{
    enum whichComm : int { worldComm = 0, selfComm, nullComm };

    int nprocs[3];
    int rank[3];

    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs[worldComm]);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank[worldComm]);

    const bool isMaster = (rank[worldComm] == 0);
    const string prefix = '[' + Foam::name(rank[worldComm]) + "] ";

    if (isMaster && nprocs[worldComm])
    {
        std::cout
            << nl << "Using MPI with "
            << nprocs[worldComm] << " procs" << nl
            << nl;
    }

    MPI_Comm worldMpiComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &worldMpiComm);
    MPI_Comm_free(&worldMpiComm);

    MPI_Comm_size(MPI_COMM_SELF, &nprocs[selfComm]);
    MPI_Comm_rank(MPI_COMM_SELF, &rank[selfComm]);

    // if (nprocs && isMaster)
    {
        std::cout
            << prefix
            << "Self: " << rank[selfComm]
            << " from " << nprocs[selfComm] << " procs" << nl;
    }

    if (isMaster)
    {
        std::cout
            << "MPI_COMM_NULL:  " << MPI_COMM_NULL << nl
            << "MPI_COMM_SELF:  " << MPI_COMM_SELF << nl
            << "MPI_COMM_WORLD: " << MPI_COMM_WORLD << nl;
    }

    return true;
}


bool stopMPI()
{
    Info<< nl << "Stopping MPI" << nl << nl;

    MPI_Finalize();

    return true;
}


string message()
{
    return
    (
        "rank " + Foam::name(Pstream::myProcNo())
      + " / " + Foam::name(Pstream::nProcs()) + "\n"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Set UPstream::debug level");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    startMPI();

    #include "setRootCase.H"

    Pout<< message().c_str();

    stopMPI();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
