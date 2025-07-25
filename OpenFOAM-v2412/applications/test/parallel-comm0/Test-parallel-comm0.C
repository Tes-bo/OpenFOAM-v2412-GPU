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

Application
    Test-parallel-comm0

Description
    Very basic checks on standard communicators

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pair.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "bitSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void printInfo(const label comm)
{
    Info<< "comm:" << comm
        << " nprocs:" << UPstream::nProcs(comm)
        << " all:" << UPstream::allProcs(comm)
        << " sub:" << UPstream::subProcs(comm) << nl;


    if (UPstream::commSelf() == comm)
    {
        Pout<< "self all:" << UPstream::allProcs(comm)
            << " sub:" << UPstream::subProcs(comm) << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Set UPstream::debug level");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    #include "setRootCase.H"

    Info<< nl
        << "nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)" << nl;

    Info<< "comm-world : ";
    printInfo(UPstream::commWorld());

    Info<< "comm-self  : ";
    printInfo(UPstream::commSelf());

    Info<< nl;

    // Reductions (using MPI intrinsics)
    {
        const label myRank = UPstream::myProcNo(UPstream::commWorld());
        const label nProcs = UPstream::nProcs(UPstream::commWorld());

        label val = myRank;

        label worldVal = returnReduce
        (
            val,
            sumOp<label>(),
            UPstream::msgType(),
            UPstream::commWorld()
        );

        label selfVal = returnReduce
        (
            val,
            sumOp<label>(),
            UPstream::msgType(),
            UPstream::commSelf()
        );

        Pout<< "value " << val
            << " (world) reduced " << worldVal
            << " (self) reduced " << selfVal << nl;

        // Identical size on all procs
        bitSet procUsed(nProcs);

        if ((myRank % 4) == 0)
        {
            procUsed.set(myRank);
        }

        Pout<< "local  procUsed " << procUsed << nl;
        reduce
        (
            procUsed.data(),
            procUsed.size_data(),
            bitOrOp<unsigned int>()
        );
        Pout<< "reduce procUsed " << procUsed << nl;

        // Identical size on all procs
        // encode as 0:empty, 1:uniform, 2:nonuniform, 3:mixed
        PackedList<2> uniformity(10);

        if ((myRank % 2) == 0)
        {
            // Every second is uniform
            uniformity.set(2, 1);
            uniformity.set(4, 1);
            uniformity.set(6, 1);
            uniformity.set(8, 1);
        }
        else if ((myRank % 3) == 0)
        {
            // Every third is nonuniform
            uniformity.set(3, 2);
            uniformity.set(6, 2);
            uniformity.set(9, 2);
        }

        Pout<< "local  uniform " << uniformity << nl;
        reduce
        (
            uniformity.data(),
            uniformity.size_data(),
            bitOrOp<unsigned int>()
        );
        Pout<< "reduce uniform " << uniformity << nl;
    }

    // Reductions (not using MPI intrinsics)
    {
        Pair<label> val
        (
            Pstream::myProcNo(UPstream::commWorld()),
            Pstream::myProcNo(UPstream::commWorld())
        );

        Pair<label> worldVal = val;

        Pstream::combineReduce
        (
            worldVal,
            minFirstEqOp<label>(),
            UPstream::msgType(),
            UPstream::commWorld()
        );

        Pair<label> selfVal = val;

        Pstream::combineReduce
        (
            worldVal,
            minFirstEqOp<label>(),
            UPstream::msgType(),
            UPstream::commSelf()
        );

        Pout<< "value " << val
            << " (world) reduced " << worldVal
            << " (self) reduced " << selfVal << nl;
    }

    Pout<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
