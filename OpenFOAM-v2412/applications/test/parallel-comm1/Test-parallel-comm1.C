/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Test-parallel-communicators

Description
    Checks communication using user-defined communicators

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar sumReduce
(
    const label comm,
    const scalar localValue
)
{
    scalar sum = 0;
    if (UPstream::parRun())
    {
        if (UPstream::master(comm))
        {
            // Add master value and all sub-procs
            sum = localValue;

            for (const int proci : UPstream::subProcs(comm))
            {
                scalar procValue;
                UIPstream::read
                (
                    UPstream::commsTypes::buffered,
                    proci,
                    reinterpret_cast<char*>(&procValue),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );

                sum += procValue;
            }

            // Send back

            for (const int proci : UPstream::subProcs(comm))
            {
                UOPstream::write
                (
                    UPstream::commsTypes::buffered,
                    proci,
                    reinterpret_cast<const char*>(&sum),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }
        }
        else
        {
            {
                UOPstream::write
                (
                    UPstream::commsTypes::buffered,
                    UPstream::masterNo(),
                    reinterpret_cast<const char*>(&localValue),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }

            {
                UIPstream::read
                (
                    UPstream::commsTypes::buffered,
                    UPstream::masterNo(),
                    reinterpret_cast<char*>(&sum),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }
        }
    }
    return sum;
}


int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();
    argList::addOption("repeat", "count");

    #include "setRootCase.H"

    const label repeat = args.getOrDefault<label>("repeat", 1);
    const int optVerbose = args.verbose();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Allocate a communicator
    label n = Pstream::nProcs(UPstream::worldComm);

    DynamicList<label> bottom(n/2);
    DynamicList<label> top(n/2);

    for (label i = 0; i < n/2; i++)
    {
        bottom.append(i);
    }
    for (label i = n/2; i < n; i++)
    {
        top.append(i);
    }


    Info<< "repeating: " << repeat << " times" << endl;

    if (optVerbose || repeat == 1)
    {
        //Pout<< "bottom:" << bottom << endl;
        Pout<< "top             :" << top << endl;
    }

    for (label count = 0; count < repeat; ++count)
    {
        label comm = UPstream::allocateCommunicator(UPstream::worldComm, top);

        scalar localValue = 111*UPstream::myProcNo(UPstream::worldComm);

        if (optVerbose || (repeat == 1 && count == 0))
        {
            Pout<< "localValue      :" << localValue << endl;
            Pout<< "allocated comm  :" << comm
                << " proci :" << UPstream::myProcNo(comm) << endl;
        }

        if (Pstream::myProcNo(comm) != -1)
        {
            scalar sum = returnReduce
            (
                localValue,
                sumOp<scalar>(),
                UPstream::msgType(),
                comm
            );

            if (optVerbose || (repeat == 1 && count == 0))
            {
                Pout<< "sum             :" << sum << endl;
            }
        }

        UPstream::freeCommunicator(comm);
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
