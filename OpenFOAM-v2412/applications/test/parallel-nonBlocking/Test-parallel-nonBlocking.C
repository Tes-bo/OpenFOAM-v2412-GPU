/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Test-parallel-nonBlocking

Description
    Test for various non-blocking parallel routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "mapDistribute.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "Tuple2.H"
#include "PstreamBuffers.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList args(argc, argv);

    // Test PstreamBuffers
    // ~~~~~~~~~~~~~~~~~~~
    if (false)
    {
        Perr<< "\nStarting transfers\n" << endl;

        vector data
        (
            Pstream::myProcNo(),
            Pstream::myProcNo(),
            Pstream::myProcNo()
        );

        PstreamBuffers pBufs;

        if (!Pstream::master())
        {
            Perr<< "sending to master" << endl;
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster << data;
        }

        pBufs.finishedGathers();

        // Consume
        DynamicList<vector> allData;
        if (Pstream::master())
        {
            // Collect my own data
            allData.push_back(data);

            for (const int proci : Pstream::subProcs())
            {
                Perr << "master receiving from " << proci << endl;
                UIPstream fromProc(proci, pBufs);
                allData.push_back(vector(fromProc));
            }
        }


        // Send allData back
        pBufs.clear();
        if (Pstream::master())
        {
            for (const int proci : Pstream::subProcs())
            {
                Perr << "master sending to " << proci << endl;
                UOPstream toProc(proci, pBufs);
                toProc << allData;
            }
        }

        pBufs.finishedScatters();

        // Consume
        if (!Pstream::master())
        {
            Perr<< "receive from master" << endl;
            UIPstream fromMaster(Pstream::masterNo(), pBufs);
            fromMaster >> allData;
            Perr<< allData << endl;
        }
    }


    // Test non-blocking reductions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalar data1 = 1.0;
    label request1 = -1;
    {
        Foam::reduce
        (
            data1,
            sumOp<scalar>(),
            UPstream::msgType(),
            UPstream::worldComm,
            request1
        );
    }

    scalar data2 = 0.1;
    UPstream::Request request2;
    {
        Foam::reduce
        (
            data2,
            sumOp<scalar>(),
            UPstream::msgType(),
            UPstream::worldComm,
            request2
        );
    }


    // Do a non-blocking send inbetween
    {
        PstreamBuffers pBufs;

        for (const int proci : Pstream::allProcs())
        {
            UOPstream toProc(proci, pBufs);
            toProc << Pstream::myProcNo();
        }

        // Start sending and receiving and block
        pBufs.finishedSends();

        // Consume
        for (const int proci : Pstream::allProcs())
        {
            UIPstream fromProc(proci, pBufs);
            label data;
            fromProc >> data;

            if (data != proci)
            {
                FatalErrorInFunction
                    << "From processor " << proci << " received " << data
                    << " but expected " << proci
                    << exit(FatalError);
            }
        }
    }


    if (request1 != -1)
    {
        Pout<< "Waiting for non-blocking reduce with request "
            << request1 << endl;
        UPstream::waitRequest(request1);
    }
    Info<< "Reduced data1:" << data1 << endl;

    if (request2.good())
    {
        Pout<< "Waiting for non-blocking reduce with request "
            << Foam::name(request2.pointer()) << endl;
        UPstream::waitRequest(request2);
    }
    Info<< "Reduced data2:" << data2 << endl;


    // Clear all outstanding requests
    UPstream::resetRequests(0);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
