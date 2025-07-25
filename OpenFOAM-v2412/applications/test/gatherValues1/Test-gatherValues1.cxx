/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Test-gatherValues1

Description
    Test list gather functionality

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    const labelList localValues
    (
        identity(2 *(UPstream::myProcNo()+1), -5*UPstream::myProcNo())
    );

    // Test resize
    {
        globalIndex globIdx(localValues.size());

        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;

        globIdx.setLocalSize(4, 0);
        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;
        globIdx.setLocalSize(3, 0);
        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;
    }


    // This now compiles (and works)
    // - reverts to normal gather for non-contiguous
    {
        const wordList sendData({"hello", "world"});

        // One-sided sizing!  master only
        const globalIndex allProcAddr
        (
            globalIndex::gatherOnly{},
            sendData.size()
        );

        Pout<< "listGather sizes: " << flatOutput(allProcAddr.sizes()) << nl;

        // Collect all values
        wordList allValues
        (
            allProcAddr.mpiGather(sendData)
        );

        Pout<< "all-data: " << allValues << endl;
    }

    // Gather all values
    {
        const auto& sendData = localValues;

        // One-sided sizing!  master only
        const globalIndex allProcAddr
        (
            globalIndex::gatherOnly{},
            sendData.size()
        );

        Pout<< "listGather sizes: " << flatOutput(allProcAddr.sizes()) << nl;

        // Collect all values
        labelList allValues
        (
            allProcAddr.mpiGather(sendData)
        );

        Pout<< "all-data: " << allValues << endl;
    }

    {
        const labelList::subList& sendData =
        (
            UPstream::master()
          ? SubList<label>(localValues, 0)  // exclude
          : SubList<label>(localValues)
        );

        const labelList sendSizes
        (
            UPstream::listGatherValues<label>(sendData.size())
        );

        const label sendSize
        (
            UPstream::listScatterValues<label>(sendSizes)
        );

        const globalIndex subProcAddr(globalIndex::gatherNone{}, sendSizes);

        Pout<< "listGather "
            << localValues.size() << " = " << flatOutput(sendSizes)
            << " offsets " << flatOutput(subProcAddr.offsets())
            << nl;

        label newLocalValue = 5 + UPstream::listScatterValues(sendSizes);

        Pout<< "listScattered: " << newLocalValue << nl;

        // Can also scatter a longer list
        Pout<< "listScatter off "
            << UPstream::listScatterValues(subProcAddr.offsets()) << nl;


        Pout<< endl << "local list [" << UPstream::myProcNo() << "] "
            << flatOutput(localValues) << nl;


        Pout<< endl << "local send [" << UPstream::myProcNo() << "] "
            << sendSize << nl;


        // Collect all off-processor values
        labelList allValues
        (
            subProcAddr.mpiGather(sendData)
        );

        Pout<< "off-proc: " << allValues << endl;

        if (UPstream::master())
        {
            Info<< "master: " << flatOutput(localValues) << nl;

            label proci = 0;
            for (const labelRange& range : subProcAddr)
            {
                Info<< proci << ": " << flatOutput(allValues.slice(range))
                    << nl;
                ++proci;
            }

            Info<< nl << "verify ranges" << nl;

            {
                globalIndex glob;
                Info<< "empty:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
            {
                globalIndex glob(labelList(Foam::one{}, 0));
                Info<< "degenerate:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
            {
                globalIndex glob
                (
                    globalIndex::gatherNone{},
                    labelList(Foam::one{}, 0)
                );
                Info<< "single:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
        }
    }

    // Non-contiguous gather - use Pstream, not UPstream!
    {
        typedef std::pair<label,vector> valueType;

        valueType sendData(UPstream::myProcNo(), vector::one);

        List<valueType> countValues
        (
            Pstream::listGatherValues(sendData)
        );

        Pout<< "listGather: " << flatOutput(countValues) << nl;
    }

    // Non-contiguous scatter - use Pstream, not UPstream!
    {
        List<fileName> allValues;

        if (UPstream::master())
        {
            allValues.resize(UPstream::nProcs());
            forAll(allValues, proci)
            {
                allValues[proci] = "processor" + Foam::name(proci);
            }
        }

        fileName procName = Pstream::listScatterValues(allValues);

        Pout<< "listScatter: " << procName << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
