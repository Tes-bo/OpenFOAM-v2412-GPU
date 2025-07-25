/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    Test-globalIndex1

Description
    Simple tests for the globalIndex class.

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "globalMeshData.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    // Global numbering of cells (proc0 elements first, then proc1, etc.)
    const globalIndex globalNumbering(mesh.nCells());

    Pout<< "local-offset: " << globalIndex::calcOffset(mesh.nCells()) << nl;
    Pout<< "local-range: " << globalIndex::calcRange(mesh.nCells()) << nl;

    Info<< "cells from:" << globalNumbering.begin_value()
        << " to:" << globalNumbering.end_value()
        << " span:" << globalNumbering.span() << nl;

    Info<< "front: " << globalNumbering.front()
        << " back: " << globalNumbering.back() << nl;

    if (globalNumbering.localSize() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Problem." << abort(FatalError);
    }


    if (!UPstream::parRun())
    {
        WarningInFunction
            << "globalIndex class is only useful in parallel code."
            << endl;
    }

    // convert from local to global and back.
    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        // to global index
        label globalCelli = globalNumbering.toGlobal(celli);

        // and back
        label proci = globalNumbering.whichProcID(globalCelli);
        label localCelli = globalNumbering.toLocal(globalCelli);

        if (proci != UPstream::myProcNo() || localCelli != celli)
        {
            FatalErrorInFunction
                << "Problem. celli:" << celli << " localCelli:" << localCelli
                << " proci:" << proci << abort(FatalError);
        }

        if (!globalNumbering.isLocal(globalCelli))
        {
            FatalErrorInFunction
                << "Problem. celli:" << celli << " globalCelli:" << globalCelli
                << " not local" << abort(FatalError);
        }
    }


    Info<< "Max number of cells: " << globalNumbering.maxSize() << nl;
    Pout<< "nCells: "
        << globalNumbering.localSize() << " Max off-processor: "
        << globalNumbering.maxNonLocalSize() << nl;

    // Try whichProcID on a few borderline cases.

    if (mesh.nCells() < 1)
    {
        FatalErrorInFunction
            << "Test needs to be run on a case with at least one"
            << " cell per processor." << abort(FatalError);
    }

    if (UPstream::myProcNo() > 0)
    {
        // We already checked that toGlobal(0) maps back correctly to myProcNo
        // so now check that the index one before maps to the previous processor
        label prevProcCelli = globalNumbering.toGlobal(0)-1;
        label proci = globalNumbering.whichProcID(prevProcCelli);

        if (proci != UPstream::myProcNo()-1)
        {
            FatalErrorInFunction
                << "Problem. global:" << prevProcCelli
                << " expected on processor:" << UPstream::myProcNo()-1
                << " but is calculated to be on proci:" << proci
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(prevProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << prevProcCelli
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(proci, prevProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << prevProcCelli
                << " not calculated as local on processor:" << proci
                << abort(FatalError);
        }
    }


    if (UPstream::myProcNo() < UPstream::nProcs()-1)
    {
        label nextProcCelli = globalNumbering.toGlobal(mesh.nCells()-1)+1;
        label proci = globalNumbering.whichProcID(nextProcCelli);

        if (proci != UPstream::myProcNo()+1)
        {
            FatalErrorInFunction
                << "Problem. global:" << nextProcCelli
                << " expected on processor:" << UPstream::myProcNo()+1
                << " but is calculated to be on proci:" << proci
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(nextProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << nextProcCelli
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(proci, nextProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << nextProcCelli
                << " not calculated as local on processor:" << proci
                << abort(FatalError);
        }
    }



    // Get a few cell indices
    const label nTotalCells = globalNumbering.totalSize();

    Random rndGen(UPstream::myProcNo());
    DynamicList<label> globalIDs;
    for (label i = 0; i < 100; i++)
    {
        globalIDs.push_back(rndGen.position(label(0), nTotalCells-1));
    }

    // Get the cell centres at those cell indices
    List<point> ccs;
    globalNumbering.get
    (
        ccs,
        globalIDs,
        [&mesh](List<point>& ccs, const labelUList& localIds)
        {
            ccs = UIndirectList<point>(mesh.cellCentres(), localIds);
        }
    );

    // Do the brute-force method as well : collect all cell centres on all
    // processors

    Info<< "Gathered/scattered cell centres:" << endl;

    labelPair inOut;
    pointField allCcs(globalNumbering.gather(mesh.cellCentres()));
    inOut[0] = allCcs.size();
    Pstream::broadcastList(allCcs);
    inOut[1] = allCcs.size();
    Pout<< "    " << inOut << endl;

    // Compare
    forAll(ccs, i)
    {
        const point& cc = ccs[i];
        const point& allCc = allCcs[globalIDs[i]];
        if (cc != allCc)
        {
            FatalErrorInFunction << "Problem cc:" << cc
                << " allCc:" << allCc << exit(FatalError);
        }
    }

    Info<< "Gather indirect list of boundary points (patch = 0)\n";
    {
        const polyPatch& pp = mesh.boundaryMesh()[0];

        // Map mesh point index to local (compact) point index
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;

        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        Info<< "local-sizes: " << globalPointsPtr().sizes() << nl;

        pointField patchPoints
        (
            globalPointsPtr().gather
            (
                UIndirectList<point>(mesh.points(), uniqueMeshPointLabels)
            )
        );

        Info<< "gathered point field = " << patchPoints.size() << " points\n";
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
