/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Test-checkDecomposePar

Group
    grpParallelUtilities

Description
    Check decomposition from kaffpa (KaHIP) output.
    foamToMetisGraph was likely used for producing the kaffpa input.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "polyMesh.H"
#include "cpuTime.H"
#include "IFstream.H"
#include "regionProperties.H"
#include "globalMeshData.H"
#include "decompositionInformation.H"
#include "decompositionModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addNote
    (
        "Check decomposition from kaffpa (KaHIP) output"
    );

    argList::noParallel();
    argList::noBanner();

    #include "addAllRegionOptions.H"

    argList::addVerboseOption
    (
        "more information about decomposition"
    );

    argList::addArgument("kaffpa-output-file");

    #include "setRootCase.H"

    const auto decompFile = args.get<fileName>(1);

    // Set time from database
    #include "createTime.H"

    // Allow override of time from specified time options, or no-op
    timeSelector::setTimeIfPresent(runTime, args);

    // Allow override of decomposeParDict location
    const fileName decompDictFile =
        args.getOrDefault<fileName>("decomposeParDict", "");

    // Get region names
    #include "getAllRegionOptions.H"

    labelList cellToProc;

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        // const word& regionDir = polyMesh::regionName(regionName);

        Info<< "\n\nDecomposing mesh " << regionName << nl << endl;
        Info<< "Create mesh..." << flush;

        polyMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        Info<< " nCells = " << mesh.nCells() << endl;

        // Expected format is a simple ASCII list
        cellToProc.resize(mesh.nCells());
        {
            IFstream is(decompFile);

            forAll(cellToProc, celli)
            {
                cellToProc[celli] = readLabel(is);
            }
        }

        const label nDomains = max(cellToProc) + 1;

        // Local mesh connectivity
        CompactListList<label> cellCells;
        globalMeshData::calcCellCells(mesh, cellCells);

        decompositionInformation info
        (
            cellCells,
            cellToProc,
            nDomains
        );

        if (args.verbose())
        {
            info.printDetails(Info);
            Info<< nl;
        }
        info.printSummary(Info);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
