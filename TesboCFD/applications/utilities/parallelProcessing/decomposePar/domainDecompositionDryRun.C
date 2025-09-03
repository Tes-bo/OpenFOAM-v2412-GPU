/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "domainDecompositionDryRun.H"
#include "volFields.H"
#include "globalMeshData.H"
#include "decompositionModel.H"
#include "decompositionInformation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecompositionDryRun::domainDecompositionDryRun
(
    const IOobject& io,
    const fileName& decompDictFile,
    const label nDomains,
    const word& methodName
)
:
    mesh_(io),
    decompDictFile_(decompDictFile),
    nDomainsOverride_(nDomains),
    methodNameOverride_(methodName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::domainDecompositionDryRun::execute
(
    const bool writeCellDist,
    const bool verbose
)
{
    cpuTime decompositionTime;

    Info<< "\nCalculating distribution of cells. nCells = "
        << mesh_.nCells() << endl;

    const decompositionModel& model = decompositionModel::New
    (
        mesh_,
        decompDictFile_
    );

    // Allow overrides for testing

    dictionary& modelDict = const_cast<decompositionModel&>(model);

    if (nDomainsOverride_ > 0)
    {
        modelDict.add
        (
            word("numberOfSubdomains"),
            nDomainsOverride_,
            true
        );
    }

    if (!methodNameOverride_.empty())
    {
        modelDict.add
        (
            word("method"),
            methodNameOverride_,
            true
        );
    }

    scalarField cellWeights;
    word weightName;
    if (model.readIfPresent("weightField", weightName))
    {
        volScalarField weights
        (
            IOobject
            (
                weightName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
        cellWeights = weights.primitiveField();
    }

    decompositionMethod& method = model.decomposer();

    // Local mesh connectivity
    CompactListList<label> cellCells;
    globalMeshData::calcCellCells(mesh_, cellCells);

    labelList cellToProc = method.decompose(mesh_, cellWeights);

    Info<< "\nFinished decomposition into "
        << method.nDomains() << " domains in "
        << decompositionTime.elapsedCpuTime() << " s" << nl << nl;

    decompositionInformation info
    (
        cellCells,
        cellToProc,
        method.nDomains()
    );

    if (writeCellDist)
    {
        // Write decomposition for visualization
        // - write as VTU to avoid any impact
        writeVTK("cellDist", cellToProc);

// Less likely that this is actually required, but may be useful...
//
//        // Write decomposition as labelList for use with 'manual'
//        // decomposition method.
//        labelIOList cellDecomposition
//        (
//            IOobject
//            (
//                "cellDecomposition",
//                mesh_.facesInstance(),
//                mesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE,
//                IOobject::NO_REGISTER
//            ),
//            std::move(cellToProc)
//        );
//        cellDecomposition.write();
//
//        Info<< nl << "Wrote decomposition to "
//            << cellDecomposition.objectRelPath()
//            << " for use in manual decomposition." << endl;

        Info<< nl;
    }

    if (verbose)
    {
        info.printDetails(Info);
        Info<< nl;
    }
    info.printSummary(Info);
}


// ************************************************************************* //
