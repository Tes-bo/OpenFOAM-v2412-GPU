/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    mapFieldsPar

Group
    grpPreProcessingUtilities

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToMesh.H"
#include "processorPolyPatch.H"
#include "MapMeshes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const word& mapMethod,
    const word& AMIMapMethod,
    const word& procMapMethod,
    const bool subtract,
    const wordRes& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Consistently creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        AMIMapMethod,
        meshToMesh::procMapMethodNames_[procMapMethod]
    );

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const word& mapMethod,
    const word& AMIMapMethod,
    const word& procMapMethod,
    const bool subtract,
    const wordRes& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        AMIMapMethod,
        patchMap,
        cuttingPatches,
        meshToMesh::procMapMethodNames_[procMapMethod]
    );

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Map volume fields from one mesh to another"
    );

    argList::addArgument("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "Specify the source time"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "Specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "Specify the target region"
    );
    argList::addBoolOption
    (
        "consistent",
        "Source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "Specify the mapping method "
        "(direct|mapNearest|cellVolumeWeight|correctedCellVolumeWeight)"
    );
    argList::addOption
    (
        "patchMapMethod",
        "word",
        "Specify the patch mapping method (direct|mapNearest|faceAreaWeight)"
    );
    argList::addOption
    (
        "procMapMethod",
        "word",
        "Specify the processor distribution map method (AABB|LOD)"
    );
    argList::addBoolOption
    (
        "subtract",
        "Subtract mapped source from target"
    );
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to reconstruct (all by default)."
        " Eg, 'T' or '(p T U \"alpha.*\")'"
    );

    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Skip mapping lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 2106});

    argList args(argc, argv);
    #include "foamDlOpenLibs.H"

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    const auto casePath = args.get<fileName>(1);
    const fileName rootDirSource = casePath.path();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << ' ' << caseDirSource;
    word sourceRegion(polyMesh::defaultRegion);
    if (args.readIfPresent("sourceRegion", sourceRegion))
    {
        Info<< " (region " << sourceRegion << ')';
    }
    Info<< endl;

    Info<< "Target: " << rootDirTarget << ' ' << caseDirTarget;
    word targetRegion(polyMesh::defaultRegion);
    if (args.readIfPresent("targetRegion", targetRegion))
    {
        Info<< " (region " << targetRegion << ')';
    }
    Info<< endl;

    const bool consistent = args.found("consistent");


    word mapMethod = meshToMesh::interpolationMethodNames_
    [
        meshToMesh::interpolationMethod::imCellVolumeWeight
    ];

    if  (args.readIfPresent("mapMethod", mapMethod))
    {
        Info<< "Mapping method: " << mapMethod << endl;
    }

    word patchMapMethod;
    if (meshToMesh::interpolationMethodNames_.found(mapMethod))
    {
        // Lookup corresponding AMI method
        meshToMesh::interpolationMethod method =
            meshToMesh::interpolationMethodNames_[mapMethod];

        patchMapMethod = meshToMesh::interpolationMethodAMI(method);
    }

    word procMapMethod =
        meshToMesh::procMapMethodNames_
        [
            meshToMesh::procMapMethod::pmAABB
        ];

    if (args.readIfPresent("procMapMethod", procMapMethod))
    {
        Info<< "Processor map method: " << procMapMethod << endl;
    }


    // Optionally override
    if (args.readIfPresent("patchMapMethod", patchMapMethod))
    {
        Info<< "Patch mapping method: " << patchMapMethod << endl;
    }


    if (patchMapMethod.empty())
    {
        FatalErrorInFunction
            << "No valid patchMapMethod for method " << mapMethod
            << ". Please supply one through the 'patchMapMethod' option"
            << exit(FatalError);
    }

    const bool subtract = args.found("subtract");
    if (subtract)
    {
        Info<< "Subtracting mapped source field from target" << endl;
    }

    // Non-mandatory
    const wordRes selectedFields(args.getList<wordRe>("fields", false));

    const bool noLagrangian = args.found("no-lagrangian");

    #include "createTimes.H"

    HashTable<word> patchMap;
    wordList cuttingPatches;

    if (!consistent)
    {
        IOdictionary mapFieldsDict
        (
            IOobject
            (
                "mapFieldsDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        mapFieldsDict.readEntry("patchMap", patchMap);
        mapFieldsDict.readEntry("cuttingPatches", cuttingPatches);
    }

    #include "setTimeIndex.H"

    Info<< "\nCreate meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.globalData().nTotalCells() << tab
        << "Target mesh size: " << meshTarget.globalData().nTotalCells()
        << nl << endl;

    if (consistent)
    {
        mapConsistentMesh
        (
            meshSource,
            meshTarget,
            mapMethod,
            patchMapMethod,
            procMapMethod,
            subtract,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        mapSubMesh
        (
            meshSource,
            meshTarget,
            patchMap,
            cuttingPatches,
            mapMethod,
            patchMapMethod,
            procMapMethod,
            subtract,
            selectedFields,
            noLagrangian
        );
    }

    runTimeSource.printExecutionTime(Info);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
