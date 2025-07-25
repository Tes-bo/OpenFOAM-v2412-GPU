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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Search for the appropriate faMeshDefinition dictionary...

Required Classes
    - Foam::polyMesh
    - Foam::IOdictionary

Required Variables
    - regionName [word]  (the polyMesh region)
    - areaRegionName [word] (the areaMesh region)
    - args [argList]
    - runTime [Time]

Provided Variables
    - meshDefDict [IOdictionary]
    - meshDictPtr [autoPtr<IOdictionary>]

\*---------------------------------------------------------------------------*/

const word dictName("faMeshDefinition");

autoPtr<IOdictionary> meshDictPtr;

{
    fileName dictPath;
    const word& regionDir = Foam::polyMesh::regionName(regionName);
    const word& areaRegionDir = Foam::polyMesh::regionName(areaRegionName);

    if (args.readIfPresent("dict", dictPath))
    {
        // Dictionary specified on the command-line ...

        if (isDir(dictPath))
        {
            dictPath /= dictName;
        }
    }
    else if
    (
        // Dictionary under system/faMeshDefinition ?
        // (v2312 and earlier)

        areaRegionDir.empty()
     && exists
        (
            runTime.path()/runTime.caseSystem()
          / regionDir/faMesh::meshSubDir/dictName
        )
    )
    {
        // Dictionary present directly in system/ (v2312 and earlier)

        dictPath = runTime.system()/regionDir/dictName;
    }
    else
    {
        // Use system/finite-area/ directory, with region qualifications

        dictPath =
        (
            runTime.system()/regionDir
          / faMesh::prefix()/areaRegionDir/dictName
        );
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER,
        true    // is globalObject
    );

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << meshDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< "Creating faMesh from definition: "
        << meshDictIO.objectRelPath() << endl;

    meshDictPtr = autoPtr<IOdictionary>::New(meshDictIO);
}

IOdictionary& meshDefDict = *meshDictPtr;


// ************************************************************************* //
