/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a dynamicFvMesh for a specified region, or the defaultRegion

Required Classes
    - Foam::dynamicFvMesh

Required Variables
    - args [argList]
    - runTime [Time]

Provided Variables
    - regionName [word]
    - mesh [dynamicFvMesh]
    - meshPtr [autoPtr<dynamicFvMesh>]

\*---------------------------------------------------------------------------*/

Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;

// "getRegionOption.H"
Foam::word regionName
(
    args.getOrDefault<word>("region", Foam::polyMesh::defaultRegion)
);

{
    Foam::Info << "Create dynamic mesh";
    if (!Foam::polyMesh::regionName(regionName).empty())
    {
        Foam::Info << ' ' << regionName;
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;

    meshPtr = Foam::dynamicFvMesh::New
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
}

Foam::dynamicFvMesh& mesh = meshPtr();

Foam::Info << Foam::endl;


// ************************************************************************* //
