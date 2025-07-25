/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a dynamicFvMesh, with init()

Required Variables
    - runTime [Time]

Provided Variables
    - mesh [dynamicFvMesh], meshPtr

\*---------------------------------------------------------------------------*/

Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;

{
    Foam::Info << "Create dynamic mesh for time = "
        << runTime.timeName() << Foam::nl;

    meshPtr = dynamicFvMesh::New
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
}


dynamicFvMesh& mesh = meshPtr();

// Calculate initial mesh-to-mesh mapping. Note that this should be
// done under the hood, e.g. as a MeshObject
mesh.update();

Foam::Info << Foam::endl;


// ************************************************************************* //
