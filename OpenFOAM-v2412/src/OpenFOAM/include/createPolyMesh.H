/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a polyMesh (defaultRegion only)

Required Variables
    - runTime [Time]

Provided Variables
    - mesh [polyMesh]

\*---------------------------------------------------------------------------*/

Foam::Info
    << "Create polyMesh for time = "
    << runTime.timeName() << Foam::nl;

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        Foam::polyMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

Foam::Info << Foam::endl;

// ************************************************************************* //
