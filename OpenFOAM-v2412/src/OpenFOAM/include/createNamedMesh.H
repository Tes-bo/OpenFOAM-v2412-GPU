/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a fvMesh for a specified region, or the defaultRegion

Required Classes
    - Foam::fvMesh

Required Variables
    - args [argList]
    - runTime [Time]

Provided Variables
    - regionName [word]
    - mesh [fvMesh]

\*---------------------------------------------------------------------------*/

// "getRegionOption.H"
Foam::word regionName
(
    args.getOrDefault<word>("region", Foam::polyMesh::defaultRegion)
);

{
    Foam::Info << "Create mesh";
    if (!Foam::polyMesh::regionName(regionName).empty())
    {
        Foam::Info << ' ' << regionName;
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;
}


Foam::fvMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);

mesh.init(true);   // Initialise all (lower levels and current)

Foam::Info << Foam::endl;


// ************************************************************************* //
