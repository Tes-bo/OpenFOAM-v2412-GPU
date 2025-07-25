/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a finite-area mesh for a specified region, or the defaultRegion

Required Classes
    - Foam::polyMesh

Required Variables
    - mesh [polyMesh]
    - runTime [Time]

Provided Variables
    - areaRegionName [word]
    - aMesh [faMesh]

\*---------------------------------------------------------------------------*/

// "getFaRegionOption.H"
Foam::word areaRegionName
(
    args.getOrDefault<word>("area-region", Foam::polyMesh::defaultRegion)
);

{
    Foam::Info << "Create finite-area mesh";
    if (!Foam::polyMesh::regionName(areaRegionName).empty())
    {
        Foam::Info << ' ' << areaRegionName;
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;
}

Foam::faMesh aMesh
(
    areaRegionName,
    mesh
);

Foam::Info << Foam::endl;


// ************************************************************************* //
