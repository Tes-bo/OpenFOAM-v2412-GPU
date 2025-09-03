/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create single or multiple regions based on the wordList 'regionNames'

Required Variables
    - runTime [Time]
    - regionNames [wordList]

Provided Variables
    - meshes [PtrList<fvMesh>]

See Also
    addAllRegionOptions.H
    getAllRegionOptions.H

\*---------------------------------------------------------------------------*/

Foam::PtrList<Foam::fvMesh> meshes(regionNames.size());

{
    forAll(regionNames, regioni)
    {
        const Foam::word& regionName = regionNames[regioni];

        Foam::Info<< "Create mesh";
        if (regionNames.size() > 1)
        {
            Foam::Info<< ' ' << regionName;
        }
        Foam::Info<< " for time = " << runTime.timeName() << Foam::nl;

        meshes.set
        (
            regioni,
            new Foam::fvMesh
            (
                Foam::IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                ),
                false  // Do not initialise
            )
        );
    }

    for (auto& mesh : meshes)
    {
        mesh.init(true);  // Initialise all (lower levels and current)
    }

    Foam::Info<< Foam::endl;
}


// ************************************************************************* //
