/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create single or multiple polyMesh regions based on the
    wordList 'regionNames'

Required Variables
    - runTime [Time]
    - regionNames [wordList]

Provided Variables
    - meshes [PtrList<polyMesh>]

See Also
    addAllRegionOptions.H
    getAllRegionOptions.H

\*---------------------------------------------------------------------------*/

Foam::PtrList<Foam::polyMesh> meshes(regionNames.size());

{
    forAll(regionNames, regioni)
    {
        const auto& regionName = regionNames[regioni];

        Foam::Info<< "Create polyMesh";
        if
        (
            regionNames.size() > 1
         || !Foam::polyMesh::regionName(regionName).empty()
        )
        {
            Foam::Info<< ' ' << regionName;
        }
        Foam::Info<< " for time = " << runTime.timeName() << Foam::nl;

        meshes.set
        (
            regioni,
            new Foam::polyMesh
            (
                Foam::IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    Foam::Info<< Foam::endl;
}


// ************************************************************************* //
