/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Removal of polyMesh directory

\*---------------------------------------------------------------------------*/

{
    // Shadows enclosing parameter (dictName)
    const word blockMeshDictName("blockMeshDict");

    const fileName polyMeshPath
    (
        runTime.path()/meshInstance/polyMesh::meshDir(regionName)
    );

    if (exists(polyMeshPath))
    {
        if (exists(polyMeshPath/blockMeshDictName))
        {
            Info<< "Not deleting polyMesh directory "
                << runTime.relativePath(polyMeshPath) << nl
                << "    because it contains " << blockMeshDictName << endl;
        }
        else
        {
            Info<< "Deleting polyMesh directory "
                << runTime.relativePath(polyMeshPath) << endl;
            rmDir(polyMeshPath);
        }
    }
}


// ************************************************************************* //
