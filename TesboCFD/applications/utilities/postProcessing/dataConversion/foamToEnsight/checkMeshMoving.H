/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Check for "points" in any of the result directories (each region)

\*---------------------------------------------------------------------------*/

// All regions moving, or no regions moving. Do not mix.
bool hasMovingMesh(false);

if (timeDirs.size() > 1)
{
    // We already loaded meshes (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    // Do all regions as moving, or all as static.

    boolList isMoving(meshes.size(), false);
    label nMoving = 0;

    Info<< "Search for moving mesh ... " << flush;
    for (const instant& inst : timeDirs)
    {
        const word& timeName = inst.name();

        forAll(meshes, regioni)
        {
            const fvMesh& mesh = meshes[regioni];

            if
            (
                !isMoving[regioni]
             && (timeName != mesh.pointsInstance())
             && IOobject
                (
                    "points",
                    timeName,
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ).typeHeaderOk<pointIOField>(true, false)
            )
            {
                isMoving[regioni] = true;
                ++nMoving;
            }
        }

        if (nMoving == isMoving.size())
        {
            break;
        }
    }

    if (nMoving)
    {
        Info<< "found " << nMoving
            << " moving regions. Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }

    hasMovingMesh = (nMoving != 0);
}

// Ensure consistency
Pstream::reduceOr(hasMovingMesh);

// ************************************************************************* //
