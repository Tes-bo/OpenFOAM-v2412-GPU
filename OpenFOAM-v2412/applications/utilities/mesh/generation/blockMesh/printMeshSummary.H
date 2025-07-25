/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Summary of mesh information (eg, after blockMesh)

\*---------------------------------------------------------------------------*/

{
    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
        << "  " << "nPoints: " << mesh.nPoints() << nl
        << "  " << "nCells: " << mesh.nCells() << nl
        << "  " << "nFaces: " << mesh.nFaces() << nl
        << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

    const auto printZone =
        [](const Foam::zone& zn)
        {
            Info<< "  " << "zone " << zn.index()
                << " (size: " << zn.size()
                << ") name: " << zn.name() << nl;
        };

    if (mesh.cellZones().size())
    {
        Info<< "----------------" << nl
            << "Cell Zones" << nl
            << "----------------" << nl;

        for (const cellZone& zn : mesh.cellZones())
        {
            printZone(zn);
        }
    }
    if (mesh.faceZones().size())
    {
        Info<< "----------------" << nl
            << "Face Zones" << nl
            << "----------------" << nl;

        for (const faceZone& zn : mesh.faceZones())
        {
            printZone(zn);
        }
    }
    if (mesh.pointZones().size())
    {
        Info<< "----------------" << nl
            << "Point Zones" << nl
            << "----------------" << nl;

        for (const pointZone& zn : mesh.pointZones())
        {
            printZone(zn);
        }
    }

    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    for (const polyPatch& p : mesh.boundaryMesh())
    {
        Info<< "  " << "patch " << p.index()
            << " (start: " << p.start()
            << " size: " << p.size()
            << ") name: " << p.name()
            << nl;
    }
}


// ************************************************************************* //
