/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    OBJ output of blockMesh topology blocks

\*---------------------------------------------------------------------------*/

{
    refPtr<polyMesh> topoMeshPtr(blocks.topology(true));
    const polyMesh& topoMesh = topoMeshPtr();

    // Write mesh as edges
    {
        OFstream os(runTime.path()/"blockTopology.obj");

        Info<< "Writing block structure in obj format: "
            << os.name().name() << endl;

        for (const point& p : topoMesh.points())
        {
            os << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }

        for (const edge& e : topoMesh.edges())
        {
            os << "l " << e.start() + 1 << ' ' << e.end() + 1 << nl;
        }
    }

    // Write centres of blocks
    {
        OFstream os(runTime.path()/"blockCentres.obj");

        Info<< "Writing block centres in obj format: "
            << os.name().name() << endl;

        for (const point& p : topoMesh.cellCentres())
        {
            os << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
    }
}


// ************************************************************************* //
