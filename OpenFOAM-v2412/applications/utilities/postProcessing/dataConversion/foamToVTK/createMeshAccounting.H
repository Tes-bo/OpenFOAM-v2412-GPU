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
    Additional mesh accounting (foamToVTK)

\*---------------------------------------------------------------------------*/

PtrList<fvMeshSubsetProxy> meshProxies(meshes.size());
PtrList<vtk::vtuCells> vtuMappings(meshes.size());

{
    forAll(meshes, regioni)
    {
        // Mesh subsetting, or pass through
        meshProxies.set
        (
            regioni,
            new fvMeshSubsetProxy
            (
                meshes[regioni],
                cellSubsetType,
                cellSelectionName
            )
        );

        // VTU sizing and decomposition information
        vtuMappings.set
        (
            regioni,
            new vtk::vtuCells(writeOpts, decomposePoly)
        );
    }
}


// ************************************************************************* //
