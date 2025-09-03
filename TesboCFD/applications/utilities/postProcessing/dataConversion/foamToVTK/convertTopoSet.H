/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for converting face and point sets - included by foamToVTK.

\*---------------------------------------------------------------------------*/

// flag to top-level code to signal early stop.
bool wroteTopoSet = false;

// Write faceSet (as PolyData)
if (faceSetName.size())
{
    // Load
    faceSet set(mesh, faceSetName);

    fileName outputName
    (
        outputDir/regionDir/"face-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    faceSet   : "
        << args.relativePath(outputName) << nl;

    vtk::writeFaceSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        UPstream::parRun()
    );

    wroteTopoSet = true;
}


// Write pointSet (as PolyData)
if (pointSetName.size())
{
    // Load
    pointSet set(mesh, pointSetName);

    fileName outputName
    (
        outputDir/regionDir/"point-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    pointSet  : "
        << args.relativePath(outputName) << nl;

    vtk::writePointSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        UPstream::parRun()
    );

    wroteTopoSet = true;
}


// ************************************************************************* //
