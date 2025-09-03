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
    Code chunk for converting volume and dimensioned fields
    included by foamToEnsight.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Volume field data output
if (doBoundary || doInternal)
{
    Info<< "    volume field (";

    writeAllVolFields(ensCase, ensMesh, objects, nearCellValue);
    writeAllDimFields(ensCase, ensMesh, objects);

    Info<< " )" << nl;

    // PointData
    // - only construct pointMesh on request (it constructs edge addressing)
    if (doPointValues)
    {
        Info<< "    point  field (";
        writeAllPointFields(ensCase, ensMesh, objects);
        Info<< " )" << nl;
    }
}


// ************************************************************************* //
