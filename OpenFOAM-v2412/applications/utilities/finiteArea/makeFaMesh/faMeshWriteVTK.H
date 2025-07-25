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
    VTK output of faMesh with some geometric or debug fields

\*---------------------------------------------------------------------------*/

{
    // finiteArea
    vtk::uindirectPatchWriter writer
    (
        aMesh.patch(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.time().globalPath() / "finiteArea"
        )
    );

    writer.writeGeometry();

    // CellData
    writer.beginCellData(3);
    writer.writeProcIDs();
    writer.write("area", aMesh.S().field());
    writer.write("normal", aMesh.faceAreaNormals());

    // PointData
    writer.beginPointData(1);
    writer.write("normal", aMesh.pointAreaNormals());

    Info<< nl
        << "Wrote faMesh in vtk format: " << writer.output().name() << nl;
}


// ************************************************************************* //
