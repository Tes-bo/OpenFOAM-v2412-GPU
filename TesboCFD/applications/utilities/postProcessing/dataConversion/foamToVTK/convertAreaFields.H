/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for converting finite-area - included by foamToVTK.

    Typedef vtkWriterType_areaMesh is declared in writeAreaFields.H

\*---------------------------------------------------------------------------*/

//
// Finite-area mesh and fields - need not exist
//

// No subsetting!
if (doFiniteArea)
{
    using reportFields = foamToVtkReportFields;

    autoPtr<faMesh> faMeshPtr;

    const label nAreaFields =
        faObjects.count(stringListOps::foundOp<word>(fieldTypes::area));

    if (nAreaFields || withMeshIds)
    {
        faMeshPtr = faMesh::TryNew(meshProxy.baseMesh());
    }

    if (faMeshPtr && (nAreaFields || withMeshIds))
    {
        const faMesh& areaMesh = faMeshPtr();

        reportFields::area(Info, faObjects);

        const auto& pp = faMeshPtr->patch();

        vtkWriterType_areaMesh writer
        (
            pp,
            writeOpts,
            (
                outputDir/regionDir/"finite-area"
              / "finiteArea" + timeDesc
            ),
            UPstream::parRun()
        );
        Info<< "    Area      : "
            << args.relativePath(writer.output()) << nl;

        writer.beginFile(areaMesh.name());

        writer.writeTimeValue(timeValue);
        writer.writeGeometry();

        // Optionally with (cellID, patchID, faceLabels, procID) fields
        writer.beginCellData
        (
            (withMeshIds ? 3 + (writer.parallel() ? 1 : 0) : 0)
          + nAreaFields
        );

        if (withMeshIds)
        {
            const globalIndex procAddr(areaMesh.nFaces());

            // Use global indexed values for the 'cell' ids
            writer.writeCellData("cellID", identity(procAddr.range()));

            // The patch ids can also be quite useful
            const polyBoundaryMesh& pbm = areaMesh.mesh().boundaryMesh();

            labelList patchIds
            (
                pbm.patchID(areaMesh.faceLabels())
            );

            writer.writeCellData("patchID", patchIds);

            // Use proc-local data for faceLabels
            // (confusing enough already without renumbering)
            writer.writeCellData("faceLabels", areaMesh.faceLabels());

            writer.writeProcIDs();  // parallel only
        }

        writeAllAreaFields
        (
            writer,
            areaMesh,
            faObjects,
            true // syncPar
        );

        fileName outputName(writer.output());

        writer.close();

        if (UPstream::master())
        {
            // Add to file-series and emit as JSON

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }
    }
}


// ************************************************************************* //
