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
    Code chunk for converting volume and dimensioned fields
    included by foamToVTK.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{
    using reportFields = foamToVtkReportFields;

    const label nVolFields =
    (
        (doInternal || doBoundary)
      ? objects.count(stringListOps::foundOp<word>(fieldTypes::volume))
      : 0
    );

    const label nDimFields =
    (
        (doInternal || doBoundary)
      ? objects.count(stringListOps::foundOp<word>(fieldTypes::internal))
      : 0
    );

    label nPointFields =
    (
        doPointValues
      ? objects.count(stringListOps::foundOp<word>(fieldTypes::point))
      : 0
    );


    if (doInternal || doBoundary)
    {
        reportFields::volume(Info, objects);
    }
    if (doInternal)
    {
        reportFields::internal(Info, objects);
    }


    // Setup for the vtm writer.
    // For legacy format, the information added is simply ignored.

    fileName vtmOutputBase
    (
        outputDir/regionDir/vtkName + timeDesc
    );

    // Combined internal + boundary in a vtm file
    vtk::vtmWriter vtmWriter;

    // Collect individual boundaries into a vtm file
    vtk::vtmWriter vtmBoundaries;

    // Setup the internal writer
    autoPtr<vtk::internalWriter> internalWriter;

    // Interpolator for volume and dimensioned fields
    autoPtr<volPointInterpolation> pInterp;

    if (doInternal)
    {
        if (vtuMeshCells.empty())
        {
            // Use the appropriate mesh (baseMesh or subMesh)
            vtuMeshCells.reset(meshProxy.mesh());

            if (doPointValues && vtuMeshCells.manifold())
            {
                doPointValues = false;
                nPointFields = 0;
                Warning
                    << nl
                    << "Manifold cells detected (eg, AMI) - disabling PointData"
                    << nl
                    << endl;
            }
        }

        internalWriter = autoPtr<vtk::internalWriter>::New
        (
            meshProxy.mesh(),
            vtuMeshCells,
            writeOpts,
            // The output base name for internal
            (
                writeOpts.legacy()
              ? vtmOutputBase
              : (vtmOutputBase / "internal")
            ),
            UPstream::parRun()
        );

        // No sub-block for internal
        vtmWriter.append_vtu
        (
            "internal",
            vtmOutputBase.name()/"internal"
        );

        Info<< "    Internal  : "
            << args.relativePath(internalWriter->output()) << nl;

        internalWriter->writeTimeValue(mesh.time().value());
        internalWriter->writeGeometry();

        if (doPointValues)
        {
            pInterp.reset(new volPointInterpolation(mesh));
        }
    }


    // Setup the patch writers

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    PtrList<vtk::patchWriter> patchWriters;
    PtrList<PrimitivePatchInterpolation<primitivePatch>> patchInterps;

    labelList patchIds;
    if (doBoundary)
    {
        patchIds = getSelectedPatches(patches, patchSelector);
    }

    if (oneBoundary && patchIds.size())
    {
        auto writer = autoPtr<vtk::patchWriter>::New
        (
            meshProxy.mesh(),
            patchIds,
            writeOpts,
            nearCellValue,
            // Output one patch: "boundary"
            (
                writeOpts.legacy()
              ?
                (
                    outputDir/regionDir/"boundary"
                  / (meshProxy.useSubMesh() ? meshProxy.name() : "boundary")
                  + timeDesc
                )
              : (vtmOutputBase / "boundary")
            ),
            UPstream::parRun()
        );

        // No sub-block for one-patch
        vtmWriter.append_vtp
        (
            "boundary",
            vtmOutputBase.name()/"boundary"
        );

        Info<< "    Boundaries: "
            << args.relativePath(writer->output()) << nl;

        writer->writeTimeValue(timeValue);
        writer->writeGeometry();

        // Transfer writer to list for later use
        patchWriters.resize(1);
        patchWriters.set(0, writer);

        // Avoid patchInterpolation for each sub-patch
        patchInterps.resize(1); // == nullptr
    }
    else if (patchIds.size())
    {
        patchWriters.resize(patchIds.size());
        if (doPointValues)
        {
            patchInterps.resize(patchIds.size());
        }

        label nPatchWriters = 0;
        label nPatchInterps = 0;

        for (const label patchId : patchIds)
        {
            const polyPatch& pp = patches[patchId];

            auto writer = autoPtr<vtk::patchWriter>::New
            (
                meshProxy.mesh(),
                labelList(one{}, pp.index()),
                writeOpts,
                nearCellValue,
                // Output patch: "boundary"/name
                (
                    writeOpts.legacy()
                  ?
                    (
                        outputDir/regionDir/pp.name()
                      / (meshProxy.useSubMesh() ? meshProxy.name() : pp.name())
                      + timeDesc
                    )
                  : (vtmOutputBase / "boundary" / pp.name())
                ),
                UPstream::parRun()
            );

            if (!nPatchWriters)
            {
                vtmWriter.beginBlock("boundary");
                vtmBoundaries.beginBlock("boundary");
            }

            vtmWriter.append_vtp
            (
                pp.name(),
                vtmOutputBase.name()/"boundary"/pp.name()
            );

            vtmBoundaries.append_vtp
            (
                pp.name(),
                "boundary"/pp.name()
            );

            Info<< "    Boundary  : "
                << args.relativePath(writer->output()) << nl;

            writer->writeTimeValue(timeValue);
            writer->writeGeometry();

            // Transfer writer to list for later use
            patchWriters.set(nPatchWriters++, writer);

            if (patchInterps.size())
            {
                patchInterps.set
                (
                    nPatchInterps++,
                    new PrimitivePatchInterpolation<primitivePatch>(pp)
                );
            }
        }

        if (nPatchWriters)
        {
            vtmWriter.endBlock("boundary");
        }

        patchWriters.resize(nPatchWriters);
        patchInterps.resize(nPatchInterps);
    }

    // With point-interpolation, cache fields to avoid multiple re-reading
    std::unique_ptr<objectRegistry> cacheFieldsPtr;

    if (doPointValues && (nVolFields || nDimFields))
    {
        cacheFieldsPtr.reset
        (
            new objectRegistry
            (
                IOobject
                (
                    "foamToVTK::volume",
                    runTime.timeName(),
                    runTime,
                    IOobject::NO_REGISTER
                )
            )
        );
    }


    // CellData
    {
        // Begin CellData
        if (internalWriter)
        {
            // Optionally with (cellID, procID) fields
            internalWriter->beginCellData
            (
                (withMeshIds ? 1 + (internalWriter->parallel() ? 1 : 0) : 0)
              + nVolFields + nDimFields
            );

            if (withMeshIds)
            {
                internalWriter->writeCellIDs();
                internalWriter->writeProcIDs(); // parallel only
            }
        }

        if (nVolFields || withMeshIds)
        {
            for (vtk::patchWriter& writer : patchWriters)
            {
                // Optionally with (patchID, procID) fields
                writer.beginCellData
                (
                    (withMeshIds ? 2 : 0)
                  + nVolFields
                );

                if (withMeshIds)
                {
                    writer.writePatchIDs();
                    writer.writeProcIDs();
                }
            }
        }

        writeAllVolFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true, // syncPar
            cacheFieldsPtr.get()
        );

        writeAllDimFields
        (
            internalWriter,
            meshProxy,
            objects,
            true, // syncPar
            cacheFieldsPtr.get()
        );

        // End CellData is implicit
    }


    // PointData
    // - only construct pointMesh on request since it constructs
    //   edge addressing
    if (doPointValues)
    {
        // Begin PointData
        if (internalWriter)
        {
            internalWriter->beginPointData
            (
                nVolFields + nDimFields + nPointFields
              + (withPointIds ? 1 : 0)
            );

            if (withPointIds)
            {
                internalWriter->writePointIDs();
            }
        }

        forAll(patchWriters, writeri)
        {
            const label nPatchFields =
            (
                (patchInterps.test(writeri) ? nVolFields : 0)
              + nPointFields
            );

            if (nPatchFields)
            {
                patchWriters[writeri].beginPointData(nPatchFields);
            }
        }

        writeAllVolFields
        (
            internalWriter, pInterp,
            patchWriters,   patchInterps,
            meshProxy,
            objects,
            true, // syncPar
            cacheFieldsPtr.get()
        );

        writeAllDimFields
        (
            internalWriter, pInterp,
            meshProxy,
            objects,
            true, // syncPar
            cacheFieldsPtr.get()
        );

        writeAllPointFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true  // syncPar
        );

        // End PointData is implicit
    }


    // Finish writers
    if (internalWriter)
    {
        internalWriter->close();
    }

    for (vtk::patchWriter& writer : patchWriters)
    {
        writer.close();
    }

    pInterp.clear();
    patchWriters.clear();
    patchInterps.clear();


    // Collective output

    if (UPstream::master())
    {
        // Naming for vtm, file series etc.
        fileName outputName(vtmOutputBase);

        if (writeOpts.legacy())
        {
            if (doInternal)
            {
                // Add to file-series and emit as JSON

                outputName.ext(vtk::legacy::fileExtension);

                fileName seriesName(vtk::seriesWriter::base(outputName));

                vtk::seriesWriter& series = vtkSeries(seriesName);

                // First time?
                // Load from file, verify against filesystem,
                // prune time >= currentTime
                if (series.empty())
                {
                    series.load(seriesName, true, timeValue);
                }

                series.append(timeValue, timeDesc);
                series.write(seriesName);
            }
        }
        else
        {
            if (vtmWriter.size())
            {
                // Emit ".vtm"

                outputName.ext(vtmWriter.ext());

                fileName seriesName(vtk::seriesWriter::base(outputName));

                vtmWriter.setTime(timeValue);
                vtmWriter.write(outputName);

                // Add to file-series and emit as JSON

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

                // Add to multi-region vtm
                vtmMultiRegion.add(regionName, regionDir, vtmWriter);
            }

            if (vtmBoundaries.size())
            {
                // Emit "boundary.vtm" with collection of boundaries

                // Naming for vtm
                fileName outputName(vtmOutputBase / "boundary");
                outputName.ext(vtmBoundaries.ext());

                vtmBoundaries.setTime(timeValue);
                vtmBoundaries.write(outputName);
            }
        }
    }
}


// ************************************************************************* //
