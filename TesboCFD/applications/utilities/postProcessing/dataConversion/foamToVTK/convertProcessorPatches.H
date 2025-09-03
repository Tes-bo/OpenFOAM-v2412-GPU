/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for converting volume fields on processor boundaries,
    included by foamToVTK.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generate processor/processorN/procBoundary

{
    using reportFields = foamToVtkReportFields;

    const label nVolFields =
    (
        objects.count(stringListOps::foundOp<word>(fieldTypes::volume))
    );

    reportFields::volume(Info, objects);

    // Setup for the vtm writer.

    fileName vtmOutputBase
    (
        outputDir/regionDir/vtkName + timeDesc
    );

    // Naming
    const auto subDirNaming =
        [](const label i) -> fileName
        { return "processor" / ("processor" + Foam::name(i)); };


    // Dummy writer.
    autoPtr<vtk::internalWriter> internalWriter;

    // Setup the patch writers
    PtrList<vtk::patchWriter> patchWriters;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList patchIds =
        identity
        (
            patches.size()-patches.nNonProcessor(),
            patches.nNonProcessor()
        );

    forAll(patchIds, i)
    {
        if (!isA<processorPolyPatch>(patches[patchIds[i]]))
        {
            patchIds.resize(i);
            break;
        }
    }

    patchWriters.resize(patchIds.size());

    label nPatchWriters = 0;

    List<wordList> procPatchNames(UPstream::nProcs());
    procPatchNames[UPstream::myProcNo()].resize(patchIds.size());

    for (const label patchId : patchIds)
    {
        const polyPatch& pp = patches[patchId];

        auto writer = autoPtr<vtk::patchWriter>::New
        (
            meshProxy.mesh(),
            labelList(one{}, pp.index()),
            writeOpts,
            nearCellValue,
            (
                vtmOutputBase
              / subDirNaming(UPstream::myProcNo())
              / pp.name()
            ),
            false  // This MUST be non-parallel (serial only)
        );

        procPatchNames[UPstream::myProcNo()][nPatchWriters] = pp.name();

        writer->writeTimeValue(timeValue);
        writer->writeGeometry();

        // Transfer writer to list for later use
        patchWriters.set(nPatchWriters++, writer);
    }
    patchWriters.resize(nPatchWriters);


    Pstream::gatherList(procPatchNames);

    // CellData
    {
        for (vtk::patchWriter& writer : patchWriters)
        {
            // Optionally with patchID, procID, neighID fields
            // - use UPstream::parRun() not writer.parallel() !!
            writer.beginCellData
            (
                (withMeshIds ? 1 + (UPstream::parRun() ? 2 : 0) : 0)
              + nVolFields
            );

            if (withMeshIds)
            {
                writer.writePatchIDs();
                writer.writeProcIDs();  // parallel only
                writer.writeNeighIDs(); // parallel only
            }
        }

        writeAllVolFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true,   // syncPar
            nullptr // no field cache (object registry)
        );

        // End CellData is implicit
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

    patchWriters.clear();


    // Collective output

    const label nProcPatches = returnReduce(nPatchWriters, sumOp<label>());

    if (UPstream::master() && nProcPatches)
    {
        Info<< "Wrote " << nProcPatches << " processor boundaries from "
            << UPstream::nProcs() << " processes" << nl;


        // Collect individual boundaries into a vtm file
        vtk::vtmWriter vtmBoundaries;

        // Naming for vtm
        fileName outputName(vtmOutputBase / "processor");
        outputName.ext(vtmBoundaries.ext());

        vtmBoundaries.setTime(timeValue);

        forAll(procPatchNames, proci)
        {
            label n = 0;

            const word blockName("proc" + Foam::name(proci));
            const fileName dirName = subDirNaming(proci);

            for (const word& patchName : procPatchNames[proci])
            {
                if (!n)
                {
                    vtmBoundaries.beginBlock(blockName);
                    ++n;
                }

                vtmBoundaries.append_vtp
                (
                    patchName,
                    dirName/patchName
                );
            }

            if (n)
            {
                vtmBoundaries.endBlock();
            }
        }


        // Emit "processor.vtm" with collection of processor boundaries
        vtmBoundaries.write(outputName);
    }
}


// ************************************************************************* //
