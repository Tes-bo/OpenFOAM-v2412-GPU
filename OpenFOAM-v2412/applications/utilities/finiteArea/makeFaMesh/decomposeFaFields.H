/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Write proc addressing and decompose area fields (parallel only).

\*---------------------------------------------------------------------------*/

// Embed do-while to support early termination
if (doDecompose && UPstream::parRun())
do
{
    faMeshReconstructor reconstructor(aMesh, IOobjectOption::LAZY_READ);

    if (!reconstructor.good())
    {
        Info<< "Missing volume proc-addressing, "
               "cannot generate area proc-addressing." << nl
            << "Also skip decomposing area fields...."
            << endl;
        break;
    }

    reconstructor.writeMesh();          // Writes on master only
    reconstructor.writeAddressing();    // Writes per-proc

    Info<< "Wrote proc-addressing and serial mesh" << nl << endl;

    // Handle area fields
    // ------------------

    faFieldDecomposer::fieldsCache areaFieldsCache;

    const faMesh& serialMesh = reconstructor.mesh();

    if (doDecompFields)
    {
        // The serial finite-area mesh exists and is identical on all
        // processors, but its fields can only reliably be read on the
        // master (eg, running with distributed roots).
        //
        // - mark mesh fields as readable on master only (haveMeshOnProc)
        // - 'subset' entire serial mesh so that a full copy will be
        //   broadcast to other ranks (subsetterPtr)
        // - scan available IOobjects on the master only

        bitSet haveMeshOnProc;
        std::unique_ptr<faMeshSubset> subsetter;
        IOobjectList objects;

        refPtr<fileOperation> newHandler(fileOperation::NewUncollated());

        const bool oldDistributed = fileHandler().distributed();
        auto oldHandler = fileOperation::fileHandler(newHandler);
        fileHandler().distributed(true);

        if (UPstream::master())
        {
            haveMeshOnProc.set(UPstream::myProcNo());
            subsetter.reset(new faMeshSubset(serialMesh));

            const bool oldParRun = UPstream::parRun(false);

            objects = IOobjectList
            (
                serialMesh.time(),
                runTime.timeName(),
                serialMesh.dbDir(),
                IOobjectOption::NO_REGISTER
            );

            UPstream::parRun(oldParRun);
        }

        // Restore settings
        (void) fileOperation::fileHandler(oldHandler);
        fileHandler().distributed(oldDistributed);

        areaFieldsCache.readAllFields
        (
            haveMeshOnProc,
            subsetter.get(),
            serialMesh,
            objects
        );
    }

    const label nAreaFields = areaFieldsCache.size();

    if (nAreaFields)
    {
        Info<< "Decomposing " << nAreaFields << " area fields" << nl;

        faFieldDecomposer fieldDecomposer
        (
            serialMesh,
            aMesh,
            reconstructor.edgeProcAddressing(),
            reconstructor.faceProcAddressing(),
            reconstructor.boundaryProcAddressing()
        );

        // Report
        areaFieldsCache.decomposeAllFields(fieldDecomposer, true);
        Info<< endl;
    }
} while (false);


// ************************************************************************* //
