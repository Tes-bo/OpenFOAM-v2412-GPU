/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for post-processing conversion of cloud(s) to VTK PolyData
    (.vtp extension).

\*---------------------------------------------------------------------------*/

if (doLagrangian)
{
    const fileName cloudPrefix = (regionDir/cloud::prefix);

    wordList cloudNames = ListOps::create<word>
    (
        readDir
        (
            mesh.time().path()/mesh.time().timeName()/cloudPrefix,
            fileName::DIRECTORY
        ),
        nameOp<fileName>()
    );

    Pstream::combineReduce(cloudNames, ListOps::uniqueEqOp<word>());
    Foam::sort(cloudNames);  // Consistent order

    for (const word& cloudName : cloudNames)
    {
        IOobjectList cloudObjs(mesh, runTime.timeName(), cloudPrefix/cloudName);

        bool isCloud = false;
        if (cloudObjs.erase("coordinates"))
        {
            isCloud = true;
        }
        if (cloudObjs.erase("positions"))
        {
            isCloud = true;
        }

        if (!returnReduceOr(isCloud))
        {
            continue;
        }

        // Limited to basic IOField types
        cloudObjs.filterClasses
        (
            stringListOps::foundOp<word>(fieldTypes::basic)
        );

        // Are there cloud fields (globally)?
        if (returnReduceAnd(cloudObjs.empty()))
        {
            continue;
        }

        vtk::lagrangianWriter writer
        (
            meshProxy.baseMesh(),
            cloudName,
            writeOpts,
            // Output name for the cloud
            (
                outputDir/regionDir/cloud::prefix
              / cloudName/cloudName + timeDesc
            ),
            UPstream::parRun()
        );

        Info<< "    Lagrangian: "
            << args.relativePath(writer.output()) << nl;

        writer.writeTimeValue(mesh.time().value());
        writer.writeGeometry();

        // Begin CellData/PointData
        writer.beginParcelData();

        writer.writeFields<label>(cloudObjs);
        writer.writeFields<scalar>(cloudObjs);
        writer.writeFields<vector>(cloudObjs);
        writer.writeFields<sphericalTensor>(cloudObjs);
        writer.writeFields<symmTensor>(cloudObjs);
        writer.writeFields<tensor>(cloudObjs);

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
