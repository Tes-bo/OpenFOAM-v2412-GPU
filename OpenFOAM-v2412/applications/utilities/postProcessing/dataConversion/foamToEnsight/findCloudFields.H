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
    Check time directories for lagrangian data.

\*---------------------------------------------------------------------------*/

// The fields for each cloud:
List<HashTable<HashTable<word>>> regionCloudFields(meshes.size());

// Identify if lagrangian data exist at any time step.
if (timeDirs.size() && doLagrangian)
{
    Info<< "Search for lagrangian ... " << flush;

    forAll(meshes, regioni)
    {
        const fvMesh& mesh = meshes[regioni];
        auto& cloudFields = regionCloudFields[regioni];

        const word& regionName = regionNames[regioni];

        const fileName cloudPrefix
        (
            polyMesh::regionName(regionName)/cloud::prefix
        );

        for (const instant& inst : timeDirs)
        {
            const word& timeName = inst.name();

            // DO NOT USE -->>  runTime.setTime(timeDirs[timeI], timeI);  <<--
            // It incurs a large overhead when done so frequently.

            fileNameList cloudDirs
            (
                readDir
                (
                    mesh.time().path()/timeName/cloudPrefix,
                    fileName::DIRECTORY
                )
            );

            for (fileName& cloudDir : cloudDirs)
            {
                const word cloudName(std::move(cloudDir));

                IOobjectList cloudObjs
                (
                    mesh,
                    timeName,
                    cloudPrefix/cloudName
                );

                // Clouds require "coordinates".
                // The "positions" are for v1706 and lower.
                // - detect and remove since these are treated specially

                bool isCloud = false;
                if (cloudObjs.erase("coordinates"))
                {
                    isCloud = true;
                }
                if (cloudObjs.erase("positions"))
                {
                    isCloud = true;
                }

                if (isCloud)
                {
                    // Save the cloud fields on a per cloud basis
                    auto& fieldsPerCloud = cloudFields(cloudName);

                    forAllConstIters(cloudObjs, fieldIter)
                    {
                        const IOobject* io = *fieldIter;

                        // Field name/type
                        fieldsPerCloud.insert
                        (
                            io->name(),
                            io->headerClassName()
                        );
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        for (auto& cloudFields : regionCloudFields)
        {
            Pstream::mapCombineReduce
            (
                cloudFields,
                HashTableOps::plusEqOp<word>()
            );
        }
    }
}

// Sorted list of cloud names
List<wordList> regionCloudNames(meshes.size());

{
    wordHashSet allRegionClouds;

    forAll(regionCloudNames, regioni)
    {
        regionCloudNames[regioni] = regionCloudFields[regioni].sortedToc();

        allRegionClouds.insert(regionCloudNames[regioni]);
    }

    const wordList allCloudNames(allRegionClouds.sortedToc());

    if (allRegionClouds.empty())
    {
        Info<< "none detected." << endl;
    }
    else
    {
        // Complete the echo information - as flatOutput
        allRegionClouds.writeList(Info) << endl;
    }
}


// ************************************************************************* //
