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
    Check field availability for last time.
    Done to avoid mapping 'undefined' when a field only exists as time 0.

Requires
    readFields.H (for the checkData function)

\*---------------------------------------------------------------------------*/

// Initially all possible objects that are available at the final time
List<wordHashSet> availableRegionObjectNames(meshes.size());
List<wordHashSet> availableFaRegionObjectNames(meshes.size());

forAll(meshes, regioni)
{
    const auto& mesh = meshes[regioni];

    IOobjectList objects;
    IOobjectList faObjects;

    if (doConvertFields && !timeDirs.empty())
    {
        // List of volume mesh objects for this instance
        objects = IOobjectList(mesh, timeDirs.back().name());

        // List of area mesh objects (assuming single region)
        faObjects = IOobjectList
        (
            mesh.time(),
            timeDirs.back().name(),
            faMesh::dbDir(mesh, word::null),
            IOobjectOption::NO_REGISTER
        );

        if (fieldSelector && !fieldSelector().empty())
        {
            objects.filterObjects(fieldSelector());
            faObjects.filterObjects(fieldSelector());
        }

        // Remove "*_0" restart fields
        objects.prune_0();
        faObjects.prune_0();

        if (!doPointValues)
        {
            // Prune point fields if disabled
            objects.filterClasses
            (
                [](const word& clsName)
                {
                    return fieldTypes::point.found(clsName);
                },
                true // prune
            );
        }
    }

    // Volume fields
    if (!objects.empty())
    {
        wordList objectNames(objects.sortedNames());

        // Check availability for all times...
        checkData
        (
            mesh.thisDb(),
            timeDirs,
            objectNames
        );

        availableRegionObjectNames[regioni] = objectNames;
    }

    // Area fields
    if (!faObjects.empty())
    {
        wordList objectNames(faObjects.sortedNames());

        // Check availability for all times... (assuming single region)
        checkData
        (
            mesh.time(),
            timeDirs,
            objectNames,
            faMesh::dbDir(mesh, word::null)
        );

        availableFaRegionObjectNames[regioni] = objectNames;
    }
}


// ************************************************************************* //
