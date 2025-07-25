/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reconstructPar

Group
    grpParallelUtilities

Description
    Reconstructs fields of a case that is decomposed for parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorMeshes.H"
#include "regionProperties.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "lagrangianReconstructor.H"

#include "faCFD.H"
#include "faMesh.H"
#include "processorFaMeshes.H"
#include "faFieldReconstructor.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

#include "hexRef8Data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool haveAllTimes
(
    const wordHashSet& masterTimeDirSet,
    const instantList& timeDirs
)
{
    // Loop over all times
    for (const instant& t : timeDirs)
    {
        if (!masterTimeDirSet.found(t.name()))
        {
            return false;
        }
    }
    return true;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reconstruct fields of a parallel case"
    );

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);  // constant(true), zero(true)
    argList::noParallel();

    #include "addAllRegionOptions.H"

    argList::addVerboseOption();
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to reconstruct (all by default)."
        " Eg, 'T' or '(p T U \"alpha.*\")'"
    );
    argList::addBoolOption
    (
        "no-fields",  // noFields
        "Skip reconstructing fields"
    );
    argList::addOptionCompat("no-fields", {"noFields", 2106});
    argList::addOption
    (
        "lagrangianFields",
        "wordRes",
        "Specify single or multiple lagrangian fields to reconstruct"
        " (all by default)."
        " Eg, '(U d)'"
        " - Positions are always included."
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Skip reconstructing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 2106});

    argList::addBoolOption
    (
        "no-sets",
        "Skip reconstructing cellSets, faceSets, pointSets"
    );
    argList::addOptionCompat("no-sets", {"noSets", 2106});

    argList::addBoolOption
    (
        "newTimes",
        "Only reconstruct new times (i.e. that do not exist already)"
    );

    #include "setRootCase.H"
    #include "createTime.H"


    const bool doFields = !args.found("no-fields");
    wordRes selectedFields;

    if (doFields)
    {
        args.readListIfPresent<wordRe>("fields", selectedFields);
    }
    else
    {
        Info<< "Skipping reconstructing fields";
        if (args.found("fields"))
        {
            Info<< ". Ignore -fields option";
        }
        Info<< nl << endl;
    }


    const bool doFiniteArea = !args.found("no-finite-area");
    if (!doFiniteArea)
    {
        Info<< "Skipping reconstructing finiteArea mesh/fields"
            << nl << endl;
    }


    const bool doLagrangian = !args.found("no-lagrangian");
    wordRes selectedLagrangianFields;

    if (doLagrangian)
    {
        args.readListIfPresent<wordRe>
        (
            "lagrangianFields", selectedLagrangianFields
        );
    }
    else
    {
        Info<< "Skipping reconstructing lagrangian positions/fields";
        if (args.found("lagrangianFields"))
        {
            Info<< ". Ignore -lagrangianFields option";
        }
        Info<< nl << endl;
    }


    const bool doReconstructSets = !args.found("no-sets");

    if (!doReconstructSets)
    {
        Info<< "Skipping reconstructing cellSets, faceSets and pointSets"
            << nl << endl;
    }

    const bool newTimes = args.found("newTimes");

    // Get region names
    #include "getAllRegionOptions.H"

    // Determine the processor count
    label nProcs{0};

    if (regionNames.empty())
    {
        FatalErrorInFunction
            << "No regions specified or detected."
            << exit(FatalError);
    }
    else if (regionNames[0] == polyMesh::defaultRegion)
    {
        nProcs = fileHandler().nProcs(args.path());
    }
    else
    {
        nProcs = fileHandler().nProcs(args.path(), regionNames[0]);

        if (regionNames.size() == 1)
        {
            Info<< "Using region: " << regionNames[0] << nl << endl;
        }
    }

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).nProcs(nProcs);

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci)),
                args.allowFunctionObjects(),
                args.allowLibs()
            )
        );
    }

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    // Note that we do not set the runTime time so it is still the
    // one set through the controlDict. The -time option
    // only affects the selected set of times from processor0.
    // - can be illogical
    // + any point motion handled through mesh.readUpdate


    if (timeDirs.empty())
    {
        WarningInFunction << "No times selected";
        exit(1);
    }


    // Get current times if -newTimes
    instantList masterTimeDirs;
    if (newTimes)
    {
        masterTimeDirs = runTime.times();
    }
    wordHashSet masterTimeDirSet(2*masterTimeDirs.size());
    for (const instant& t : masterTimeDirs)
    {
        masterTimeDirSet.insert(t.name());
    }


    // Set all times on processor meshes equal to reconstructed mesh
    forAll(databases, proci)
    {
        databases[proci].setTime(runTime);
    }


    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = polyMesh::regionName(regionName);

        Info<< "\n\nReconstructing fields" << nl
            << "region=" << regionName << nl << endl;

        if
        (
            newTimes
         && regionNames.size() == 1
         && regionDir.empty()
         && haveAllTimes(masterTimeDirSet, timeDirs)
        )
        {
            Info<< "Skipping region " << regionName
                << " since already have all times"
                << endl << endl;
            continue;
        }


        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );


        // Read all meshes and addressing to reconstructed mesh
        processorMeshes procMeshes(databases, regionName);

        // Loop over all times
        forAll(timeDirs, timei)
        {
            if (newTimes && masterTimeDirSet.found(timeDirs[timei].name()))
            {
                Info<< "Skipping time " << timeDirs[timei].name()
                    << endl << endl;
                continue;
            }


            // Set time for global database
            runTime.setTime(timeDirs[timei], timei);

            Info<< "Time = " << runTime.timeName() << endl << endl;

            // Set time for all databases
            forAll(databases, proci)
            {
                databases[proci].setTime(timeDirs[timei], timei);
            }

            // Check if any new meshes need to be read.
            polyMesh::readUpdateState meshStat = mesh.readUpdate();

            polyMesh::readUpdateState procStat = procMeshes.readUpdate();

            if (procStat == polyMesh::POINTS_MOVED)
            {
                // Reconstruct the points for moving mesh cases and write
                // them out
                procMeshes.reconstructPoints(mesh);
            }
            else if (meshStat != procStat)
            {
                WarningInFunction
                    << "readUpdate for the reconstructed mesh:"
                    << meshStat << nl
                    << "readUpdate for the processor meshes  :"
                    << procStat << nl
                    << "These should be equal or your addressing"
                    << " might be incorrect."
                    << " Please check your time directories for any "
                    << "mesh directories." << endl;
            }


            // Get list of objects from processor0 database
            IOobjectList objects
            (
                procMeshes.meshes()[0],
                databases[0].timeName(),
                IOobjectOption::NO_REGISTER
            );

            IOobjectList faObjects;

            if (doFiniteArea && doFields)
            {
                // List of area mesh objects (assuming single region)
                // - scan on processor0
                faObjects = IOobjectList
                (
                    procMeshes.meshes()[0],
                    databases[0].timeName(),
                    faMesh::dbDir(word::null),  // local relative to mesh
                    IOobjectOption::NO_REGISTER
                );
            }

            if (doFields)
            {
                // If there are any FV fields, reconstruct them
                Info<< "Reconstructing FV fields" << nl << endl;

                fvFieldReconstructor reconstructor
                (
                    mesh,
                    procMeshes.meshes(),
                    procMeshes.faceProcAddressing(),
                    procMeshes.cellProcAddressing(),
                    procMeshes.boundaryProcAddressing()
                );

                reconstructor.reconstructAllFields(objects, selectedFields);

                if (reconstructor.nReconstructed() == 0)
                {
                    Info<< "No FV fields" << nl << endl;
                }
            }

            if (doFields)
            {
                Info<< "Reconstructing point fields" << nl << endl;

                const pointMesh& pMesh = pointMesh::New
                (
                    mesh,
                    IOobject::READ_IF_PRESENT
                );

                pointFieldReconstructor reconstructor
                (
                    pMesh,
                    procMeshes.pointMeshes(),
                    procMeshes.pointProcAddressing(),
                    procMeshes.pointMeshBoundaryProcAddressing()
                );

                reconstructor.reconstructAllFields(objects, selectedFields);

                if (reconstructor.nReconstructed() == 0)
                {
                    Info<< "No point fields" << nl << endl;
                }
            }


            // If there are any clouds, reconstruct them.
            // The problem is that a cloud of size zero will not get written so
            // in pass 1 we determine the cloud names and per cloud name the
            // fields. Note that the fields are stored as IOobjectList from
            // the first processor that has them. They are in pass2 only used
            // for name and type (scalar, vector etc).

            if (doLagrangian)
            {
                HashTable<IOobjectList> allCloudObjects;

                forAll(databases, proci)
                {
                    fileName lagrangianDir
                    (
                        fileHandler().filePath
                        (
                            databases[proci].timePath()
                          / regionDir
                          / cloud::prefix
                        )
                    );

                    fileNameList cloudDirs;
                    if (!lagrangianDir.empty())
                    {
                        cloudDirs = fileHandler().readDir
                        (
                            lagrangianDir,
                            fileName::DIRECTORY
                        );
                    }

                    for (const fileName& cloudDir : cloudDirs)
                    {
                        // Check if we already have cloud objects for this
                        // cloudname
                        if (!allCloudObjects.found(cloudDir))
                        {
                            // Do local scan for valid cloud objects
                            IOobjectList localObjs
                            (
                                procMeshes.meshes()[proci],
                                databases[proci].timeName(),
                                cloud::prefix/cloudDir
                            );

                            if
                            (
                                localObjs.found("coordinates")
                             || localObjs.found("positions")
                            )
                            {
                                allCloudObjects.insert(cloudDir, localObjs);
                            }
                        }
                    }
                }


                if (allCloudObjects.size())
                {
                    lagrangianReconstructor reconstructor
                    (
                        mesh,
                        procMeshes.meshes(),
                        procMeshes.faceProcAddressing(),
                        procMeshes.cellProcAddressing()
                    );

                    // Pass2: reconstruct the cloud
                    forAllConstIters(allCloudObjects, iter)
                    {
                        const word cloudName = word::validate(iter.key());

                        // Objects (on arbitrary processor)
                        const IOobjectList& cloudObjs = iter.val();

                        Info<< "Reconstructing lagrangian fields for cloud "
                            << cloudName << nl << endl;

                        reconstructor.reconstructPositions(cloudName);

                        reconstructor.reconstructAllFields
                        (
                            cloudName,
                            cloudObjs,
                            selectedLagrangianFields
                        );
                    }
                }
                else
                {
                    Info<< "No lagrangian fields" << nl << endl;
                }
            }


            // If there are any FA fields, reconstruct them

            if (!doFiniteArea)
            {
            }
            else if
            (
                faObjects.count<areaScalarField>()
             || faObjects.count<areaVectorField>()
             || faObjects.count<areaSphericalTensorField>()
             || faObjects.count<areaSymmTensorField>()
             || faObjects.count<areaTensorField>()
             || faObjects.count<edgeScalarField>()
            )
            {
                Info << "Reconstructing FA fields" << nl << endl;

                faMesh aMesh(mesh);

                processorFaMeshes procFaMeshes(procMeshes.meshes());

                faFieldReconstructor reconstructor
                (
                    aMesh,
                    procFaMeshes.meshes(),
                    procFaMeshes.edgeProcAddressing(),
                    procFaMeshes.faceProcAddressing(),
                    procFaMeshes.boundaryProcAddressing()
                );

                reconstructor.reconstructAllFields(faObjects);
            }
            else
            {
                Info << "No FA fields" << nl << endl;
            }

            if (doReconstructSets)
            {
                // Scan to find all sets
                HashTable<label> cSetNames;
                HashTable<label> fSetNames;
                HashTable<label> pSetNames;

                forAll(procMeshes.meshes(), proci)
                {
                    const fvMesh& procMesh = procMeshes.meshes()[proci];

                    // Note: look at sets in current time only or between
                    // mesh and current time?. For now current time. This will
                    // miss out on sets in intermediate times that have not
                    // been reconstructed.
                    IOobjectList objects
                    (
                        procMesh,
                        databases[0].timeName(),    //procMesh.facesInstance()
                        polyMesh::meshSubDir/"sets"
                    );

                    for (const IOobject& io : objects.csorted<cellSet>())
                    {
                        cSetNames.insert(io.name(), cSetNames.size());
                    }

                    for (const IOobject& io : objects.csorted<faceSet>())
                    {
                        fSetNames.insert(io.name(), fSetNames.size());
                    }

                    for (const IOobject& io : objects.csorted<pointSet>())
                    {
                        pSetNames.insert(io.name(), pSetNames.size());
                    }
                }

                if (cSetNames.size() || fSetNames.size() || pSetNames.size())
                {
                    // Construct all sets
                    PtrList<cellSet> cellSets(cSetNames.size());
                    PtrList<faceSet> faceSets(fSetNames.size());
                    PtrList<pointSet> pointSets(pSetNames.size());

                    Info<< "Reconstructing sets:" << endl;
                    if (cSetNames.size())
                    {
                        Info<< "    cellSets "
                            << cSetNames.sortedToc() << endl;
                    }
                    if (fSetNames.size())
                    {
                        Info<< "    faceSets "
                            << fSetNames.sortedToc() << endl;
                    }
                    if (pSetNames.size())
                    {
                        Info<< "    pointSets "
                            << pSetNames.sortedToc() << endl;
                    }

                    // Load sets
                    forAll(procMeshes.meshes(), proci)
                    {
                        const fvMesh& procMesh = procMeshes.meshes()[proci];

                        IOobjectList objects
                        (
                            procMesh,
                            databases[0].timeName(),
                            polyMesh::meshSubDir/"sets"
                        );

                        // cellSets
                        const labelList& cellMap =
                            procMeshes.cellProcAddressing()[proci];

                        for (const IOobject& io : objects.csorted<cellSet>())
                        {
                            // Load cellSet
                            const cellSet procSet(io);
                            const label seti = cSetNames[io.name()];
                            if (!cellSets.set(seti))
                            {
                                cellSets.set
                                (
                                    seti,
                                    new cellSet
                                    (
                                        mesh,
                                        io.name(),
                                        procSet.size()
                                    )
                                );
                            }
                            cellSet& cSet = cellSets[seti];
                            cSet.instance() = runTime.timeName();

                            for (const label celli : procSet)
                            {
                                cSet.insert(cellMap[celli]);
                            }
                        }

                        // faceSets
                        const labelList& faceMap =
                            procMeshes.faceProcAddressing()[proci];

                        for (const IOobject& io : objects.csorted<faceSet>())
                        {
                            // Load faceSet
                            const faceSet procSet(io);
                            const label seti = fSetNames[io.name()];
                            if (!faceSets.set(seti))
                            {
                                faceSets.set
                                (
                                    seti,
                                    new faceSet
                                    (
                                        mesh,
                                        io.name(),
                                        procSet.size()
                                    )
                                );
                            }
                            faceSet& fSet = faceSets[seti];
                            fSet.instance() = runTime.timeName();

                            for (const label facei : procSet)
                            {
                                fSet.insert(mag(faceMap[facei])-1);
                            }
                        }
                        // pointSets
                        const labelList& pointMap =
                            procMeshes.pointProcAddressing()[proci];

                        for (const IOobject& io : objects.csorted<pointSet>())
                        {
                            // Load pointSet
                            const pointSet procSet(io);
                            const label seti = pSetNames[io.name()];
                            if (!pointSets.set(seti))
                            {
                                pointSets.set
                                (
                                    seti,
                                    new pointSet
                                    (
                                        mesh,
                                        io.name(),
                                        procSet.size()
                                    )
                                );
                            }
                            pointSet& pSet = pointSets[seti];
                            pSet.instance() = runTime.timeName();

                            for (const label pointi : procSet)
                            {
                                pSet.insert(pointMap[pointi]);
                            }
                        }
                    }

                    // Write sets

                    for (const auto& set : cellSets)
                    {
                        set.write();
                    }
                    for (const auto& set : faceSets)
                    {
                        set.write();
                    }
                    for (const auto& set : pointSets)
                    {
                        set.write();
                    }
                }


            // Reconstruct refinement data
            {
                PtrList<hexRef8Data> procData(procMeshes.meshes().size());

                forAll(procMeshes.meshes(), procI)
                {
                    const fvMesh& procMesh = procMeshes.meshes()[procI];

                    procData.set
                    (
                        procI,
                        new hexRef8Data
                        (
                            IOobject
                            (
                                "dummy",
                                procMesh.time().timeName(),
                                polyMesh::meshSubDir,
                                procMesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE,
                                IOobject::NO_REGISTER
                            )
                        )
                    );
                }

                // Combine individual parts

                const PtrList<labelIOList>& cellAddr =
                    procMeshes.cellProcAddressing();

                UPtrList<const labelList> cellMaps(cellAddr.size());
                forAll(cellAddr, i)
                {
                    cellMaps.set(i, &cellAddr[i]);
                }

                const PtrList<labelIOList>& pointAddr =
                    procMeshes.pointProcAddressing();

                UPtrList<const labelList> pointMaps(pointAddr.size());
                forAll(pointAddr, i)
                {
                    pointMaps.set(i, &pointAddr[i]);
                }

                UPtrList<const hexRef8Data> procRefs(procData.size());
                forAll(procData, i)
                {
                    procRefs.set(i, &procData[i]);
                }

                hexRef8Data
                (
                    IOobject
                    (
                        "dummy",
                        mesh.time().timeName(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    cellMaps,
                    pointMaps,
                    procRefs
                ).write();
            }
            }


            // Reconstruct refinement data
            {
                PtrList<hexRef8Data> procData(procMeshes.meshes().size());

                forAll(procMeshes.meshes(), procI)
                {
                    const fvMesh& procMesh = procMeshes.meshes()[procI];

                    procData.set
                    (
                        procI,
                        new hexRef8Data
                        (
                            IOobject
                            (
                                "dummy",
                                procMesh.time().timeName(),
                                polyMesh::meshSubDir,
                                procMesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE,
                                IOobject::NO_REGISTER
                            )
                        )
                    );
                }

                // Combine individual parts

                const PtrList<labelIOList>& cellAddr =
                    procMeshes.cellProcAddressing();

                UPtrList<const labelList> cellMaps(cellAddr.size());
                forAll(cellAddr, i)
                {
                    cellMaps.set(i, &cellAddr[i]);
                }

                const PtrList<labelIOList>& pointAddr =
                    procMeshes.pointProcAddressing();

                UPtrList<const labelList> pointMaps(pointAddr.size());
                forAll(pointAddr, i)
                {
                    pointMaps.set(i, &pointAddr[i]);
                }

                UPtrList<const hexRef8Data> procRefs(procData.size());
                forAll(procData, i)
                {
                    procRefs.set(i, &procData[i]);
                }

                hexRef8Data
                (
                    IOobject
                    (
                        "dummy",
                        mesh.time().timeName(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    cellMaps,
                    pointMaps,
                    procRefs
                ).write();
            }

            // If there is a "uniform" directory in the time region
            // directory copy from the master processor
            {
                fileName uniformDir0
                (
                    fileHandler().filePath
                    (
                        databases[0].timePath()/regionDir/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp(uniformDir0, runTime.timePath()/regionDir);
                }
            }

            // For the first region of a multi-region case additionally
            // copy the "uniform" directory in the time directory
            if (regioni == 0 && !regionDir.empty())
            {
                fileName uniformDir0
                (
                    fileHandler().filePath
                    (
                        databases[0].timePath()/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp(uniformDir0, runTime.timePath());
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
