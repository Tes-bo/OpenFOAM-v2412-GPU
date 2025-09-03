/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022,2024 OpenCFD Ltd.
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
    decomposePar

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTIONS]

    Options:
      - \par -allRegions
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -case \<dir\>
        Specify case directory to use (instead of the cwd).

      - \par -cellDist
        Write the cell distribution as a labelList, for use with 'manual'
        decomposition method and as a VTK or volScalarField for visualization.

      - \par -constant
        Include the 'constant/' dir in the times list.

      - \par -copyUniform
        Copy any \a uniform directories too.

      - \par -copyZero
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -debug-switch \<name=val\>
        Specify the value of a registered debug switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -decomposeParDict \<file\>
        Use specified file for decomposePar dictionary.

      - \par -dry-run
        Test without writing the decomposition. Changes -cellDist to
        only write VTK output.

      - \par -fields
        Use existing geometry decomposition and convert fields only.

      - \par fileHandler \<handler\>
        Override the file handler type.

      - \par -force
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

      - \par -ifRequired
        Only decompose the geometry if the number of domains has changed from a
        previous decomposition. No \a processor subdirectories will be removed
        unless the \a -force option is also specified. This option can be used
        to avoid redundant geometry decomposition (eg, in scripts), but should
        be used with caution when the underlying (serial) geometry or the
        decomposition method etc. have been changed between decompositions.

      - \par -info-switch \<name=val\>
        Specify the value of a registered info switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -latestTime
        Select the latest time.

      - \par -lib \<name\>
        Additional library or library list to load (can be used multiple times).

      - \par -noFunctionObjects
        Do not execute function objects.

      - \par -noSets
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -noZero
        Exclude the \a 0 dir from the times list.

      - \par -opt-switch \<name=val\>
        Specify the value of a registered optimisation switch (int/bool).
        Default is 1 if the value is omitted. (Can be used multiple times)

      - \par -region \<regionName\>
        Decompose named region. Does not check for existence of processor*.

      - \par -time \<ranges\>
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases. Eg, ':10,20 40:70 1000:', 'none', etc.

      - \par -verbose
        Additional verbosity.

      - \par -doc
        Display documentation in browser.

      - \par -doc-source
        Display source code in browser.

      - \par -help
        Display short help and exit.

      - \par -help-man
        Display full help (manpage format) and exit.

      - \par -help-notes
        Display help notes (description) and exit.

      - \par -help-full
        Display full help and exit.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "OSspecific.H"
#include "IOobjectList.H"

#include "decompositionModel.H"
#include "domainDecomposition.H"
#include "domainDecompositionDryRun.H"

#include "regionProperties.H"

#include "fieldsDistributor.H"

#include "fvFieldDecomposer.H"
#include "pointFields.H"
#include "pointFieldDecomposer.H"

#include "lagrangianFieldDecomposer.H"

#include "emptyFaPatch.H"
#include "faFieldDecomposer.H"
#include "faMeshDecomposition.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read proc addressing at specific instance.
// Uses polyMesh/fvMesh meshSubDir by default
autoPtr<labelIOList> procAddressing
(
    const objectRegistry& procRegistry,
    const word& name,
    const word& instance,
    const word& local = polyMesh::meshSubDir
)
{
    return autoPtr<labelIOList>::New
    (
        IOobject
        (
            name,
            instance,
            local,
            procRegistry,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
}


// Read proc addressing at specific instance.
// Uses the finiteArea meshSubDir
autoPtr<labelIOList> faProcAddressing
(
    const objectRegistry& procRegistry,
    const word& name,
    const word& instance,
    const word& local = faMesh::meshSubDir
)
{
    return procAddressing(procRegistry, name, instance, local);
}


// Return cached or read proc addressing from facesInstance
const labelIOList& procAddressing
(
    const PtrList<fvMesh>& procMeshList,
    const label proci,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const auto& procMesh = procMeshList[proci];

    return procAddressingList.try_emplace
    (
        proci,
        IOobject
        (
            name,
            procMesh.facesInstance(),
            polyMesh::meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
}


void decomposeUniform
(
    const bool copyUniform,
    const domainDecomposition& mesh,
    const Time& processorDb,
    const word& regionDir = word::null
)
{
    const Time& runTime = mesh.time();

    // Any uniform data to copy/link?
    const fileName uniformDir(regionDir/"uniform");

    if (fileHandler().isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.timePath()/uniformDir
            << endl;

        // Bit of trickery to synthesise the correct directory base,
        // e.g. processors4/0.01
        const fileName timePath = fileHandler().objectPath
        (
            IOobject
            (
                "dummy",
                runTime.timeName(),
                processorDb
            ),
            word::null
        ).path();

        // If no fields have been decomposed the destination
        // directory will not have been created so make sure.
        mkDir(timePath);

        if (copyUniform || mesh.distributed())
        {
            if (!fileHandler().exists(timePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
        }
        else
        {
            // Link with relative paths
            string parentPath = string("..")/"..";

            if (!regionDir.empty())
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());
            chDir(timePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTime.timeName()/uniformDir,
                    uniformDir
                );
            }
            chDir(currentDir);
        }
    }
}

} // End namespace Foam


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Alternative decomposePar dictionary file"
    );

    #include "addAllRegionOptions.H"

    argList::addDryRunOption
    (
        "Test without writing the decomposition. "
        "Changes -cellDist to only write VTK output."
    );
    argList::addVerboseOption();
    argList::addOption
    (
        "domains",
        "N",
        "Override numberOfSubdomains (-dry-run only)",
        true  // Advanced option
    );
    argList::addOption
    (
        "method",
        "name",
        "Override decomposition method (-dry-run only)",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress finiteArea mesh/field decomposition",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "no-lagrangian",
        "Suppress lagrangian (cloud) decomposition",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method and as a volScalarField for visualization."
    );
    argList::addBoolOption
    (
        "copyZero",
        "Copy 0/ directory to processor*/ rather than decompose the fields"
    );
    argList::addBoolOption
    (
        "copyUniform",
        "Copy any uniform/ directories too"
    );
    argList::addBoolOption
    (
        "fields",
        "Use existing geometry decomposition and convert fields only"
    );
    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields (volume, finite-area, lagrangian)"
    );

    argList::addBoolOption
    (
        "no-sets",
        "Skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addOptionCompat("no-sets", {"noSets", 2106});

    argList::addBoolOption
    (
        "force",
        "Remove existing processor*/ subdirs before decomposing the geometry"
    );
    argList::addBoolOption
    (
        "ifRequired",
        "Only decompose geometry if the number of domains has changed"
    );

    // Allow explicit -constant, have zero from time range
    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    #include "setRootCase.H"

    const bool writeCellDist    = args.found("cellDist");

    // Most of these are ignored for dry-run (not triggered anywhere)
    const bool copyZero         = args.found("copyZero");
    const bool copyUniform      = args.found("copyUniform");
    const bool decomposeSets    = !args.found("no-sets");

    const bool decomposeIfRequired = args.found("ifRequired");

    const bool doDecompFields = !args.found("no-fields");
    const bool doFiniteArea = !args.found("no-finite-area");
    const bool doLagrangian = !args.found("no-lagrangian");

    bool decomposeFieldsOnly = args.found("fields");
    bool forceOverwrite      = args.found("force");


    // Set time from database
    #include "createTime.H"

    // Allow override of time (unless dry-run)
    instantList times;
    if (args.dryRun())
    {
        Info<< "\ndry-run: ignoring -copy*, -fields, -force, time selection"
            << nl;
    }
    else
    {
        if (decomposeFieldsOnly && !doDecompFields)
        {
            FatalErrorIn(args.executable())
                << "Options -fields and -no-fields are mutually exclusive"
                << " ... giving up" << nl
                << exit(FatalError);
        }

        if (!doDecompFields)
        {
            Info<< "Skip decompose of all fields" << nl;
        }
        if (!doFiniteArea)
        {
            Info<< "Skip decompose of finiteArea mesh/fields" << nl;
        }
        if (!doLagrangian)
        {
            Info<< "Skip decompose of lagrangian positions/fields" << nl;
        }

        times = timeSelector::selectIfPresent(runTime, args);
    }


    // Allow override of decomposeParDict location
    fileName decompDictFile;
    if
    (
        args.readIfPresent("decomposeParDict", decompDictFile)
     && !decompDictFile.empty() && !decompDictFile.isAbsolute()
    )
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    // Get region names
    #include "getAllRegionOptions.H"

    const bool optRegions =
        (regionNames.size() != 1 || regionNames[0] != polyMesh::defaultRegion);

    if (regionNames.size() == 1 && regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Using region: " << regionNames[0] << nl << endl;
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = polyMesh::regionName(regionName);

        if (args.dryRun())
        {
            Info<< "dry-run: decomposing mesh " << regionName << nl << nl
                << "Create mesh..." << flush;

            domainDecompositionDryRun decompTest
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                decompDictFile,
                args.getOrDefault<label>("domains", 0),
                args.getOrDefault<word>("method", word::null)
            );

            decompTest.execute(writeCellDist, args.verbose());
            continue;
        }

        Info<< "\n\nDecomposing mesh";
        if (!regionDir.empty())
        {
            Info<< ' ' << regionName;
        }
        Info<< nl << endl;

        // Determine the existing processor count directly
        const label nProcsOld =
            fileHandler().nProcs(runTime.path(), regionDir);

        // Get requested numberOfSubdomains directly from the dictionary.
        // Note: have no mesh yet so cannot use decompositionModel::New
        const label nDomains = decompositionMethod::nDomains
        (
            IOdictionary
            (
                IOobject::selectIO
                (
                    IOobject
                    (
                        decompositionModel::canonicalName,
                        runTime.time().system(),
                        regionDir,  // region (if non-default)
                        runTime,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    decompDictFile
                )
            )
        );

        // Give file handler a chance to determine the output directory
        const_cast<fileOperation&>(fileHandler()).nProcs(nDomains);

        if (decomposeFieldsOnly)
        {
            // Sanity check on previously decomposed case
            if (nProcsOld != nDomains)
            {
                FatalErrorIn(args.executable())
                    << "Specified -fields, but the case was decomposed with "
                    << nProcsOld << " domains"
                    << nl
                    << "instead of " << nDomains
                    << " domains as specified in decomposeParDict" << nl
                    << exit(FatalError);
            }
        }
        else if (nProcsOld)
        {
            bool procDirsProblem = true;

            if (decomposeIfRequired && nProcsOld == nDomains)
            {
                // We can reuse the decomposition
                decomposeFieldsOnly = true;
                procDirsProblem = false;
                forceOverwrite = false;

                Info<< "Using existing processor directories" << nl;
            }

            if (optRegions)
            {
                procDirsProblem = false;
                forceOverwrite = false;
            }

            if (forceOverwrite)
            {
                Info<< "Removing " << nProcsOld
                    << " existing processor directories" << endl;

                // Remove existing processors directory
                fileNameList dirs
                (
                    fileHandler().readDir
                    (
                        runTime.path(),
                        fileName::Type::DIRECTORY
                    )
                );
                forAllReverse(dirs, diri)
                {
                    const fileName& d = dirs[diri];

                    label proci = -1;

                    if
                    (
                        d.starts_with("processor")
                     &&
                        (
                            // Collated is "processors"
                            d[9] == 's'

                            // Uncollated has integer(s) after 'processor'
                         || Foam::read(d.substr(9), proci)
                        )
                    )
                    {
                        if (fileHandler().exists(d))
                        {
                            fileHandler().rmDir(d);
                        }
                    }
                }

                procDirsProblem = false;
            }

            if (procDirsProblem)
            {
                FatalErrorIn(args.executable())
                    << "Case is already decomposed with " << nProcsOld
                    << " domains, use the -force option or manually" << nl
                    << "remove processor directories before decomposing. e.g.,"
                    << nl
                    << "    rm -rf " << runTime.path().c_str() << "/processor*"
                    << nl
                    << exit(FatalError);
            }
        }

        Info<< "Create mesh" << endl;
        domainDecomposition mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            decompDictFile
        );
        // Make sure pointMesh gets read as well
        (void)pointMesh::New(mesh, IOobject::READ_IF_PRESENT);


        // Decompose the mesh
        if (!decomposeFieldsOnly)
        {
            mesh.decomposeMesh();
            mesh.writeDecomposition(decomposeSets);

            if (writeCellDist)
            {
                const labelList& procIds = mesh.cellToProc();

                // Write decomposition for visualization
                mesh.writeVolField("cellDist");
                //TBD: mesh.writeVTK("cellDist");

                // Write decomposition as labelList for use with 'manual'
                // decomposition method.
                labelIOList cellDecomposition
                (
                    IOobject
                    (
                        "cellDecomposition",
                        mesh.facesInstance(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    procIds
                );
                cellDecomposition.write();

                Info<< nl << "Wrote decomposition to "
                    << cellDecomposition.objectRelPath()
                    << " for use in manual decomposition." << endl;
            }

            fileHandler().flush();
        }


        if (copyZero)
        {
            // Copy the 0/ directory into each of the processor directories
            // with fallback of 0.orig/ directory if necessary.

            fileName inputDir(runTime.path()/"0");

            bool canCopy = fileHandler().isDir(inputDir);
            if (!canCopy)
            {
                // Try with "0.orig" instead
                inputDir.ext("orig");
                canCopy = fileHandler().isDir(inputDir);
            }

            if (canCopy)
            {
                // Avoid copying into the same directory multiple times
                // (collated format). Don't need a hash here.
                fileName prevOutputDir;
                for (label proci = 0; proci < mesh.nProcs(); ++proci)
                {
                    Time processorDb
                    (
                        Time::controlDictName,
                        args.rootPath(),
                        args.caseName()/("processor" + Foam::name(proci)),
                        false,  // No function objects
                        false   // No extra controlDict libs
                    );
                    // processorDb.setTime(runTime);

                    // Get corresponding directory name
                    // (to handle processors/)
                    const fileName outputDir
                    (
                        fileHandler().objectPath
                        (
                            IOobject
                            (
                                word::null, // name
                                "0",        // instance (time == 0)
                                processorDb
                            ),
                            word::null
                        )
                    );

                    if (outputDir != prevOutputDir)
                    {
                        Info<< "Processor " << proci
                            << ": copying \""
                            << inputDir.name() << "/\" to "
                            << runTime.relativePath(outputDir)
                            << endl;

                        fileHandler().cp(inputDir, outputDir);
                        prevOutputDir = outputDir;
                    }
                }
            }
            else
            {
                Info<< "No 0/ or 0.orig/ directory to copy" << nl;
            }
        }
        else
        {
            // Decompose field files, lagrangian, finite-area

            // Cached processor meshes and maps. These are only preserved if
            // running with multiple times.
            PtrList<Time> processorDbList(mesh.nProcs());
            PtrList<fvMesh> procMeshList(mesh.nProcs());
            PtrList<labelIOList> faceProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> cellProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> boundaryProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> pointProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> pointBoundaryProcAddressingList(mesh.nProcs());

            PtrList<fvFieldDecomposer> fieldDecomposerList(mesh.nProcs());
            PtrList<pointFieldDecomposer> pointFieldDecomposerList
            (
                mesh.nProcs()
            );


            // Loop over all times
            forAll(times, timei)
            {
                runTime.setTime(times[timei], timei);

                Info<< "Time = " << runTime.timeName() << endl;

                // Field objects at this time
                IOobjectList objects;
                IOobjectList faObjects;

                if (doDecompFields)
                {
                    // List of volume mesh objects for this instance
                    objects = IOobjectList(mesh, runTime.timeName());

                    // List of area mesh objects (assuming single region)
                    faObjects = IOobjectList
                    (
                        mesh.time(),
                        runTime.timeName(),
                        faMesh::dbDir(mesh, word::null),
                        IOobjectOption::NO_REGISTER
                    );

                    // Ignore generated fields: (cellDist)
                    objects.remove("cellDist");
                }

                // Finite area handling
                autoPtr<faMeshDecomposition> faMeshDecompPtr;
                if (doFiniteArea)
                {
                    const word boundaryInst =
                        mesh.time().findInstance(mesh.meshDir(), "boundary");

                    IOobject io
                    (
                        "faBoundary",
                        boundaryInst,
                        faMesh::meshDir(mesh, word::null),
                        mesh.time(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    );

                    if (io.typeHeaderOk<faBoundaryMesh>(true))
                    {
                        // Always based on the volume decomposition!
                        faMeshDecompPtr.reset
                        (
                            new faMeshDecomposition
                            (
                                mesh,
                                mesh.nProcs(),
                                mesh.model()
                            )
                        );
                    }
                }


                // Volume/surface/internal fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                fvFieldDecomposer::fieldsCache volumeFieldCache;

                if (doDecompFields)
                {
                    volumeFieldCache.readAllFields(mesh, objects);
                }


                // Point fields
                // ~~~~~~~~~~~~

                // Read decomposed pointMesh
                const pointMesh& pMesh =
                    pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

                pointFieldDecomposer::fieldsCache pointFieldCache;

                if (doDecompFields)
                {
                    pointFieldCache.readAllFields(pMesh, objects);
                }


                // Lagrangian fields
                // ~~~~~~~~~~~~~~~~~

                fileNameList cloudDirs;

                if (doDecompFields && doLagrangian)
                {
                    cloudDirs = fileHandler().readDir
                    (
                        runTime.timePath()/cloud::prefix,
                        fileName::DIRECTORY
                    );
                }

                // Particles
                PtrList<Cloud<indexedParticle>> lagrangianPositions
                (
                    cloudDirs.size()
                );
                // Particles per cell
                PtrList<List<SLList<indexedParticle*>*>> cellParticles
                (
                    cloudDirs.size()
                );

                lagrangianFieldDecomposer::fieldsCache lagrangianFieldCache
                (
                    cloudDirs.size()
                );

                label cloudI = 0;

                for (const fileName& cloudDir : cloudDirs)
                {
                    IOobjectList cloudObjects
                    (
                        mesh,
                        runTime.timeName(),
                        cloud::prefix/cloudDir,
                        IOobject::NO_REGISTER
                    );

                    // Note: look up "positions" for backwards compatibility
                    if
                    (
                        cloudObjects.found("coordinates")
                     || cloudObjects.found("positions")
                    )
                    {
                        // Read lagrangian particles
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~

                        Info<< "Identified lagrangian data set: "
                            << cloudDir << endl;

                        lagrangianPositions.set
                        (
                            cloudI,
                            new Cloud<indexedParticle>
                            (
                                mesh,
                                cloudDir,
                                false
                            )
                        );


                        // Sort particles per cell
                        // ~~~~~~~~~~~~~~~~~~~~~~~

                        cellParticles.set
                        (
                            cloudI,
                            new List<SLList<indexedParticle*>*>
                            (
                                mesh.nCells(),
                                static_cast<SLList<indexedParticle*>*>(nullptr)
                            )
                        );

                        label i = 0;

                        for (indexedParticle& p : lagrangianPositions[cloudI])
                        {
                            p.index() = i++;

                            label celli = p.cell();

                            // Check
                            if (celli < 0 || celli >= mesh.nCells())
                            {
                                FatalErrorIn(args.executable())
                                    << "Illegal cell number " << celli
                                    << " for particle with index "
                                    << p.index()
                                    << " at position "
                                    << p.position() << nl
                                    << "Cell number should be between 0 and "
                                    << mesh.nCells()-1 << nl
                                    << "On this mesh the particle should"
                                    << " be in cell "
                                    << mesh.findCell(p.position())
                                    << exit(FatalError);
                            }

                            if (!cellParticles[cloudI][celli])
                            {
                                cellParticles[cloudI][celli] =
                                    new SLList<indexedParticle*>();
                            }

                            cellParticles[cloudI][celli]->append(&p);
                        }

                        // Read fields
                        // ~~~~~~~~~~~

                        IOobjectList lagrangianObjects
                        (
                            mesh,
                            runTime.timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            IOobject::NO_REGISTER
                        );

                        lagrangianFieldCache.readAllFields
                        (
                            cloudI,
                            lagrangianObjects
                        );

                        ++cloudI;
                    }
                }

                lagrangianPositions.resize(cloudI);
                cellParticles.resize(cloudI);
                lagrangianFieldCache.resize(cloudI);

                Info<< endl;

                // split the fields over processors
                for
                (
                    label proci = 0;
                    doDecompFields && proci < mesh.nProcs();
                    ++proci
                )
                {
                    Info<< "Processor " << proci << ": field transfer" << endl;

                    // open the database
                    if (!processorDbList.set(proci))
                    {
                        processorDbList.set
                        (
                            proci,
                            new Time
                            (
                                Time::controlDictName,
                                args.rootPath(),
                                args.caseName()
                              / ("processor" + Foam::name(proci)),
                                args.allowFunctionObjects(),
                                args.allowLibs()
                            )
                        );
                    }
                    Time& processorDb = processorDbList[proci];


                    processorDb.setTime(runTime);

                    // read the mesh
                    if (!procMeshList.set(proci))
                    {
                        procMeshList.set
                        (
                            proci,
                            new fvMesh
                            (
                                IOobject
                                (
                                    regionName,
                                    processorDb.timeName(),
                                    processorDb
                                )
                            )
                        );
                    }
                    const fvMesh& procMesh = procMeshList[proci];

                    const labelIOList& faceProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "faceProcAddressing",
                        faceProcAddressingList
                    );

                    const labelIOList& cellProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "cellProcAddressing",
                        cellProcAddressingList
                    );

                    const labelIOList& boundaryProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "boundaryProcAddressing",
                        boundaryProcAddressingList
                    );


                    // FV fields: volume, surface, internal
                    {
                        if (!fieldDecomposerList.set(proci))
                        {
                            fieldDecomposerList.set
                            (
                                proci,
                                new fvFieldDecomposer
                                (
                                    mesh,
                                    procMesh,
                                    faceProcAddressing,
                                    cellProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }

                        volumeFieldCache.decomposeAllFields
                        (
                            fieldDecomposerList[proci]
                        );

                        if (times.size() == 1)
                        {
                            // Clear cached decomposer
                            fieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // Point fields
                    if (!pointFieldCache.empty())
                    {
                        const labelIOList& pointProcAddressing = procAddressing
                        (
                            procMeshList,
                            proci,
                            "pointProcAddressing",
                            pointProcAddressingList
                        );

                        const pointMesh& procPMesh =
                            pointMesh::New(procMesh, IOobject::READ_IF_PRESENT);

                        if (!pointBoundaryProcAddressingList.set(proci))
                        {
                            pointBoundaryProcAddressingList.set
                            (
                                proci,
                                autoPtr<labelIOList>::New
                                (
                                    IOobject
                                    (
                                        "boundaryProcAddressing",
                                        procMesh.facesInstance(),
                                        polyMesh::meshSubDir
                                       /pointMesh::meshSubDir,
                                        procPMesh.thisDb(),
                                        IOobject::READ_IF_PRESENT,
                                        IOobject::NO_WRITE,
                                        IOobject::NO_REGISTER
                                    ),
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const auto& pointBoundaryProcAddressing =
                            pointBoundaryProcAddressingList[proci];


                        if (!pointFieldDecomposerList.set(proci))
                        {
                            pointFieldDecomposerList.set
                            (
                                proci,
                                new pointFieldDecomposer
                                (
                                    pMesh,
                                    procPMesh,
                                    pointProcAddressing,
                                    pointBoundaryProcAddressing
                                )
                            );
                        }

                        pointFieldCache.decomposeAllFields
                        (
                            pointFieldDecomposerList[proci]
                        );

                        if (times.size() == 1)
                        {
                            // Early deletion
                            pointBoundaryProcAddressingList.set
                            (
                                proci,
                                nullptr
                            );
                            pointProcAddressingList.set(proci, nullptr);
                            pointFieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // If there is lagrangian data write it out
                    forAll(lagrangianPositions, cloudi)
                    {
                        if (lagrangianPositions[cloudi].size())
                        {
                            lagrangianFieldDecomposer fieldDecomposer
                            (
                                mesh,
                                procMesh,
                                faceProcAddressing,
                                cellProcAddressing,
                                cloudDirs[cloudi],
                                lagrangianPositions[cloudi],
                                cellParticles[cloudi]
                            );

                            // Lagrangian fields
                            lagrangianFieldCache.decomposeAllFields
                            (
                                cloudi,
                                cloudDirs[cloudi],
                                fieldDecomposer
                            );
                        }
                    }

                    if (doDecompFields)
                    {
                        // Decompose "uniform" directory in the time region
                        // directory
                        decomposeUniform
                        (
                            copyUniform, mesh, processorDb, regionDir
                        );

                        // For a multi-region case, also decompose "uniform"
                        // directory in the time directory
                        if (regionNames.size() > 1 && regioni == 0)
                        {
                            decomposeUniform(copyUniform, mesh, processorDb);
                        }
                    }


                    // We have cached all the constant mesh data for the current
                    // processor. This is only important if running with
                    // multiple times, otherwise it is just extra storage.
                    if (times.size() == 1)
                    {
                        boundaryProcAddressingList.set(proci, nullptr);
                        cellProcAddressingList.set(proci, nullptr);
                        faceProcAddressingList.set(proci, nullptr);
                        procMeshList.set(proci, nullptr);
                        processorDbList.set(proci, nullptr);
                    }
                }


                // Finite area mesh and field decomposition
                if (faMeshDecompPtr)
                {
                    Info<< "\nFinite area mesh decomposition" << endl;

                    faMeshDecomposition& aMesh = faMeshDecompPtr();

                    aMesh.decomposeMesh();
                    aMesh.writeDecomposition();


                    // Area/edge fields
                    // ~~~~~~~~~~~~~~~~

                    faFieldDecomposer::fieldsCache areaFieldCache;

                    if (doDecompFields)
                    {
                        areaFieldCache.readAllFields(aMesh, faObjects);
                    }

                    const label nAreaFields = areaFieldCache.size();

                    Info<< endl;
                    Info<< "Finite area field transfer: "
                        << nAreaFields << " fields" << endl;

                    // Split the fields over processors
                    for
                    (
                        label proci = 0;
                        nAreaFields && proci < mesh.nProcs();
                        ++proci
                    )
                    {
                        Info<< "    Processor " << proci << endl;

                        // open the database
                        Time processorDb
                        (
                            Time::controlDictName,
                            args.rootPath(),
                            args.caseName()/("processor" + Foam::name(proci)),
                            false,  // No function objects
                            false   // No extra controlDict libs
                        );

                        processorDb.setTime(runTime);

                        // Read the mesh
                        fvMesh procFvMesh
                        (
                            IOobject
                            (
                                regionName,
                                processorDb.timeName(),
                                processorDb
                            )
                        );

                        faMesh procMesh(procFvMesh);

                        // // Does not work.  HJ, 15/Aug/2017
                        // const labelIOList& faceProcAddressing =
                        //     procAddressing
                        //     (
                        //         procMeshList,
                        //         proci,
                        //         "faceProcAddressing",
                        //         faceProcAddressingList
                        //     );

                        // const labelIOList& boundaryProcAddressing =
                        //     procAddressing
                        //     (
                        //         procMeshList,
                        //         proci,
                        //         "boundaryProcAddressing",
                        //         boundaryProcAddressingList
                        //     );

                        // Addressing from faMesh (not polyMesh) meshSubDir

                        autoPtr<labelIOList> tfaceProcAddr =
                            faProcAddressing
                            (
                                procMesh,
                                "faceProcAddressing",
                                runTime.constant()
                            );
                        auto& faceProcAddressing = *tfaceProcAddr;

                        autoPtr<labelIOList> tboundaryProcAddr =
                            faProcAddressing
                            (
                                procMesh,
                                "boundaryProcAddressing",
                                runTime.constant()
                            );
                        auto& boundaryProcAddressing = *tboundaryProcAddr;

                        autoPtr<labelIOList> tedgeProcAddr =
                            faProcAddressing
                            (
                                procMesh,
                                "edgeProcAddressing",
                                runTime.constant()
                            );
                        const auto& edgeProcAddressing = *tedgeProcAddr;

                        faFieldDecomposer fieldDecomposer
                        (
                            aMesh,
                            procMesh,
                            edgeProcAddressing,
                            faceProcAddressing,
                            boundaryProcAddressing
                        );

                        areaFieldCache.decomposeAllFields(fieldDecomposer);
                    }
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
