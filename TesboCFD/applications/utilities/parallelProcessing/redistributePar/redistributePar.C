/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    redistributePar

Group
    grpParallelUtilities

Description
    Redistributes existing decomposed mesh and fields according to the current
    settings in the decomposeParDict file.

    Must be run on maximum number of source and destination processors.
    Balances mesh and writes new mesh to new time directory.

    Can optionally run in decompose/reconstruct mode to decompose/reconstruct
    mesh and fields.

Usage
    \b redistributePar [OPTION]

    Options:
      - \par -decompose
        Remove any existing \a processor subdirectories and decomposes the
        mesh. Equivalent to running without processor subdirectories.

      - \par -reconstruct
        Reconstruct mesh and fields (like reconstructParMesh+reconstructPar).

      - \par -newTimes
        (in combination with -reconstruct) reconstruct only new times.

      - \par -dry-run
        (not in combination with -reconstruct) Test without actually
        decomposing.

      - \par -cellDist
        not in combination with -reconstruct) Write the cell distribution
        as a labelList, for use with 'manual'
        decomposition method and as a volScalarField for visualization.

      - \par -region \<regionName\>
        Distribute named region.

      - \par -allRegions
        Distribute all regions in regionProperties. Does not check for
        existence of processor*.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "sigFpe.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvMeshTools.H"
#include "fvMeshDistribute.H"
#include "fieldsDistributor.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "timeSelector.H"
#include "PstreamReduceOps.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmapDistributePolyMesh.H"
#include "IOobjectList.H"
#include "globalIndex.H"
#include "loadOrCreateMesh.H"
#include "processorFvPatchField.H"
#include "topoSet.H"
#include "regionProperties.H"

#include "parFvFieldDistributor.H"
#include "parPointFieldDistributor.H"
#include "hexRef8Data.H"
#include "meshRefinement.H"
#include "pointFields.H"

#include "faMeshSubset.H"
#include "faMeshTools.H"
#include "faMeshDistributor.H"
#include "faMeshesRegistry.H"
#include "parFaFieldDistributorCache.H"

#include "redistributeLagrangian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Use -verbose -verbose to switch on debug info. TBD.
int debug(::Foam::debug::debugSwitch("redistributePar", 0));
#define InfoOrPout (::debug ? Pout : Info.stream())


// Allocate a new file handler on valid processors only
// retaining the original IO ranks if possible
autoPtr<fileOperation>
getNewHandler(const boolUList& useProc, const bool verbose = true)
{
    autoPtr<fileOperation> handler
    (
        fileOperation::New(fileHandler(), useProc, verbose)
    );

    if (::debug && handler)
    {
        Pout<< "Allocated " << handler().info()
            << " myProcNo:" << UPstream::myProcNo(handler().comm())
            << " ptr:" << Foam::name(handler.get()) << endl;
    }

    return handler;
}


// Allocate a new file handler on valid processors only
// retaining the original IO ranks if possible
void newHandler(const boolUList& useProc, refPtr<fileOperation>& handler)
{
    if (!handler)
    {
        handler = getNewHandler(useProc);
    }
}


void createTimeDirs(const fileName& path)
{
    // Get current set of local processor's time directories. Uses
    // fileHandler
    instantList localTimeDirs(Time::findTimes(path, "constant"));

    instantList masterTimeDirs;
    if (Pstream::master())
    {
        masterTimeDirs = localTimeDirs;
    }
    Pstream::broadcast(masterTimeDirs);

    // Sync any cached times (e.g. masterUncollatedFileOperation::times_)
    // since only master would have done the findTimes
    for (const instant& t : masterTimeDirs)
    {
        if (!localTimeDirs.contains(t))
        {
            const fileName timePath(path/t.name());

            //Pout<< "Time:" << t << nl
            //    << "    raw       :" << timePath << nl
            //    << endl;
            // Bypass fileHandler
            Foam::mkDir(timePath);
        }
    }

    // Just to make sure remove all state and re-scan
    fileHandler().flush();
    (void)Time::findTimes(path, "constant");
}


void copyUniform
(
    refPtr<fileOperation>& readHandler,
    refPtr<fileOperation>& writeHandler,

    const bool reconstruct,
    const bool decompose,

    const word& readTimeName,
    const fileName& readCaseName,

    const objectRegistry& readDb,
    const objectRegistry& writeDb
)
{
    // 3 modes: reconstruct, decompose, redistribute

    // In reconstruct mode (separate reconstructed mesh):
    // - read using readDb + readHandler
    // - write using writeDb + writeHandler

    // In decompose mode (since re-using processor0 mesh):
    // - read using readDb + readCaseName + readHandler
    // - write using writeDb + writeHandler

    // In redistribute mode:
    // - read using readDb + readHandler
    // - write using writeDb + writeHandler

    fileName readPath;

    if (readHandler)
    {
        const label oldComm = UPstream::commWorld(readHandler().comm());

        Time& readTime = const_cast<Time&>(readDb.time());
        bool oldProcCase = readTime.processorCase();
        string oldCaseName;
        if (decompose)
        {
            //Pout<< "***Setting caseName to " << readCaseName
            //    << " to read undecomposed uniform" << endl;
            oldCaseName = readTime.caseName();
            readTime.caseName() = readCaseName;
            oldProcCase = readTime.processorCase(false);
        }

        // Detect uniform/ at original database + time
        readPath = readHandler().dirPath
        (
            false,          // local directory
            IOobject("uniform", readTimeName, readDb),
            false           // do not search in time
        );


        UPstream::commWorld(oldComm);

        if (decompose)
        {
            // Reset caseName on master
            //Pout<< "***Restoring caseName to " << oldCaseName << endl;
            readTime.caseName() = oldCaseName;
            readTime.processorCase(oldProcCase);
        }
    }
    Pstream::broadcast(readPath, UPstream::worldComm);

    if (!readPath.empty())
    {
        InfoOrPout
            << "Detected additional non-decomposed files in "
            << readPath << endl;

        // readPath: searching is the same for all file handlers. Typical:
        //  <case>/0.1/uniform   (parent dir, decompose mode)
        //  <case>/processor1/0.1/uniform   (redistribute/reconstruct mode)
        //  <case>/processors2/0.1/uniform  ,,
        // writePath:
        //  uncollated : <case>/0.1/uniform (reconstruct mode). Should only
        //               be done by master
        //  uncollated : <case>/processorXXX/0.1/uniform. Should be done by all.
        //  collated   : <case>/processors2/0.1/uniform. Should be done by
        //               local master only.

        const IOobject writeIO
        (
            "uniform",
            writeDb.time().timeName(),
            writeDb
        );

        // Switch to writeHandler
        if (writeHandler)
        {
            auto oldHandler = fileOperation::fileHandler(writeHandler);

            // Check: fileHandler.comm() is size 1 for uncollated
            const label writeComm = fileHandler().comm();

            if (reconstruct)
            {
                const bool oldParRun = UPstream::parRun(false);
                const label oldNumProcs(fileHandler().nProcs());
                const fileName writePath
                (
                    fileHandler().objectPath
                    (
                        writeIO,
                        word::null
                    )
                );
                fileHandler().cp(readPath, writePath);
                const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
                UPstream::parRun(oldParRun);
            }
            else
            {
                const fileName writePath
                (
                    fileHandler().objectPath
                    (
                        writeIO,
                        word::null
                    )
                );

                if (::debug)
                {
                    Pout<< "    readPath   :" << readPath << endl;
                    Pout<< "    writePath  :" << writePath << endl;
                }

                fileHandler().broadcastCopy
                (
                    writeComm,                   // send to all in writeComm
                    UPstream::master(writeComm), // to use ioranks. Check!
                    readPath,
                    writePath
                );
            }
            writeHandler = fileOperation::fileHandler(oldHandler);
        }
    }
}


void printMeshData(const polyMesh& mesh)
{
    // Collect all data on master

    labelListList patchNeiProcNo(Pstream::nProcs());
    labelListList patchSize(Pstream::nProcs());
    const labelList& pPatches = mesh.globalData().processorPatches();
    patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    forAll(pPatches, i)
    {
        const processorPolyPatch& ppp = refCast<const processorPolyPatch>
        (
            mesh.boundaryMesh()[pPatches[i]]
        );
        patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
        patchSize[Pstream::myProcNo()][i] = ppp.size();
    }
    Pstream::gatherList(patchNeiProcNo);
    Pstream::gatherList(patchSize);


    // Print stats

    const globalIndex globalCells(mesh.nCells());
    const globalIndex globalBoundaryFaces(mesh.nBoundaryFaces());

    label maxProcCells = 0;
    label maxProcFaces = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;

    for (const int proci : Pstream::allProcs())
    {
        const label nLocalCells = globalCells.localSize(proci);
        const label nBndFaces = globalBoundaryFaces.localSize(proci);

        InfoOrPout<< nl
            << "Processor " << proci;

        if (!nLocalCells)
        {
            InfoOrPout<< " (empty)" << endl;
            continue;
        }
        else
        {
            InfoOrPout<< nl
                << "    Number of cells = " << nLocalCells << endl;
        }

        label nProcFaces = 0;
        const labelList& nei = patchNeiProcNo[proci];

        forAll(patchNeiProcNo[proci], i)
        {
            InfoOrPout
                << "    Number of faces shared with processor "
                << patchNeiProcNo[proci][i] << " = "
                << patchSize[proci][i] << nl;

            nProcFaces += patchSize[proci][i];
        }

        {
            InfoOrPout
                << "    Number of processor patches = " << nei.size() << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = "
                << nBndFaces-nProcFaces << endl;
        }

        maxProcCells = max(maxProcCells, nLocalCells);
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcFaces = max(maxProcFaces, nProcFaces);
        maxProcPatches = max(maxProcPatches, nei.size());
    }

    // Summary stats

    InfoOrPout
        << nl
        << "Number of processor faces = " << (totProcFaces/2) << nl
        << "Max number of cells = " << maxProcCells;

    if (maxProcCells != globalCells.totalSize())
    {
        scalar avgValue = scalar(globalCells.totalSize())/Pstream::nProcs();

        InfoOrPout
            << " (" << 100.0*(maxProcCells-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    InfoOrPout<< nl;

    InfoOrPout<< "Max number of processor patches = " << maxProcPatches;
    if (totProcPatches)
    {
        scalar avgValue = scalar(totProcPatches)/Pstream::nProcs();

        InfoOrPout
            << " (" << 100.0*(maxProcPatches-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    InfoOrPout<< nl;

    InfoOrPout<< "Max number of faces between processors = " << maxProcFaces;
    if (totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/Pstream::nProcs();

        InfoOrPout
            << " (" << 100.0*(maxProcFaces-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    InfoOrPout<< nl << endl;
}


// Debugging: write volScalarField with decomposition for post processing.
void writeDecomposition
(
    const word& name,
    const fvMesh& mesh,
    const labelUList& decomp
)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    IOListRef<label>
    (
        IOobject
        (
            "cellDecomposition",
            mesh.facesInstance(),  // mesh read from facesInstance
            mesh,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        decomp
    ).write();

    InfoOrPout
        << "Writing wanted cell distribution to volScalarField " << name
        << " for postprocessing purposes." << nl << endl;

    volScalarField procCells
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("", dimless, -1),
        fvPatchFieldBase::zeroGradientType()
    );

    forAll(procCells, celli)
    {
        procCells[celli] = decomp[celli];
    }

    procCells.correctBoundaryConditions();
    procCells.write();
}


void determineDecomposition
(
    refPtr<fileOperation>& readHandler,
    const Time& baseRunTime,
    const fileName& decompDictFile, // optional location for decomposeParDict
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const fileName& proc0CaseName,
    const fvMesh& mesh,
    const bool writeCellDist,

    label& nDestProcs,
    labelList& decomp
)
{
    // Switch to readHandler since decomposition method might do IO
    // (e.g. read decomposeParDict)
    auto oldHandler = fileOperation::fileHandler(readHandler);

    // Read decomposeParDict (on all processors)
    const decompositionModel& method = decompositionModel::New
    (
        mesh,
        decompDictFile
    );

    decompositionMethod& decomposer = method.decomposer();

    if (!decomposer.parallelAware())
    {
        WarningInFunction
            << "You have selected decomposition method \""
            << decomposer.type() << "\n"
            << "    which does not synchronise decomposition across"
               " processor patches.\n"
               "    You might want to select a decomposition method"
               " that is aware of this. Continuing...." << endl;
    }

    Time& tm = const_cast<Time&>(mesh.time());

    const bool oldProcCase = tm.processorCase();
    if (decompose)
    {
        InfoOrPout
            << "Setting caseName to " << baseRunTime.caseName()
            << " to read decomposeParDict" << endl;
        tm.caseName() = baseRunTime.caseName();
        tm.processorCase(false);
    }

    scalarField cellWeights;
    if (method.found("weightField"))
    {
        word weightName = method.get<word>("weightField");

        volScalarField weights
        (
            IOobject
            (
                weightName,
                tm.timeName(),
                mesh,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            mesh
        );
        cellWeights = weights.internalField();
    }

    nDestProcs = decomposer.nDomains();
    decomp = decomposer.decompose(mesh, cellWeights);

    readHandler = fileOperation::fileHandler(oldHandler);


    if (decompose)
    {
        InfoOrPout
            << "Restoring caseName" << endl;
        tm.caseName() = proc0CaseName;
        tm.processorCase(oldProcCase);
    }

    // Dump decomposition to volScalarField
    if (writeCellDist)
    {
        const label oldNumProcs =
            const_cast<fileOperation&>(fileHandler()).nProcs(nDestProcs);

        // Note: on master make sure to write to processor0
        if (decompose)
        {
            if (UPstream::master())
            {
                const bool oldParRun = UPstream::parRun(false);

                InfoOrPout
                    << "Setting caseName to " << baseRunTime.caseName()
                    << " to write undecomposed cellDist" << endl;

                tm.caseName() = baseRunTime.caseName();
                tm.processorCase(false);
                writeDecomposition("cellDist", mesh, decomp);
                InfoOrPout<< "Restoring caseName" << endl;
                tm.caseName() = proc0CaseName;
                tm.processorCase(oldProcCase);

                UPstream::parRun(oldParRun);
            }
        }
        else
        {
            writeDecomposition("cellDist", mesh, decomp);
        }

        const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
    }
}


// Variant of GeometricField::correctBoundaryConditions that only
// evaluates selected patch fields
template<class CoupledPatchType, class GeoField>
void correctCoupledBoundaryConditions(fvMesh& mesh)
{
    for (GeoField& fld : mesh.sorted<GeoField>())
    {
        fld.boundaryFieldRef().template evaluateCoupled<CoupledPatchType>();
    }
}


// Inplace redistribute mesh and any fields
autoPtr<mapDistributePolyMesh> redistributeAndWrite
(
    refPtr<fileOperation>& readHandler,
    refPtr<fileOperation>& writeHandler,
    const Time& baseRunTime,
    const fileName& proc0CaseName,

    // Controls
    const bool doReadFields,
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const bool reconstruct,
    const bool overwrite,

    // Decomposition information
    const label nDestProcs,
    const labelList& decomp,

    // Mesh information
    const boolList& volMeshOnProc,
    const fileName& volMeshInstance,
    fvMesh& mesh
)
{
    Time& runTime = const_cast<Time&>(mesh.time());
    const bool oldProcCase = runTime.processorCase();

    //// Print some statistics
    //Pout<< "Before distribution:" << endl;
    //printMeshData(mesh);


    // Storage of fields

    PtrList<volScalarField> volScalarFields;
    PtrList<volVectorField> volVectorFields;
    PtrList<volSphericalTensorField> volSphereTensorFields;
    PtrList<volSymmTensorField> volSymmTensorFields;
    PtrList<volTensorField> volTensorFields;

    PtrList<surfaceScalarField> surfScalarFields;
    PtrList<surfaceVectorField> surfVectorFields;
    PtrList<surfaceSphericalTensorField> surfSphereTensorFields;
    PtrList<surfaceSymmTensorField> surfSymmTensorFields;
    PtrList<surfaceTensorField> surfTensorFields;

    PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
    PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
    PtrList<DimensionedField<sphericalTensor, volMesh>> dimSphereTensorFields;
    PtrList<DimensionedField<symmTensor, volMesh>> dimSymmTensorFields;
    PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;

    PtrList<pointScalarField> pointScalarFields;
    PtrList<pointVectorField> pointVectorFields;
    PtrList<pointTensorField> pointTensorFields;
    PtrList<pointSphericalTensorField> pointSphTensorFields;
    PtrList<pointSymmTensorField> pointSymmTensorFields;

    // Self-contained pointMesh for reading pointFields
    const pointMesh oldPointMesh(mesh);

    // Track how many (if any) pointFields are read/mapped
    label nPointFields = 0;

    refPtr<fileOperation> noWriteHandler;

    parPointFieldDistributor pointDistributor
    (
        oldPointMesh,   // source mesh
        false,          // savePoints=false (ie, delay until later)
        //false           // Do not write
        noWriteHandler    // Do not write
    );


    if (doReadFields)
    {
        // Create 0 sized mesh to do all the generation of zero sized
        // fields on processors that have zero sized meshes. Note that this is
        // only necessary on master but since polyMesh construction with
        // Pstream::parRun does parallel comms we have to do it on all
        // processors
        autoPtr<fvMeshSubset> subsetterPtr;

        // Missing a volume mesh somewhere?
        if (volMeshOnProc.found(false))
        {
            // A zero-sized mesh with boundaries.
            // This is used to create zero-sized fields.
            subsetterPtr.reset(new fvMeshSubset(mesh, zero{}));
            subsetterPtr().subMesh().init(true);
            subsetterPtr().subMesh().globalData();
            subsetterPtr().subMesh().tetBasePtIs();
            subsetterPtr().subMesh().geometricD();
        }


        // Get original objects (before incrementing time!)
        if (decompose)
        {
            InfoOrPout
                << "Setting caseName to " << baseRunTime.caseName()
                << " to read IOobjects" << endl;
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }

        //IOobjectList objects(mesh, runTime.timeName());
        // Swap to reading fileHandler and read IOobjects
        IOobjectList objects;
        if (readHandler)
        {
            auto oldHandler = fileOperation::fileHandler(readHandler);
            const label oldComm = UPstream::commWorld(fileHandler().comm());

            objects = IOobjectList(mesh, runTime.timeName());
            readHandler = fileOperation::fileHandler(oldHandler);
            UPstream::commWorld(oldComm);
        }

        if (decompose)
        {
            InfoOrPout
                << "Restoring caseName" << endl;
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }

        InfoOrPout
            << "From time " << runTime.timeName()
            << " mesh:" << mesh.objectRegistry::objectRelPath()
            << " have objects:" << objects.names() << endl;

        // We don't want to map the decomposition (mapping already tested when
        // mapping the cell centre field)
        (void)objects.erase("cellDist");


        if (decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }

        // Field reading

        #undef  doFieldReading
        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                volMeshOnProc, readHandler, mesh, subsetterPtr,               \
                objects, Storage                                              \
            );                                                                \
        }

        // volField
        doFieldReading(volScalarFields);
        doFieldReading(volVectorFields);
        doFieldReading(volSphereTensorFields);
        doFieldReading(volSymmTensorFields);
        doFieldReading(volTensorFields);

        // surfaceField
        doFieldReading(surfScalarFields);
        doFieldReading(surfVectorFields);
        doFieldReading(surfSphereTensorFields);
        doFieldReading(surfSymmTensorFields);
        doFieldReading(surfTensorFields);

        // Dimensioned internal fields
        doFieldReading(dimScalarFields);
        doFieldReading(dimVectorFields);
        doFieldReading(dimSphereTensorFields);
        doFieldReading(dimSymmTensorFields);
        doFieldReading(dimTensorFields);

        // pointFields
        nPointFields = 0;

        #undef  doFieldReading
        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                volMeshOnProc, readHandler, oldPointMesh,                     \
                subsetterPtr, objects, Storage,                               \
                true  /* (deregister field) */                                \
            );                                                                \
            nPointFields += Storage.size();                                   \
        }

        doFieldReading(pointScalarFields);
        doFieldReading(pointVectorFields);
        doFieldReading(pointSphTensorFields);
        doFieldReading(pointSymmTensorFields);
        doFieldReading(pointTensorFields);
        #undef doFieldReading


        // Done reading

        if (decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }
    }

    // Save pointMesh information before any topology changes occur!
    if (nPointFields)
    {
        pointDistributor.saveMeshPoints();
    }


    // Read refinement data
    autoPtr<hexRef8Data> refDataPtr;
    {
        // Read refinement data
        if (decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }
        IOobject io
        (
            "dummy",
            volMeshInstance,    //mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobjectOption::LAZY_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );

        if (readHandler)
        {
            auto oldHandler = fileOperation::fileHandler(readHandler);
            const label oldComm = UPstream::commWorld(fileHandler().comm());

            // Read
            refDataPtr.reset(new hexRef8Data(io));

            UPstream::commWorld(oldComm);
            readHandler = fileOperation::fileHandler(oldHandler);
        }
        else
        {
            io.readOpt(IOobjectOption::NO_READ);
            refDataPtr.reset(new hexRef8Data(io));
        }

        if (decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }

        // Make sure all processors have valid data (since only some will
        // read)
        refDataPtr().sync(io);

        // Now we've read refinement data we can remove it
        meshRefinement::removeFiles(mesh);
    }


    // If faMeshesRegistry exists, it is also owned by the polyMesh and will
    // be destroyed by clearGeom() in fvMeshDistribute::distribute()
    //
    // Rescue faMeshesRegistry from destruction by temporarily moving
    // it to be locally owned.
    std::unique_ptr<faMeshesRegistry> faMeshesRegistry_saved
    (
        faMeshesRegistry::Release(mesh)
    );

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do all the distribution of mesh and fields
    autoPtr<mapDistributePolyMesh> distMap = distributor.distribute(decomp);

    // Restore ownership onto the polyMesh
    faMeshesRegistry::Store(std::move(faMeshesRegistry_saved));


    // Print some statistics
    InfoOrPout<< "After distribution:" << endl;
    printMeshData(mesh);

    // Get other side of processor boundaries
    do
    {
        #undef  doCorrectCoupled
        #define doCorrectCoupled(FieldType)  \
        correctCoupledBoundaryConditions<processorFvPatch, FieldType>(mesh);

        doCorrectCoupled(volScalarField);
        doCorrectCoupled(volVectorField);
        doCorrectCoupled(volSphericalTensorField);
        doCorrectCoupled(volSymmTensorField);
        doCorrectCoupled(volTensorField);
        #undef doCorrectCoupled
    }
    while (false);

    // No update surface fields


    // Map pointFields
    if (nPointFields)
    {
        // Construct new pointMesh from distributed mesh
        const pointMesh& newPointMesh = pointMesh::New(mesh);

        pointDistributor.resetTarget(newPointMesh, distMap());

        pointDistributor.distributeAndStore(pointScalarFields);
        pointDistributor.distributeAndStore(pointVectorFields);
        pointDistributor.distributeAndStore(pointSphTensorFields);
        pointDistributor.distributeAndStore(pointSymmTensorFields);
        pointDistributor.distributeAndStore(pointTensorFields);
    }


    // More precision (for points data)
    IOstream::minPrecision(10);


    if (!overwrite)
    {
        ++runTime;
        mesh.setInstance(runTime.timeName());
    }
    else
    {
        mesh.setInstance(volMeshInstance);
    }


    // Register mapDistributePolyMesh for automatic writing...
    IOmapDistributePolyMeshRef distMapRef
    (
        IOobject
        (
            "procAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::AUTO_WRITE
        ),
        distMap()
    );


    if (reconstruct)
    {
        auto oldHandler = fileOperation::fileHandler(writeHandler);

        if (UPstream::master())
        {
            InfoOrPout
                << "Setting caseName to " << baseRunTime.caseName()
                << " to write reconstructed mesh (and fields)." << endl;
            runTime.caseName() = baseRunTime.caseName();
            const bool oldProcCase(runTime.processorCase(false));
            const label oldNumProcs
            (
                const_cast<fileOperation&>(fileHandler()).nProcs(nDestProcs)
            );
            const bool oldParRun = UPstream::parRun(false);

            mesh.write();
            topoSet::removeFiles(mesh);

            // Now we've written all. Reset caseName on master
            InfoOrPout<< "Restoring caseName" << endl;
            UPstream::parRun(oldParRun);
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }

        writeHandler = fileOperation::fileHandler(oldHandler);
    }
    else
    {
        auto oldHandler = fileOperation::fileHandler(writeHandler);

        const label oldNumProcs
        (
            const_cast<fileOperation&>(fileHandler()).nProcs(nDestProcs)
        );
        mesh.write();
        const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);

        writeHandler = fileOperation::fileHandler(oldHandler);

        topoSet::removeFiles(mesh);
    }
    InfoOrPout
        << "Written redistributed mesh to "
        << mesh.facesInstance() << nl << endl;


    if (decompose)
    {
        // Decompose (1 -> N)
        // so {boundary,cell,face,point}ProcAddressing have meaning
        fvMeshTools::writeProcAddressing
        (
            mesh,
            distMap(),
            decompose,
            mesh.facesInstance(),    //oldFacesInstance,
            writeHandler             // to write *ProcAddressing
        );
    }
    else if (reconstruct)
    {
        // Reconstruct (N -> 1)
        // so {boundary,cell,face,point}ProcAddressing have meaning. Make sure
        // to write these to meshes containing the source meshes (i.e. using
        // the read handler)
        fvMeshTools::writeProcAddressing
        (
            mesh,
            distMap(),
            decompose,
            volMeshInstance,    //oldFacesInstance,
            readHandler //writeHandler
        );
    }
    else
    {
        // Redistribute (N -> M)
        // {boundary,cell,face,point}ProcAddressing would be incorrect
        // - can either remove or redistribute previous
        removeProcAddressing(mesh);
    }


    // Refinement data
    if (refDataPtr)
    {
        auto& refData = refDataPtr();

        // Set instance
        IOobject io
        (
            "dummy",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );
        refData.sync(io);


        // Distribute
        refData.distribute(distMap());


        auto oldHandler = fileOperation::fileHandler(writeHandler);

        if (reconstruct)
        {
            if (UPstream::master())
            {
                const bool oldParRun = UPstream::parRun(false);
                const label oldNumProcs
                (
                    const_cast<fileOperation&>(fileHandler()).nProcs(nDestProcs)
                );

                InfoOrPout
                    << "Setting caseName to " << baseRunTime.caseName()
                    << " to write reconstructed refinement data." << endl;
                runTime.caseName() = baseRunTime.caseName();
                const bool oldProcCase(runTime.processorCase(false));

                refData.write();

                // Now we've written all. Reset caseName on master
                InfoOrPout<< "Restoring caseName" << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);

                const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
                UPstream::parRun(oldParRun);
            }
        }
        else
        {
            const label oldNumProcs
            (
                const_cast<fileOperation&>(fileHandler()).nProcs(nDestProcs)
            );
            refData.write();
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
        }

        writeHandler = fileOperation::fileHandler(oldHandler);
    }

    //// Sets. Disabled for now.
    //{
    //    // Read sets
    //    if (decompose)
    //    {
    //        runTime.caseName() = baseRunTime.caseName();
    //        runTime.processorCase(false);
    //    }
    //    IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
    //
    //    PtrList<cellSet> cellSets;
    //    ReadFields(objects, cellSets);
    //
    //    if (decompose)
    //    {
    //        runTime.caseName() = proc0CaseName;
    //        runTime.processorCase(oldProcCase);
    //    }
    //
    //    forAll(cellSets, i)
    //    {
    //        cellSets[i].distribute(distMap());
    //    }
    //
    //    if (reconstruct)
    //    {
    //        if (Pstream::master())
    //        {
    //            Pout<< "Setting caseName to " << baseRunTime.caseName()
    //                << " to write reconstructed refinement data." << endl;
    //            runTime.caseName() = baseRunTime.caseName();
    //            const bool oldProcCase(runTime.processorCase(false));
    //
    //            forAll(cellSets, i)
    //            {
    //                cellSets[i].distribute(distMap());
    //            }
    //
    //            // Now we've written all. Reset caseName on master
    //            Pout<< "Restoring caseName" << endl;
    //            runTime.caseName() = proc0CaseName;
    //            runTime.processorCase(oldProcCase);
    //        }
    //    }
    //    else
    //    {
    //        forAll(cellSets, i)
    //        {
    //            cellSets[i].distribute(distMap());
    //        }
    //    }
    //}


    return distMap;
}


/*---------------------------------------------------------------------------*\
                                     main
\*---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Redistribute decomposed mesh and fields according"
        " to the decomposeParDict settings.\n"
        "Optionally run in decompose/reconstruct mode"
    );

    argList::noFunctionObjects();  // Never use function objects

    // enable -constant ... if someone really wants it
    // enable -zeroTime to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);

    #include "addAllRegionOptions.H"

    #include "addOverwriteOption.H"
    argList::addBoolOption("decompose", "Decompose case");
    argList::addBoolOption("reconstruct", "Reconstruct case");
    argList::addVerboseOption
    (
        "Additional verbosity. (Can be used multiple times for debugging)"
    );
    argList::addDryRunOption
    (
        "Test without writing the decomposition. "
        "Changes -cellDist to only write volScalarField."
    );
    argList::addVerboseOption("Additional verbosity");
    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );
    argList::addBoolOption
    (
        "newTimes",
        "Only reconstruct new times (i.e. that do not exist already)"
    );
    argList::addVerboseOption
    (
        "Additional verbosity. (Can be used multiple times)"
    );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress finiteArea mesh/field handling",
        true  // Advanced option
    );


    //- Disable caching of times/processor dirs etc. Cause massive parallel
    //  problems when e.g decomposing.
    fileOperation::cacheLevel(0);


    // Handle arguments
    // ~~~~~~~~~~~~~~~~
    // (replacement for setRootCase that does not abort)

    argList args(argc, argv);


    // As much as possible avoid synchronised operation. To be looked at more
    // closely for the three scenarios:
    // - decompose - reads on master (and from parent directory) and sends
    //               dictionary to slaves
    // - distribute - reads on potentially a different number of processors
    //                than it writes to
    // - reconstruct - reads parallel, write on master only and to parent
    //                 directory

    // Running data-distributed?
    // (processors cannot see all other processors' files)
    const bool hasDistributedFiles(fileHandler().distributed());


    // Set up loose processorsXXX directory matching (for collated) so e.g.
    // when checking for -latestTime we don't miss anything. Once we know
    // the time, actual number of processors etc we switch back to strict
    // matching.
    fileOperation::nProcsFilter(0);

    // Need this line since we don't include "setRootCase.H"
    #include "foamDlOpenLibs.H"

    const bool reconstruct = args.found("reconstruct");
    const bool writeCellDist = args.found("cellDist");
    const bool dryrun = args.dryRun();
    const bool newTimes = args.found("newTimes");
    const int optVerbose = args.verbose();

    const bool doFiniteArea = !args.found("no-finite-area");
    bool decompose = args.found("decompose");
    bool overwrite = args.found("overwrite");

    // Field restrictions...
    const wordRes selectedFields;
    const wordRes selectedLagrangianFields;


    if (!UPstream::parRun())
    {
        FatalErrorInFunction
            << ": This utility can only be run parallel"
            << exit(FatalError);
    }

    if (decompose && reconstruct)
    {
        FatalErrorInFunction
            << "Cannot specify both -decompose and -reconstruct"
            << exit(FatalError);
    }

    if (optVerbose)
    {
        // Report on output
        faMeshDistributor::verbose_ = 1;
        parPointFieldDistributor::verbose_ = 1;

        if (optVerbose > 1)
        {
            Info<< "Additional debugging enabled" << nl << endl;
            ::debug = 1;
        }
    }

    // Disable NaN setting and floating point error trapping. This is to avoid
    // any issues inside the field redistribution inside fvMeshDistribute
    // which temporarily moves processor faces into existing patches. These
    // will now not have correct values. After all bits have been assembled
    // the processor fields will be restored but until then there might
    // be uninitialised values which might upset some patch field constructors.
    // Workaround by disabling floating point error trapping. TBD: have
    // actual field redistribution instead of subsetting inside
    // fvMeshDistribute.
    Foam::sigFpe::unset(true);


    // File handlers (read/write)
    // ==========================

    // Read handler on processors with a volMesh
    refPtr<fileOperation> volMeshReadHandler;

    // Read handler on processors with an areaMesh
    refPtr<fileOperation> areaMeshReadHandler;

    // Handler for master-only operation (read/writing from/to undecomposed)
    // - only the 'real' master, not io-rank masters
    refPtr<fileOperation> masterOnlyHandler;

    if (UPstream::master(UPstream::worldComm))
    {
        const bool oldParRun = UPstream::parRun(false);

        masterOnlyHandler = fileOperation::NewUncollated();

        UPstream::parRun(oldParRun);
    }

    if (decompose)
    {
        InfoOrPout<< "Decomposing case (like decomposePar)"
            << nl << endl;
    }
    else if (reconstruct)
    {
        InfoOrPout<< "Reconstructing case (like reconstructParMesh)"
            << nl << endl;
    }

    if (decompose || reconstruct)
    {
        // The UPstream::nProcs is either the source or destination procs
        fileOperation::nProcsFilter(UPstream::nProcs());
        InfoOrPout<< "Switching to exact matching for "
            << fileOperation::processorsBaseDir + Foam::name(UPstream::nProcs())
            << " processor directories"
            << nl << endl;
    }
    else
    {
        // Redistribute mode. Accept any processorsXXX naming since we don't
        // know yet what the source/target number of processors is
        fileOperation::nProcsFilter(0);
        InfoOrPout<< "Switching to matching any "
            << fileOperation::processorsBaseDir << " directory" << nl << endl;
    }


    if ((decompose || reconstruct) && !overwrite)
    {
        overwrite = true;

        Warning
            << nl << "    "
            << "Added implicit -overwrite for decompose/reconstruct modes"
            << nl << endl;
    }

    if
    (
        fileHandler().ioRanks().contains(UPstream::myProcNo())
    && !Foam::isDir(args.rootPath())
    )
    {
        //FatalErrorInFunction
        WarningInFunction
            << ": cannot open root directory " << args.rootPath()
            << endl;
            //<< exit(FatalError);
    }


    if (hasDistributedFiles)
    {
        InfoOrPout<< "Detected multiple roots i.e. non-nfs running"
            << nl << endl;

        fileHandler().mkDir(args.globalPath());
    }


    // Check if we have processor directories. Ideally would like to
    // use fileHandler().dirPath here but we don't have runTime yet and
    // want to delay constructing runTime until we've synced all time
    // directories...
    const fileName procDir(fileHandler().filePath(args.path()));
    if (Foam::isDir(procDir))
    {
        if (decompose)
        {
            InfoOrPout<< "Removing existing processor directory:"
                << args.relativePath(procDir) << endl;
            fileHandler().rmDir(procDir);
        }
    }
    else
    {
        // Directory does not exist. If this happens on master -> decompose mode
        if (UPstream::master() && !reconstruct && !decompose)
        {
           decompose = true;
           InfoOrPout
                << "No processor directories; switching on decompose mode"
                << nl << endl;
        }
    }
    // If master changed to decompose mode make sure all nodes know about it
    Pstream::broadcast(decompose);
    if (decompose)
    {
        // The UPstream::nProcs is either the source or destination procs
        fileOperation::nProcsFilter(UPstream::nProcs());
        InfoOrPout
            << "Switching to exact matching for "
            << fileOperation::processorsBaseDir + Foam::name(UPstream::nProcs())
            << " processor directories"
            << nl << endl;
     }



    // If running distributed we have problem of new processors not finding
    // a system/controlDict. However if we switch on the master-only reading
    // the problem becomes that the time directories are differing sizes and
    // e.g. latestTime will pick up a different time (which causes createTime.H
    // to abort). So for now make sure to have master times on all
    // processors
    if (!reconstruct)
    {
        InfoOrPout<< "Creating time directories on all processors"
            << nl << endl;
        createTimeDirs(args.path());
    }

    // Construct time
    // ~~~~~~~~~~~~~~

    // Replace #include "createTime.H" with our own version
    // that has MUST_READ instead of READ_MODIFIED

    Info<< "Create time\n" << endl;
    Time runTime
    (
        Time::controlDictName,
        args,
        false,  // no enableFunctionObjects
        true,   // permit enableLibs
        IOobjectOption::MUST_READ  // Without watching
    );


    refPtr<fileOperation> writeHandler;


    runTime.functionObjects().off();  // Extra safety?
    // Make sure that no runTime checking is done since fileOperation::addWatch
    // etc. does broadcast over world, even if constructed only on a subset
    // of procs
    runTime.runTimeModifiable(false);

    // Save local processor0 casename
    const fileName proc0CaseName = runTime.caseName();
    const bool oldProcCase = runTime.processorCase();


    // Construct undecomposed Time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This will read the same controlDict but might have a different
    // set of times so enforce same times

    if (hasDistributedFiles)
    {
        InfoOrPout<< "Creating time directories for undecomposed Time"
            << " on all processors" << nl << endl;
        createTimeDirs(args.globalPath());
    }


    InfoOrPout<< "Create undecomposed database" << nl << endl;
    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false,  // No function objects
        false,  // No extra controlDict libs
        IOobjectOption::MUST_READ  // Without watching
    );
    // Make sure that no runTime checking is done since fileOperation::addWatch
    // etc. does broadcast over world, even if constructed only on a subset
    // of procs
    baseRunTime.runTimeModifiable(false);


    wordHashSet masterTimeDirSet;
    if (newTimes)
    {
        instantList baseTimeDirs(baseRunTime.times());
        for (const instant& t : baseTimeDirs)
        {
            masterTimeDirSet.insert(t.name());
        }
    }


    // Allow override of decomposeParDict location
    const fileName decompDictFile =
        args.getOrDefault<fileName>("decomposeParDict", "");

    // Get region names
    #include "getAllRegionOptions.H"

    if (regionNames.size() == 1 && regionNames[0] != polyMesh::defaultRegion)
    {
        InfoOrPout<< "Using region: " << regionNames[0] << nl << endl;
    }


    // Demand driven lagrangian mapper
    autoPtr<parLagrangianDistributor> lagrangianDistributorPtr;

    if (reconstruct)
    {
        // use the times list from the master processor
        // and select a subset based on the command-line options
        instantList timeDirs = timeSelector::select(runTime.times(), args);
        Pstream::broadcast(timeDirs);

        if (timeDirs.empty())
        {
            FatalErrorInFunction
                << "No times selected"
                << exit(FatalError);
        }


        // Pass1 : reconstruct mesh and addressing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        InfoOrPout
            << nl
            << "Reconstructing mesh and addressing" << nl << endl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];

            const fileName volMeshSubDir
            (
                polyMesh::meshDir(regionName)
            );
            const fileName areaMeshSubDir
            (
                // Assume single-region area mesh
                faMesh::meshDir(regionName, word::null)
            );

            InfoOrPout
                << nl
                << "Reconstructing mesh:"
                << polyMesh::regionName(regionName) << nl << endl;

            bool areaMeshDetected = false;

            // Loop over all times
            forAll(timeDirs, timeI)
            {
                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                InfoOrPout<< "Time = " << runTime.timeName() << endl << endl;

                // Where meshes are
                fileName volMeshInstance;
                fileName areaMeshInstance;

                volMeshInstance = runTime.findInstance
                (
                    volMeshSubDir,
                    "faces",
                    IOobjectOption::READ_IF_PRESENT
                );

                if (doFiniteArea)
                {
                    areaMeshInstance = runTime.findInstance
                    (
                        areaMeshSubDir,
                        "faceLabels",
                        IOobjectOption::READ_IF_PRESENT
                    );
                }

                Pstream::broadcasts
                (
                    UPstream::worldComm,
                    volMeshInstance,
                    areaMeshInstance
                );


                // Check processors have meshes
                // - check for 'faces' file (polyMesh)
                // - check for 'faceLabels' file (faMesh)
                boolList volMeshOnProc;
                boolList areaMeshOnProc(UPstream::nProcs(), false);

                volMeshOnProc = haveMeshFile
                (
                    runTime,
                    volMeshInstance/volMeshSubDir,
                    "faces"
                );

                // Create handler for reading
                newHandler(volMeshOnProc, volMeshReadHandler);

                if (doFiniteArea)
                {
                    areaMeshOnProc = haveMeshFile
                    (
                        runTime,
                        areaMeshInstance/areaMeshSubDir,
                        "faceLabels"
                    );
                    areaMeshDetected = areaMeshOnProc.found(true);

                    if (areaMeshOnProc == volMeshOnProc)
                    {
                        if (volMeshReadHandler)
                        {
                            // Use same reader for faMesh as for fvMesh
                            areaMeshReadHandler.ref(volMeshReadHandler.ref());
                        }
                    }
                    else
                    {
                        newHandler(areaMeshOnProc, areaMeshReadHandler);
                    }
                }


                // Addressing back to reconstructed mesh as xxxProcAddressing.
                // - all processors have consistent faceProcAddressing
                // - processors without a mesh don't need faceProcAddressing


                // Note: filePath searches up on processors that don't have
                //       processor if instance = constant so explicitly check
                //       found filename.
                bool haveVolAddressing = false;
                if (volMeshOnProc[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    faceCompactIOList faces
                    (
                        IOobject
                        (
                            "faces",
                            volMeshInstance,
                            volMeshSubDir,
                            runTime,
                            IOobjectOption::MUST_READ
                        )
                    );

                    // Check faceProcAddressing
                    labelIOList faceProcAddressing
                    (
                        IOobject
                        (
                            "faceProcAddressing",
                            volMeshInstance,
                            volMeshSubDir,
                            runTime,
                            IOobjectOption::READ_IF_PRESENT
                        )
                    );

                    haveVolAddressing =
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faces.size()
                    );
                }
                else
                {
                    // Have no mesh. Don't need addressing
                    haveVolAddressing = true;
                }

                bool haveAreaAddressing = false;
                if (areaMeshOnProc[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    labelIOList faceLabels
                    (
                        IOobject
                        (
                            "faceLabels",
                            areaMeshInstance,
                            areaMeshSubDir,
                            runTime,
                            IOobjectOption::MUST_READ
                        )
                    );

                    // Check faceProcAddressing
                    labelIOList faceProcAddressing
                    (
                        IOobject
                        (
                            "faceProcAddressing",
                            areaMeshInstance,
                            areaMeshSubDir,
                            runTime,
                            IOobjectOption::READ_IF_PRESENT
                        )
                    );

                    haveAreaAddressing =
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faceLabels.size()
                    );
                }
                else if (areaMeshDetected)
                {
                    // Have no mesh. Don't need addressing
                    haveAreaAddressing = true;
                }


                // Additionally check for master faces being readable. Could
                // do even more checks, e.g. global number of cells same
                // as cellProcAddressing

                bool volMeshHaveUndecomposed = false;
                bool areaMeshHaveUndecomposed = false;

                if (Pstream::master())
                {
                    InfoOrPout
                        << "Checking " << baseRunTime.caseName()
                        << " for undecomposed volume and area meshes..."
                        << endl;

                    const bool oldParRun = Pstream::parRun(false);
                    const label oldNumProcs(fileHandler().nProcs());

                    // Volume
                    {
                        faceCompactIOList facesIO
                        (
                            IOobject
                            (
                                "faces",
                                volMeshInstance,
                                volMeshSubDir,
                                baseRunTime,
                                IOobjectOption::NO_READ
                            ),
                            label(0)
                        );
                        volMeshHaveUndecomposed = facesIO.headerOk();
                    }

                    // Area
                    if (doFiniteArea)
                    {
                        labelIOList labelsIO
                        (
                            IOobject
                            (
                                "faceLabels",
                                areaMeshInstance,
                                areaMeshSubDir,
                                baseRunTime,
                                IOobjectOption::NO_READ
                            )
                        );
                        areaMeshHaveUndecomposed = labelsIO.headerOk();
                    }

                    const_cast<fileOperation&>
                    (
                        fileHandler()
                    ).nProcs(oldNumProcs);
                    Pstream::parRun(oldParRun);  // Restore parallel state
                }

                Pstream::broadcasts
                (
                    UPstream::worldComm,
                    volMeshHaveUndecomposed,
                    areaMeshHaveUndecomposed
                );

                // Report
                {
                    InfoOrPout
                        << "    volume mesh ["
                        << volMeshHaveUndecomposed << "] : "
                        << volMeshInstance << nl
                        << "    area   mesh ["
                        << areaMeshHaveUndecomposed << "] : "
                        << areaMeshInstance << nl
                        << endl;
                }


                if
                (
                    !volMeshHaveUndecomposed
                 || !returnReduceAnd(haveVolAddressing)
                )
                {
                    InfoOrPout
                        << "No undecomposed mesh. Creating from: "
                        << volMeshInstance << endl;

                    if (areaMeshHaveUndecomposed)
                    {
                        areaMeshHaveUndecomposed = false;
                        InfoOrPout
                            << "Also ignore any undecomposed area mesh"
                            << endl;
                    }

                    autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            volMeshInstance,
                            runTime,
                            IOobjectOption::MUST_READ
                        ),
                        volMeshReadHandler
                    );
                    fvMeshTools::setBasicGeometry(volMeshPtr());
                    fvMesh& mesh = volMeshPtr();


                    InfoOrPout<< nl << "Reconstructing mesh" << nl << endl;

                    // Reconstruct (1 processor)
                    const label nDestProcs(1);
                    const labelList finalDecomp(mesh.nCells(), Zero);

                    redistributeAndWrite
                    (
                        volMeshReadHandler,
                        masterOnlyHandler,  //writeHandler,
                        baseRunTime,
                        proc0CaseName,

                        // Controls
                        false,      // do not read fields
                        false,      // do not read undecomposed case on proc0
                        true,       // write redistributed files to proc0
                        overwrite,

                        // Decomposition information
                        nDestProcs,
                        finalDecomp,

                        // For finite-volume
                        volMeshOnProc,
                        volMeshInstance,
                        mesh
                    );
                }


                // Similarly for finiteArea
                // - may or may not have undecomposed mesh
                // - may or may not have decomposed meshes

                if
                (
                    areaMeshOnProc.found(true)  // ie, areaMeshDetected
                 &&
                    (
                        !areaMeshHaveUndecomposed
                     || !returnReduceAnd(haveAreaAddressing)
                    )
                )
                {
                    InfoOrPout
                        << "Loading area mesh from "
                        << areaMeshInstance << endl;

                    InfoOrPout<< "    getting volume mesh support" << endl;

                    autoPtr<fvMesh> baseMeshPtr = fvMeshTools::newMesh
                    (
                        IOobject
                        (
                            regionName,
                            baseRunTime.timeName(),
                            baseRunTime,
                            IOobjectOption::MUST_READ
                        ),
                        true            // read on master only
                    );
                    fvMeshTools::setBasicGeometry(baseMeshPtr());

                    autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            baseMeshPtr().facesInstance(),
                            runTime,
                            IOobjectOption::MUST_READ
                        ),
                        volMeshReadHandler
                    );
                    fvMesh& mesh = volMeshPtr();

                    // Read volume proc addressing back to base mesh
                    autoPtr<mapDistributePolyMesh> distMap
                    (
                        fvMeshTools::readProcAddressing(mesh, baseMeshPtr)
                    );


                    //autoPtr<faMesh> areaMeshPtr =
                    // faMeshTools::loadOrCreateMesh
                    //(
                    //    IOobject
                    //    (
                    //        regionName,
                    //        areaMeshInstance,
                    //        runTime,
                    //        IOobjectOption::MUST_READ
                    //    ),
                    //    mesh,  // <- The referenced polyMesh (from above)
                    //    decompose
                    //);
                    autoPtr<faMesh> areaMeshPtr = faMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            areaMeshInstance,
                            runTime,
                            IOobjectOption::MUST_READ
                        ),
                        mesh,  // <- The referenced polyMesh (from above)
                        areaMeshReadHandler
                    );
                    faMesh& areaMesh = areaMeshPtr();

                    faMeshTools::forceDemandDriven(areaMesh);
                    faMeshTools::unregisterMesh(areaMesh);

                    autoPtr<faMesh> areaBaseMeshPtr;

                    // Reconstruct using polyMesh distribute map
                    mapDistributePolyMesh faDistMap
                    (
                        faMeshDistributor::distribute
                        (
                            areaMesh,
                            distMap(),      // The polyMesh distMap
                            baseMeshPtr(),  // Target polyMesh
                            areaBaseMeshPtr
                        )
                    );

                    faMeshTools::forceDemandDriven(areaBaseMeshPtr());
                    faMeshTools::unregisterMesh(areaBaseMeshPtr());


                    if (Pstream::master())
                    {
                        InfoOrPout
                            << "Setting caseName to " << baseRunTime.caseName()
                            << " to write reconstructed area mesh." << endl;
                        runTime.caseName() = baseRunTime.caseName();
                        const bool oldProcCase(runTime.processorCase(false));
                        const bool oldParRun(Pstream::parRun(false));
                        const label oldNumProcs(fileHandler().nProcs());

                        areaBaseMeshPtr().write();

                        // Now we've written all. Reset caseName on master
                        InfoOrPout<< "Restoring caseName" << endl;
                        const_cast<fileOperation&>
                        (
                            fileHandler()
                        ).nProcs(oldNumProcs);
                        Pstream::parRun(oldParRun);
                        runTime.caseName() = proc0CaseName;
                        runTime.processorCase(oldProcCase);
                    }

                    // Update for the reconstructed procAddressing
                    faMeshTools::writeProcAddressing
                    (
                        areaBaseMeshPtr(),  // Reconstruct location
                        faDistMap,
                        false,              // decompose=false
                        writeHandler,       //writeHandler,
                        areaMeshPtr.get()   // procMesh
                    );
                }
            }

            // Make sure all is finished writing until re-reading in pass2
            // below
            fileHandler().flush();


            // Pass2 : read mesh and addressing and reconstruct fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            InfoOrPout
                << nl
                << "Reconstructing fields" << nl << endl;

            runTime.setTime(timeDirs[0], 0);
            baseRunTime.setTime(timeDirs[0], 0);
            InfoOrPout<< "Time = " << runTime.timeName() << endl << endl;


            // Read undecomposed mesh on master and 'empty' mesh
            // (zero faces, point, cells but valid patches and zones) on slaves.
            // This is a bit of tricky code and hidden inside fvMeshTools for
            // now.
            InfoOrPout<< "Reading undecomposed mesh (on master)" << endl;
            //autoPtr<fvMesh> baseMeshPtr = fvMeshTools::newMesh
            //(
            //    IOobject
            //    (
            //        regionName,
            //        baseRunTime.timeName(),
            //        baseRunTime,
            //        IOobjectOption::MUST_READ
            //    ),
            //    true            // read on master only
            //);
            fileName facesInstance;
            fileName pointsInstance;
            masterMeshInstance
            (
                IOobject
                (
                    regionName,
                    baseRunTime.timeName(),
                    baseRunTime,
                    IOobjectOption::MUST_READ
                ),
                facesInstance,
                pointsInstance
            );

            autoPtr<fvMesh> baseMeshPtr = fvMeshTools::loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    facesInstance,              //baseRunTime.timeName(),
                    baseRunTime,
                    IOobjectOption::MUST_READ
                ),
                masterOnlyHandler               // read on master only
            );

            if (::debug)
            {
                Pout<< "Undecomposed mesh :"
                    << " instance:" << baseMeshPtr().facesInstance()
                    << " nCells:" << baseMeshPtr().nCells()
                    << " nFaces:" << baseMeshPtr().nFaces()
                    << " nPoints:" << baseMeshPtr().nPoints()
                    << endl;
            }

            InfoOrPout<< "Reading local, decomposed mesh" << endl;
            autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    baseMeshPtr().facesInstance(),
                    runTime,
                    IOobjectOption::MUST_READ
                ),
                volMeshReadHandler              // read on fvMesh processors
            );
            fvMesh& mesh = volMeshPtr();


            // Similarly for finiteArea
            autoPtr<faMesh> areaBaseMeshPtr;
            autoPtr<faMesh> areaMeshPtr;
            autoPtr<faMeshDistributor> faDistributor;
            mapDistributePolyMesh areaDistMap;

            if (areaMeshDetected)
            {
                //areaBaseMeshPtr = faMeshTools::newMesh
                //(
                //    IOobject
                //    (
                //        regionName,
                //        baseRunTime.timeName(),
                //        baseRunTime,
                //        IOobjectOption::MUST_READ
                //    ),
                //    baseMeshPtr(),
                //    true            // read on master only
                //);
                areaBaseMeshPtr = faMeshTools::loadOrCreateMesh
                (
                    IOobject
                    (
                        regionName,
                        baseMeshPtr().facesInstance(),
                        baseRunTime,
                        IOobjectOption::MUST_READ
                    ),
                    baseMeshPtr(),
                    masterOnlyHandler
                );

                //areaMeshPtr = faMeshTools::loadOrCreateMesh
                //(
                //    IOobject
                //    (
                //        regionName,
                //        areaBaseMeshPtr().facesInstance(),
                //        runTime,
                //        IOobjectOption::MUST_READ
                //    ),
                //    mesh,
                //    decompose
                //);
                areaMeshPtr = faMeshTools::loadOrCreateMesh
                (
                    IOobject
                    (
                        regionName,
                        areaBaseMeshPtr().facesInstance(),
                        runTime,
                        IOobjectOption::MUST_READ
                    ),
                    mesh,
                    areaMeshReadHandler
                );

                areaDistMap =
                    faMeshTools::readProcAddressing
                    (
                        areaMeshPtr(),
                        areaBaseMeshPtr
                    );

                faMeshTools::forceDemandDriven(areaMeshPtr());

                // Create an appropriate field distributor
                faDistributor.reset
                (
                    new faMeshDistributor
                    (
                        areaMeshPtr(),      // source
                        areaBaseMeshPtr(),  // target
                        areaDistMap,
                        masterOnlyHandler   // only write on master
                    )
                );
                // Report some messages. Tbd.
                faMeshDistributor::verbose_ = 1;
            }


            // Read addressing back to base mesh
            autoPtr<mapDistributePolyMesh> distMap;
            distMap = fvMeshTools::readProcAddressing(mesh, baseMeshPtr);

            // Construct field mapper
            auto fvDistributorPtr =
                autoPtr<parFvFieldDistributor>::New
                (
                    mesh,               // source
                    baseMeshPtr(),      // target
                    distMap(),
                    masterOnlyHandler   // Write on master only
                );

            // Construct point field mapper
            const auto& basePointMesh = pointMesh::New(baseMeshPtr());
            const auto& procPointMesh = pointMesh::New(mesh);

            auto pointFieldDistributorPtr =
                autoPtr<parPointFieldDistributor>::New
                (
                    procPointMesh,   // source
                    basePointMesh,   // target
                    distMap(),
                    false,           // delay
                    //UPstream::master()  // Write reconstructed on master
                    masterOnlyHandler   // Write on master only
                );


            // Since we start from Times[0] and not runTime.timeName() we
            // might overlook point motion in the first timestep
            // (since mesh.readUpdate() below will not be triggered). Instead
            // detect points by hand
            if (mesh.pointsInstance() != mesh.facesInstance())
            {
                InfoOrPout
                    << "    Detected initial mesh motion;"
                    << " reconstructing points" << nl
                    << endl;
                fvDistributorPtr().reconstructPoints();
            }


            // Loop over all times
            forAll(timeDirs, timeI)
            {
                if (newTimes && masterTimeDirSet.found(timeDirs[timeI].name()))
                {
                    InfoOrPout
                        << "Skipping time " << timeDirs[timeI].name()
                        << nl << endl;
                    continue;
                }

                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                InfoOrPout<< "Time = " << runTime.timeName() << endl << endl;


                // Check if any new meshes need to be read.
                polyMesh::readUpdateState procStat = mesh.readUpdate();

                if (procStat == polyMesh::POINTS_MOVED)
                {
                    InfoOrPout
                        << "    Detected mesh motion; reconstructing points"
                        << nl << endl;
                    fvDistributorPtr().reconstructPoints();
                }
                else if
                (
                    procStat == polyMesh::TOPO_CHANGE
                 || procStat == polyMesh::TOPO_PATCH_CHANGE
                )
                {
                    InfoOrPout
                        << "    Detected topology change;"
                        << " reconstructing addressing" << nl << endl;

                    if (baseMeshPtr)
                    {
                        // Cannot do a baseMesh::readUpdate() since not all
                        // processors will have mesh files. So instead just
                        // recreate baseMesh
                        baseMeshPtr.clear();
                        //baseMeshPtr = fvMeshTools::newMesh
                        //(
                        //    IOobject
                        //    (
                        //        regionName,
                        //        baseRunTime.timeName(),
                        //        baseRunTime,
                        //        IOobjectOption::MUST_READ
                        //    ),
                        //    true            // read on master only
                        //);
                        baseMeshPtr = fvMeshTools::loadOrCreateMesh
                        (
                            IOobject
                            (
                                regionName,
                                baseRunTime.timeName(),
                                baseRunTime,
                                IOobjectOption::MUST_READ
                            ),
                            masterOnlyHandler   // read on master only
                        );
                        if (::debug)
                        {
                            Pout<< "Undecomposed mesh :"
                                << " nCells:" << baseMeshPtr().nCells()
                                << " nFaces:" << baseMeshPtr().nFaces()
                                << " nPoints:" << baseMeshPtr().nPoints()
                                << endl;
                        }
                    }

                    // Re-read procAddressing
                    distMap =
                        fvMeshTools::readProcAddressing(mesh, baseMeshPtr);

                    // Reset field mappers

                    fvDistributorPtr.reset
                    (
                        new parFvFieldDistributor
                        (
                            mesh,           // source
                            baseMeshPtr(),  // target
                            distMap(),
                            masterOnlyHandler   // Write on master only
                        )
                    );

                    // Construct point field mapper
                    const auto& basePointMesh = pointMesh::New(baseMeshPtr());
                    const auto& procPointMesh = pointMesh::New(mesh);

                    pointFieldDistributorPtr.reset
                    (
                        new parPointFieldDistributor
                        (
                            procPointMesh,  // source
                            basePointMesh,  // target
                            distMap(),
                            false,          // delay until later
                            masterOnlyHandler   // Write on master only
                        )
                    );

                    lagrangianDistributorPtr.reset();

                    if (areaMeshPtr)
                    {
                        InfoOrPout
                            << "    Discarding finite-area addressing"
                            << " (TODO)" << nl << endl;

                        areaBaseMeshPtr.reset();
                        areaMeshPtr.reset();
                        faDistributor.reset();
                        areaDistMap.clear();
                    }
                }


                // Get list of objects
                IOobjectList objects(mesh, runTime.timeName());

                // Mesh fields (vol, surface, volInternal)
                fvDistributorPtr()
                    .distributeAllFields(objects, selectedFields);

                // pointfields
                // - distribute and write (verbose)
                pointFieldDistributorPtr()
                    .distributeAllFields(objects, selectedFields);


                // Clouds (note: might not be present on all processors)
                reconstructLagrangian
                (
                    lagrangianDistributorPtr,
                    baseMeshPtr(),
                    mesh,
                    distMap(),
                    selectedLagrangianFields
                    //masterOnlyHandler
                );

                if (faDistributor)
                {
                    faDistributor()
                        .distributeAllFields(objects, selectedFields);
                }


                // If there are any "uniform" directories copy them from
                // the master processor

                copyUniform
                (
                    volMeshReadHandler, //masterOnlyHandler,  // readHandler
                    masterOnlyHandler,  // writeHandler

                    true,               // reconstruct
                    false,              // decompose

                    mesh.time().timeName(),
                    word::null,         // optional caseName for reading
                    mesh,
                    baseMeshPtr()
                );
                // Non-region specific. Note: should do outside region loop
                // but would then have to replicate the whole time loop ...
                copyUniform
                (
                    volMeshReadHandler, //masterOnlyHandler,  // readHandler,
                    masterOnlyHandler,  // writeHandler

                    true,               // reconstruct
                    false,              // decompose

                    mesh.time().timeName(),
                    word::null,         // optional caseName for reading
                    mesh.time(),            // runTime
                    baseMeshPtr().time()    // baseRunTime
                );
            }
        }
    }
    else
    {
        // decompose or redistribution mode.
        //  decompose : master : read from parent dir
        //              slave  : dummy mesh
        //  redistribute : all read mesh or dummy mesh

        Time& readRunTime =
        (
            (decompose)
          ? baseRunTime
          : runTime
        );

        // Time coming from processor0 (or undecomposed if no processor0)
        scalar masterTime = timeSelector::selectIfPresent
        (
            readRunTime,
            args
        )[0].value();
        Pstream::broadcast(masterTime);
        InfoOrPout
            << "Setting time to that of master or undecomposed case : "
            << masterTime << endl;
        runTime.setTime(masterTime, 0);
        baseRunTime.setTime(masterTime, 0);




        // Save old time name (since might be incremented)
        const word oldTimeName(runTime.timeName());

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];

            const fileName volMeshSubDir
            (
                polyMesh::meshDir(regionName)
            );
            const fileName areaMeshSubDir
            (
                // Assume single-region area mesh
                faMesh::meshDir(regionName, word::null)
            );

            InfoOrPout
                << nl << nl
                << (decompose ? "Decomposing" : "Redistributing")
                << " mesh:" << polyMesh::regionName(regionName) << nl << endl;


            // Get time instance directory
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // At this point we should be able to read at least a mesh on
            // processor0. Note the changing of the processor0 casename to
            // enforce it to read/write from the undecomposed case

            fileName volMeshMasterInstance;
            fileName areaMeshMasterInstance;

            // Assume to be true
            bool volMeshHaveUndecomposed = true;
            bool areaMeshHaveUndecomposed = doFiniteArea;

            if (Pstream::master())
            {
                if (decompose)
                {
                    InfoOrPout
                        << "Checking undecomposed mesh in case: "
                        << baseRunTime.caseName() << endl;
                    runTime.caseName() = baseRunTime.caseName();
                    runTime.processorCase(false);
                }

                const bool oldParRun = Pstream::parRun(false);
                const label oldNumProcs(fileHandler().nProcs());
                volMeshMasterInstance = readRunTime.findInstance
                (
                    volMeshSubDir,
                    "faces",
                    IOobjectOption::READ_IF_PRESENT
                );

                if (doFiniteArea)
                {
                    areaMeshMasterInstance = readRunTime.findInstance
                    (
                        areaMeshSubDir,
                        "faceLabels",
                        IOobjectOption::READ_IF_PRESENT
                    );

                    // Note: findInstance returns "constant" even if not found,
                    // so recheck now for a false positive.

                    if ("constant" == areaMeshMasterInstance)
                    {
                        const boolList areaMeshOnProc
                        (
                            haveMeshFile
                            (
                                readRunTime,
                                areaMeshMasterInstance/areaMeshSubDir,
                                "faceLabels",
                                false  // verbose=false
                            )
                        );

                        if (areaMeshOnProc.empty() || !areaMeshOnProc[0])
                        {
                            areaMeshHaveUndecomposed = false;
                        }
                    }
                }

                const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
                Pstream::parRun(oldParRun);  // Restore parallel state

                if (decompose)
                {
                    InfoOrPout
                        << "    volume mesh ["
                        << volMeshHaveUndecomposed << "] : "
                        << volMeshMasterInstance << nl
                        << "    area   mesh ["
                        << areaMeshHaveUndecomposed << "] : "
                        << areaMeshMasterInstance << nl
                        << nl << nl;

                    // Restoring caseName
                    InfoOrPout<< "Restoring caseName" << endl;
                    runTime.caseName() = proc0CaseName;
                    runTime.processorCase(oldProcCase);
                }
            }

            Pstream::broadcasts
            (
                UPstream::worldComm,
                volMeshHaveUndecomposed,
                areaMeshHaveUndecomposed,
                volMeshMasterInstance,
                areaMeshMasterInstance
            );

            // Check processors have meshes
            // - check for 'faces' file (polyMesh)
            // - check for 'faceLabels' file (faMesh)
            boolList volMeshOnProc;
            boolList areaMeshOnProc;

            if (decompose)
            {
                // Already determined above that master can read 'faces' file.
                // This avoids doing all the casename setting/restoring again.
                volMeshOnProc.setSize(UPstream::nProcs(), false);
                volMeshOnProc[UPstream::masterNo()] = volMeshHaveUndecomposed;
            }
            else
            {
                // All check if can read 'faces' file
                volMeshOnProc = haveMeshFile
                (
                    runTime,
                    volMeshMasterInstance/volMeshSubDir,
                    "faces"
                );
            }

            // Create handler for reading
            newHandler(volMeshOnProc, volMeshReadHandler);


            // Now we've determined which processors are reading switch back
            // to exact matching of 'processorsXXX' directory names.
            // - this determines the processorsXXX fileNames
            // - the XXX comes from the number of read processors
            // - also adapt the masterOnlyReader (used in copyUniform)

            fileOperation::nProcsFilter
            (
                findIndices(volMeshOnProc, true).size()
            );


            if (doFiniteArea)
            {
                if (decompose)
                {
                    // Already determined above that master can read
                    // 'faceLabels' file.
                    areaMeshOnProc.setSize(UPstream::nProcs(), false);
                    areaMeshOnProc[UPstream::masterNo()] =
                    areaMeshHaveUndecomposed;
                }
                else
                {
                    areaMeshOnProc = haveMeshFile
                    (
                        runTime,
                        areaMeshMasterInstance/areaMeshSubDir,
                        "faceLabels"
                    );
                }

                // Create handler for reading
                if (areaMeshOnProc == volMeshOnProc)
                {
                    if (volMeshReadHandler)
                    {
                        // Use same reader for faMesh as for fvMesh
                        areaMeshReadHandler.ref(volMeshReadHandler.ref());
                    }
                }
                else
                {
                    newHandler(areaMeshOnProc, areaMeshReadHandler);
                }
            }


            // Prior to loadOrCreateMesh, note which meshes already exist
            // for the current file handler.
            // - where mesh would be written if it didn't exist already.
            fileNameList volMeshDir(Pstream::nProcs());
            {
                volMeshDir[Pstream::myProcNo()] =
                (
                    fileHandler().objectPath
                    (
                        IOobject
                        (
                            "faces",
                            volMeshMasterInstance/volMeshSubDir,
                            runTime
                        ),
                        word::null
                    ).path()
                );

                Pstream::allGatherList(volMeshDir);

                if (optVerbose && Pstream::master())
                {
                    Info<< "Per processor faces dirs:" << nl
                        << '(' << nl;

                    for (const int proci : Pstream::allProcs())
                    {
                        Info<< "    "
                            << runTime.relativePath(volMeshDir[proci]);

                        if (!volMeshOnProc[proci])
                        {
                            Info<< " [missing]";
                        }
                        Info<< nl;
                    }
                    Info<< ')' << nl << endl;
                }
            }

            fileNameList areaMeshDir(Pstream::nProcs());
            if (doFiniteArea)
            {
                areaMeshDir[Pstream::myProcNo()] =
                (
                    fileHandler().objectPath
                    (
                        IOobject
                        (
                            "faceLabels",
                            areaMeshMasterInstance/areaMeshSubDir,
                            runTime
                        ),
                        word::null
                    ).path()
                );

                Pstream::allGatherList(areaMeshDir);

                if (optVerbose && Pstream::master())
                {
                    Info<< "Per processor faceLabels dirs:" << nl
                        << '(' << nl;

                    for (const int proci : Pstream::allProcs())
                    {
                        Info<< "    "
                            << runTime.relativePath(areaMeshDir[proci]);

                        if (!areaMeshOnProc[proci])
                        {
                            Info<< " [missing]";
                        }
                        Info<< nl;
                    }
                    Info<< ')' << nl << endl;
                }
            }


            // Load mesh (or create dummy one)
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (decompose)
            {
                InfoOrPout
                    << "Setting caseName to " << baseRunTime.caseName()
                    << " to read undecomposed mesh" << endl;
                runTime.caseName() = baseRunTime.caseName();
                runTime.processorCase(false);
            }

            autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    volMeshMasterInstance,
                    runTime,
                    IOobjectOption::MUST_READ
                ),
                volMeshReadHandler
            );
            fvMesh& mesh = volMeshPtr();


            // Area mesh

            autoPtr<faMesh> areaMeshPtr;

            // Decomposing: must have an undecomposed mesh
            // Redistributing: have any proc mesh
            if
            (
                doFiniteArea
             &&
                (
                    decompose
                  ? areaMeshHaveUndecomposed
                  : areaMeshOnProc.found(true)
                )
            )
            {
                areaMeshPtr = faMeshTools::loadOrCreateMesh
                (
                    IOobject
                    (
                        regionName,
                        areaMeshMasterInstance,
                        runTime,
                        IOobjectOption::MUST_READ
                    ),
                    mesh,  // <- The referenced polyMesh (from above)
                    areaMeshReadHandler
                );

                faMeshTools::forceDemandDriven(*areaMeshPtr);
                faMeshTools::unregisterMesh(*areaMeshPtr);
            }


            if (decompose)
            {
                InfoOrPout<< "Restoring caseName" << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }

            const label nOldCells = mesh.nCells();


            // Determine decomposition
            // ~~~~~~~~~~~~~~~~~~~~~~~

            label nDestProcs;
            labelList finalDecomp;
            determineDecomposition
            (
                volMeshReadHandler,         // how to read decomposeParDict
                baseRunTime,
                decompDictFile,
                decompose,
                proc0CaseName,
                mesh,
                writeCellDist,

                nDestProcs,
                finalDecomp
            );

            if (dryrun)
            {
                continue;
            }

            if (!writeHandler && nDestProcs < fileHandler().nProcs())
            {
                boolList isWriteProc(UPstream::nProcs(), false);
                isWriteProc.slice(0, nDestProcs) = true;
                InfoOrPout
                    << "    dest procs ["
                    << isWriteProc << "]" << nl
                    << endl;
                newHandler(isWriteProc, writeHandler);
            }

            // Area fields first. Read and deregister
            parFaFieldDistributorCache areaFields;
            if (areaMeshPtr)
            {
                areaFields.read
                (
                    baseRunTime,
                    proc0CaseName,
                    decompose,

                    areaMeshOnProc,
                    areaMeshReadHandler,
                    areaMeshMasterInstance,
                    (*areaMeshPtr)
                );
            }


            // Detect lagrangian fields
            if (decompose)
            {
                InfoOrPout
                    << "Setting caseName to " << baseRunTime.caseName()
                    << " to read lagrangian" << endl;
                if (UPstream::master())
                {
                    // Change case name but only on the master - this will
                    // hopefully cause the slaves to not read.
                    runTime.caseName() = baseRunTime.caseName();
                }
                else
                {
                    // Explicitly make sure that casename is not recognised as
                    // a processor case since that has special handling for
                    // caching processor directories etc.
                    runTime.caseName() = "#invalid-name#";
                }
                runTime.processorCase(false);
            }

            // Read lagrangian fields and store on cloud (objectRegistry)
            PtrList<unmappedPassivePositionParticleCloud> clouds
            (
                readLagrangian
                (
                    mesh,
                    selectedLagrangianFields
                )
            );

            if (decompose)
            {
                InfoOrPout<< "Restoring caseName" << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }


            // Load fields, do all distribution (mesh and fields)
            // - but not lagrangian fields; these are done later
            autoPtr<mapDistributePolyMesh> distMap = redistributeAndWrite
            (
                volMeshReadHandler,         // readHandler
                writeHandler,               // writeHandler,
                baseRunTime,
                proc0CaseName,

                // Controls
                true,           // read fields
                decompose,      // decompose, i.e. read from undecomposed case
                false,          // no reconstruction
                overwrite,

                // Decomposition information
                nDestProcs,
                finalDecomp,

                // For finite volume
                volMeshOnProc,
                volMeshMasterInstance,
                mesh
            );

            // Redistribute any clouds
            redistributeLagrangian
            (
                lagrangianDistributorPtr,
                mesh,
                nOldCells,
                distMap(),
                clouds
            );


            // Redistribute area fields

            mapDistributePolyMesh faDistMap;
            autoPtr<faMesh> areaProcMeshPtr;

            if (areaMeshPtr)
            {
                faDistMap = faMeshDistributor::distribute
                (
                    areaMeshPtr(),
                    distMap(),
                    areaProcMeshPtr
                );

                // Force recreation of everything that might vaguely
                // be used by patches:

                faMeshTools::forceDemandDriven(areaProcMeshPtr());


                if (reconstruct)
                {
                    if (Pstream::master())
                    {
                        InfoOrPout
                            << "Setting caseName to " << baseRunTime.caseName()
                            << " to write reconstructed mesh (and fields)."
                            << endl;
                        runTime.caseName() = baseRunTime.caseName();
                        const bool oldProcCase(runTime.processorCase(false));
                        const bool oldParRun = UPstream::parRun(false);
                        const label oldNumProcs(fileHandler().nProcs());

                        areaProcMeshPtr->write();

                        // Now we've written all. Reset caseName on master
                        InfoOrPout<< "Restoring caseName" << endl;
                        const_cast<fileOperation&>
                        (
                            fileHandler()
                        ).nProcs(oldNumProcs);
                        UPstream::parRun(oldParRun);
                        runTime.caseName() = proc0CaseName;
                        runTime.processorCase(oldProcCase);
                    }
                }
                else
                {
                    auto oldHandler = fileOperation::fileHandler(writeHandler);

                    IOmapDistributePolyMeshRef
                    (
                        IOobject
                        (
                            "procAddressing",
                            areaProcMeshPtr->facesInstance(),
                            faMesh::meshSubDir,
                            areaProcMeshPtr->thisDb(),
                            IOobjectOption::NO_READ,
                            IOobjectOption::NO_WRITE,
                            IOobjectOption::NO_REGISTER
                        ),
                        faDistMap
                    ).write();

                    areaProcMeshPtr->write();

                    writeHandler = fileOperation::fileHandler(oldHandler);

                    if (decompose)
                    {
                        faMeshTools::writeProcAddressing
                        (
                            areaProcMeshPtr(),
                            faDistMap,
                            decompose,
                            writeHandler
                        );
                    }
                }

                InfoOrPout
                    << "Written redistributed mesh to "
                    << areaProcMeshPtr->facesInstance() << nl << endl;

                faMeshDistributor distributor
                (
                    areaMeshPtr(),      // source
                    areaProcMeshPtr(),  // target
                    faDistMap,
                    writeHandler
                );

                areaFields.redistributeAndWrite(distributor, true);
            }


            // Get reference to standard write handler
            refPtr<fileOperation> defaultHandler;
            if (writeHandler)
            {
                defaultHandler.ref(writeHandler.ref());
            }
            else
            {
                defaultHandler.ref(const_cast<fileOperation&>(fileHandler()));
            }


            copyUniform
            (
                volMeshReadHandler, // read handler
                defaultHandler,     //TBD: should be all IOranks

                reconstruct,        // reconstruct
                decompose,          // decompose

                oldTimeName,
                (decompose ? baseRunTime.caseName() : proc0CaseName),
                mesh,               // read location is mesh (but oldTimeName)
                mesh                // write location is mesh
            );
        }


        // Get reference to standard write handler
        refPtr<fileOperation> defaultHandler;
        if (writeHandler)
        {
            defaultHandler.ref(writeHandler.ref());
        }
        else
        {
            defaultHandler.ref(const_cast<fileOperation&>(fileHandler()));
        }

        copyUniform
        (
            volMeshReadHandler, // read handler
            defaultHandler,     //TBD: should be all IOranks

            reconstruct,        // reconstruct (=false)
            decompose,          // decompose

            oldTimeName,        // provided read time
            (decompose ? baseRunTime.caseName() : proc0CaseName),
            readRunTime,
            runTime             // writing location
        );
    }


    InfoOrPout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
