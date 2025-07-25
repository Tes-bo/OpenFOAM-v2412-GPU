/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "faMeshTools.H"
#include "faBoundaryMeshEntries.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "fileOperation.H"
#include "BitOps.H"
#include "polyMesh.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::faMeshTools::unregisterMesh(const faMesh& mesh)
{
    auto& obr = const_cast<objectRegistry&>(mesh.thisDb());

    // Checkout by name (casting ambiguity)
    obr.checkOut(faMesh::typeName);
    obr.checkOut("faBoundaryMesh");
    obr.checkOut("faSchemes");
    obr.checkOut("faSolution");
}


void Foam::faMeshTools::forceDemandDriven(faMesh& mesh)
{
    (void)mesh.globalData();

    (void)mesh.Le();
    (void)mesh.magLe();
    (void)mesh.areaCentres();
    (void)mesh.edgeCentres();

    (void)mesh.faceAreaNormals();
    (void)mesh.edgeAreaNormals();
    (void)mesh.pointAreaNormals();
    (void)mesh.faceCurvatures();
    (void)mesh.edgeTransformTensors();

    mesh.syncGeom();
}


Foam::autoPtr<Foam::faMesh>
Foam::faMeshTools::newMesh
(
    const IOobject& io,
    const polyMesh& pMesh,
    const bool masterOnlyReading,
    const bool verbose
)
{
    // The mesh directory (assuming single area region), relative to Time
    const fileName meshSubDir
    (
        faMesh::meshDir(pMesh, word::null)
    );


    fileName facesInstance;

    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (UPstream::master())
    {
        const bool oldParRun = UPstream::parRun(false);

        facesInstance = io.time().findInstance
        (
            meshSubDir,
            "faceLabels",
            IOobject::MUST_READ
        );

        patchEntries = faBoundaryMeshEntries
        (
            IOobject
            (
                "faBoundary",
                facesInstance,
                meshSubDir,
                io.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        UPstream::parRun(oldParRun);
    }

    // Broadcast information to all
    Pstream::broadcasts
    (
        UPstream::worldComm,
        patchEntries,
        facesInstance
    );


    // Dummy meshes
    // ~~~~~~~~~~~~

    // Fake read-if-present behaviour to obtain the faceLabels

    IOobject meshIO(io);
    meshIO.instance() = facesInstance;
    meshIO.readOpt(IOobject::READ_IF_PRESENT);

    // For mesh components (faceLabels, ...)
    IOobject cmptIO(io.time(), "faceLabels", meshSubDir);
    cmptIO.instance() = facesInstance;
    cmptIO.readOpt(IOobject::MUST_READ);
    cmptIO.writeOpt(IOobject::NO_WRITE);
    cmptIO.registerObject(IOobject::NO_REGISTER);


    // Check who has a mesh

    const fileName meshDir = io.time().path()/facesInstance/meshSubDir;
    bool haveMesh = isDir(meshDir);
    if (masterOnlyReading && !UPstream::master())
    {
        haveMesh = false;
        meshIO.readOpt(IOobject::NO_READ);
    }

    if (!haveMesh)
    {
        cmptIO.readOpt(IOobject::NO_READ);
    }


    // Read mesh
    // ~~~~~~~~~
    // Now all processors use supplied points,faces etc
    // Note: solution, schemes are also using the supplied IOobject so
    //       on slave will be NO_READ, on master READ_IF_PRESENT. This will
    //       conflict with e.g. timeStampMaster reading so switch off.

    const auto oldCheckType = IOobject::fileModificationChecking;
    IOobject::fileModificationChecking = IOobject::timeStamp;


    // faceLabels
    cmptIO.rename("faceLabels");
    labelIOList faceLabels(cmptIO);


    auto meshPtr = autoPtr<faMesh>::New
    (
        pMesh,
        std::move(faceLabels),
        meshIO
    );
    auto& mesh = *meshPtr;

    IOobject::fileModificationChecking = oldCheckType;


    // Some processors without patches? - add patches

    if (returnReduceOr(mesh.boundary().empty()))
    {
        // Use patchEntries, which were read on master and broadcast

        faPatchList patches(patchEntries.size());
        label nPatches = 0;

        const bool isEmptyMesh = (mesh.faceLabels().empty());

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
            )
            {
                // Stop at the first processor patch.
                // - logic will not work with inter-mixed proc-patches anyhow
                break;
            }
            else
            {
                dictionary patchDict(e.dict());

                if (isEmptyMesh)
                {
                    patchDict.set("edgeLabels", labelList());
                }

                patches.set
                (
                    patchi,
                    faPatch::New
                    (
                        name,
                        patchDict,
                        nPatches++,
                        mesh.boundary()
                    )
                );
            }
        }

        patches.resize(nPatches);
        mesh.addFaPatches(patches, false);  // No parallel comms
    }

    // Recreate basic geometry, globalMeshData etc.
    mesh.init(false);
    (void)mesh.globalData();

    return meshPtr;
}


Foam::autoPtr<Foam::faMesh>
Foam::faMeshTools::loadOrCreateMeshImpl
(
    const IOobject& io,
    refPtr<fileOperation>* readHandlerPtr,  // Can be nullptr
    const polyMesh& pMesh,
    const bool decompose,
    const bool verbose
)
{
    // The mesh directory (assuming single area region), relative to Time
    const fileName meshSubDir
    (
        faMesh::meshDir(pMesh, word::null)
    );


    // Patch types
    // ~~~~~~~~~~~
    // Read and scatter master patches (without reading master mesh!)

    PtrList<entry> patchEntries;
    if (UPstream::master())
    {
        const bool oldParRun = UPstream::parRun(false);
        const label oldNumProcs = fileHandler().nProcs();
        const int oldCache = fileOperation::cacheLevel(0);

        patchEntries = faBoundaryMeshEntries
        (
            IOobject
            (
                "faBoundary",
                io.instance(),
                meshSubDir,
                io.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        fileOperation::cacheLevel(oldCache);
        if (oldParRun)
        {
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
        }
        UPstream::parRun(oldParRun);
    }

    // Broadcast: send patches to all
    Pstream::broadcast(patchEntries, UPstream::worldComm);


    // Check who has or needs a mesh.
    bool haveLocalMesh = false;

    if (readHandlerPtr)
    {
        // Non-null reference when a mesh exists on given processor
        haveLocalMesh = (*readHandlerPtr).good();
    }
    else
    {
        // No file handler.
        // For 'decompose', only need mesh on master.
        // Otherwise check for presence of the "faceLabels" file

        haveLocalMesh =
        (
            decompose
          ? UPstream::master()
          : fileHandler().isFile
            (
                fileHandler().filePath
                (
                    io.time().path()/io.instance()/meshSubDir/"faceLabels"
                )
            )
        );
    }


    // Globally consistent information about who has a mesh
    boolList haveMesh
    (
        UPstream::allGatherValues<bool>(haveLocalMesh)
    );


    autoPtr<faMesh> meshPtr;

    if (!haveLocalMesh)
    {
        // No local mesh - need to synthesize one

        const bool oldParRun = UPstream::parRun(false);
        const label oldNumProcs = fileHandler().nProcs();
        const int oldCache = fileOperation::cacheLevel(0);

        // Create dummy mesh - on procs that don't already have a mesh
        meshPtr.reset
        (
            new faMesh
            (
                pMesh,
                labelList(),  // Similar to Foam::zero{}
                IOobject(io, IOobject::NO_READ, IOobject::AUTO_WRITE)
            )
        );
        faMesh& mesh = *meshPtr;

        // Add patches
        faPatchList patches(patchEntries.size());
        label nPatches = 0;

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
            )
            {
                // Stop at the first processor patch.
                // - logic will not work with inter-mixed proc-patches anyhow
                break;
            }
            else
            {
                dictionary patchDict(e.dict());
                patchDict.set("edgeLabels", labelList());

                patches.set
                (
                    patchi,
                    faPatch::New
                    (
                        name,
                        patchDict,
                        nPatches++,
                        mesh.boundary()
                    )
                );
            }
        }
        patches.resize(nPatches);
        mesh.addFaPatches(patches, false);  // No parallel comms

        if (!readHandlerPtr)
        {
            // The 'old' way of doing things.
            // Write the dummy mesh to disk for subsequent re-reading.
            //
            // This is not particularly elegant.

            // Bad hack, but the underlying polyMesh is NO_WRITE
            // so it does not create the faMesh subDir for us...
            Foam::mkDir(mesh.boundary().path());

            //Pout<< "Writing dummy mesh to " << mesh.boundary().path() << nl;
            mesh.write();

            // Discard - it will be re-read later
            meshPtr.reset(nullptr);
        }

        fileOperation::cacheLevel(oldCache);
        if (oldParRun)
        {
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
        }
        UPstream::parRun(oldParRun);  // Restore parallel state
    }
    else if (readHandlerPtr && haveLocalMesh)
    {
        const labelList meshProcIds(BitOps::sortedToc(haveMesh));

        UPstream::communicator newCommunicator;
        const label oldWorldComm = UPstream::commWorld();

        auto& readHandler = *readHandlerPtr;
        auto oldHandler = fileOperation::fileHandler(readHandler);

        // With IO ranks the communicator of the fileOperation will
        // only include the ranks for the current IO rank.
        // Instead allocate a new communicator for everyone with a mesh

        const auto& handlerProcIds = UPstream::procID(fileHandler().comm());

        // Comparing global ranks in the communicator.
        // Use std::equal for the List<label> vs List<int> comparison

        if
        (
            meshProcIds.size() == handlerProcIds.size()
         && std::equal
            (
                meshProcIds.cbegin(),
                meshProcIds.cend(),
                handlerProcIds.cbegin()
            )
        )
        {
            // Can use the handler communicator as is.
            UPstream::commWorld(fileHandler().comm());
        }
        else if
        (
            UPstream::nProcs(fileHandler().comm())
         != UPstream::nProcs(UPstream::worldComm)
        )
        {
            // Need a new communicator for the fileHandler.

            // Warning: MS-MPI currently uses MPI_Comm_create() instead of
            // MPI_Comm_create_group() so it will block here!

            newCommunicator.reset(UPstream::worldComm, meshProcIds);
            UPstream::commWorld(newCommunicator.comm());
        }

        // Load but do not initialise
        meshPtr = autoPtr<faMesh>::New(pMesh, false);

        readHandler = fileOperation::fileHandler(oldHandler);
        UPstream::commWorld(oldWorldComm);

        // Reset mesh communicator to the real world comm
        meshPtr().comm() = UPstream::commWorld();
    }


    if (!meshPtr)
    {
        // Using the 'old' way of doing things (writing to disk and re-reading).

        // Read mesh from disk
        //
        // Now all processors have a (possibly zero size) mesh so can
        // read in parallel

        /// Pout<< "Reading area mesh from " << io.objectRelPath() << endl;
        // Load but do not initialise
        meshPtr = autoPtr<faMesh>::New(pMesh, false);
    }

    faMesh& mesh = meshPtr();

    // Check patches
    // ~~~~~~~~~~~~~

    #if 0
    if (!UPstream::master() && haveLocalMesh)
    {
        // Check master names against mine

        const faBoundaryMesh& patches = mesh.boundary();

        forAll(patchEntries, patchi)
        {
            const entry& e = patchEntries[patchi];
            const word type(e.dict().get<word>("type"));
            const word& name = e.keyword();

            if
            (
                type == processorFaPatch::typeName
            )
            {
                break;
            }

            if (patchi >= patches.size())
            {
                FatalErrorInFunction
                    << "Non-processor patches not synchronised." << endl
                    << "Processor " << UPstream::myProcNo()
                    << " has only " << patches.size()
                    << " patches, master has "
                    << patchi
                    << exit(FatalError);
            }

            if
            (
                type != patches[patchi].type()
             || name != patches[patchi].name()
            )
            {
                FatalErrorInFunction
                    << "Non-processor patches not synchronised." << endl
                    << "Master patch " << patchi
                    << " name:" << type
                    << " type:" << type << endl
                    << "Processor " << UPstream::myProcNo()
                    << " patch " << patchi
                    << " has name:" << patches[patchi].name()
                    << " type:" << patches[patchi].type()
                    << exit(FatalError);
            }
        }
    }
    #endif


    // Recreate basic geometry, globalMeshData etc.
    mesh.init(false);
    (void)mesh.globalData();

    /// #if 0
    /// faMeshTools::forceDemandDriven(mesh);
    /// faMeshTools::unregisterMesh(mesh);
    /// #endif

    // Do some checks.

    // Check if the boundary definition is unique
    // and processor patches are correct
    mesh.boundary().checkDefinition(verbose);
    mesh.boundary().checkParallelSync(verbose);

    return meshPtr;
}


Foam::autoPtr<Foam::faMesh>
Foam::faMeshTools::loadOrCreateMesh
(
    const IOobject& io,
    const polyMesh& pMesh,
    const bool decompose,
    const bool verbose
)
{
    return faMeshTools::loadOrCreateMeshImpl
    (
        io,
        nullptr,  // fileOperation (ignore)
        pMesh,
        decompose,
        verbose
    );
}


Foam::autoPtr<Foam::faMesh>
Foam::faMeshTools::loadOrCreateMesh
(
    const IOobject& io,
    const polyMesh& pMesh,
    refPtr<fileOperation>& readHandler,
    const bool verbose
)
{
    return faMeshTools::loadOrCreateMeshImpl
    (
        io,
       &readHandler,
        pMesh,
        false,  // decompose (ignored)
        verbose
    );
}


// ************************************************************************* //
