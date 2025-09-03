/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    reconstructParMesh

Group
    grpParallelUtilities

Description
    Reconstructs a mesh using geometric information only.

    Writes point/face/cell procAddressing so afterwards reconstructPar can be
    used to reconstruct fields.

Usage
    \b reconstructParMesh [OPTION]

    Options:
      - \par -fullMatch
        Does geometric matching on all boundary faces. Assumes no point
        ordering

      - \par -procMatch
        Assumes processor patches already in face order but not point order.
        This is the pre v2106 default behaviour but might be removed if the new
        topological method works well

      - \par -mergeTol \<tol\>
        Specifies non-default merge tolerance (fraction of mesh bounding box)
        for above options

    The default is to assume all processor boundaries are correctly ordered
    (both faces and points) in which case no merge tolerance is needed.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "IOobjectList.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapAddedPolyMesh.H"
#include "polyMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMeshAdder.H"
#include "polyTopoChange.H"
#include "topoSet.H"
#include "regionProperties.H"
#include "fvMeshTools.H"

#include "faMeshReconstructor.H"
#include "ignoreFaPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
// usually meshes get written with limited precision (6 digits)
static const scalar defaultMergeTol = 1e-7;


// Determine which faces are coupled. Uses geometric merge distance.
// Looks either at all boundaryFaces (fullMatch) or only at the
// procBoundaries for proci. Assumes that masterMesh contains already merged
// all the processors < proci.
autoPtr<faceCoupleInfo> determineCoupledFaces
(
    const bool fullMatch,
    const label masterMeshProcStart,
    const label masterMeshProcEnd,
    const polyMesh& masterMesh,
    const label meshToAddProcStart,
    const label meshToAddProcEnd,
    const polyMesh& meshToAdd,
    const scalar mergeDist
)
{
    if (fullMatch || masterMesh.nCells() == 0)
    {
        return autoPtr<faceCoupleInfo>::New
        (
            masterMesh,
            meshToAdd,
            mergeDist,      // Absolute merging distance
            true            // Matching faces identical
        );
    }
    else
    {
        // Pick up all patches on masterMesh ending in "toDDD" where DDD is
        // the processor number proci.

        const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();


        DynamicList<label> masterFaces
        (
            masterMesh.nFaces()
          - masterMesh.nInternalFaces()
        );


        forAll(masterPatches, patchi)
        {
            const polyPatch& pp = masterPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                for
                (
                    label proci=meshToAddProcStart;
                    proci<meshToAddProcEnd;
                    proci++
                )
                {
                    const string toProcString("to" + name(proci));
                    if (
                        pp.name().rfind(toProcString)
                     == (pp.name().size()-toProcString.size())
                    )
                    {
                        label meshFacei = pp.start();
                        forAll(pp, i)
                        {
                            masterFaces.append(meshFacei++);
                        }
                        break;
                    }
                }

            }
        }
        masterFaces.shrink();


        // Pick up all patches on meshToAdd ending in "procBoundaryDDDtoYYY"
        // where DDD is the processor number proci and YYY is < proci.

        const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

        DynamicList<label> addFaces
        (
            meshToAdd.nFaces()
          - meshToAdd.nInternalFaces()
        );

        forAll(addPatches, patchi)
        {
            const polyPatch& pp = addPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                bool isConnected = false;

                for
                (
                    label mergedProci=masterMeshProcStart;
                    !isConnected && (mergedProci < masterMeshProcEnd);
                    mergedProci++
                )
                {
                    for
                    (
                        label proci = meshToAddProcStart;
                        proci < meshToAddProcEnd;
                        proci++
                    )
                    {
                        const word fromProcString
                        (
                            processorPolyPatch::newName(proci, mergedProci)
                        );

                        if (pp.name() == fromProcString)
                        {
                            isConnected = true;
                            break;
                        }
                    }
                }

                if (isConnected)
                {
                    label meshFacei = pp.start();
                    forAll(pp, i)
                    {
                        addFaces.append(meshFacei++);
                    }
                }
            }
        }
        addFaces.shrink();

        return autoPtr<faceCoupleInfo>::New
        (
            masterMesh,
            masterFaces,
            meshToAdd,
            addFaces,
            mergeDist,      // Absolute merging distance
            true,           // Matching faces identical?
            false,          // If perfect match are faces already ordered
                            // (e.g. processor patches)
            false           // are faces each on separate patch?
        );
    }
}


autoPtr<mapPolyMesh> mergeSharedPoints
(
    const scalar mergeDist,
    polyMesh& mesh,
    labelListList& pointProcAddressing
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.
    Map<label> pointToMaster
    (
        fvMeshAdder::findSharedPoints
        (
            mesh,
            mergeDist
        )
    );

    Info<< "mergeSharedPoints : detected " << pointToMaster.size()
        << " points that are to be merged." << endl;

    if (returnReduceAnd(pointToMaster.empty()))
    {
        return nullptr;
    }

    polyTopoChange meshMod(mesh);

    fvMeshAdder::mergePoints(mesh, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields. No inflation, parallel sync.
    mesh.updateMesh(map());

    // pointProcAddressing give indices into the master mesh so adapt them
    // for changed point numbering.

    // Adapt constructMaps for merged points.
    forAll(pointProcAddressing, proci)
    {
        labelList& constructMap = pointProcAddressing[proci];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            // New label of point after changeMesh.
            label newPointi = map().reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }

    return map;
}


boundBox procBounds
(
    const PtrList<Time>& databases,
    const word& regionName
)
{
    boundBox bb;

    for (const Time& procDb : databases)
    {
        fileName pointsInstance
        (
            procDb.findInstance(polyMesh::meshDir(regionName), "points")
        );

        if (pointsInstance != procDb.timeName())
        {
            FatalErrorInFunction
                << "Your time was specified as " << procDb.timeName()
                << " but there is no polyMesh/points in that time." << nl
                << "(points file at " << pointsInstance << ')' << nl
                << "Please rerun with the correct time specified"
                << " (through the -constant, -time or -latestTime "
                << "(at your option)."
                << endl << exit(FatalError);
        }

        Info<< "Reading points from "
            << procDb.caseName()
            << " for time = " << procDb.timeName()
            << nl << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                pointsInstance,
                polyMesh::meshDir(regionName),
                procDb,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        bb.add(points);
    }

    return bb;
}


void writeDistribution
(
    Time& runTime,
    const fvMesh& masterMesh,
    const labelListList& cellProcAddressing
)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            masterMesh.facesInstance(),
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        masterMesh.nCells()
    );

    forAll(cellProcAddressing, proci)
    {
        const labelList& pCells = cellProcAddressing[proci];
        labelUIndList(cellDecomposition, pCells) = proci;
    }

    cellDecomposition.write();

    Info<< nl << "Wrote decomposition to "
        << cellDecomposition.objectRelPath()
        << " for use in manual decomposition." << endl;

    // Write as volScalarField for postprocessing.
    // Change time to 0 if was 'constant'
    {
        const scalar oldTime = runTime.value();
        const label oldIndex = runTime.timeIndex();
        if (runTime.timeName() == runTime.constant() && oldIndex == 0)
        {
            runTime.setTime(0, oldIndex+1);
        }

        volScalarField cellDist
        (
            IOobject
            (
                "cellDist",
                runTime.timeName(),
                masterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            masterMesh,
            dimensionedScalar(word::null, dimless, -1),
            fvPatchFieldBase::zeroGradientType()
        );

        forAll(cellDecomposition, celli)
        {
            cellDist[celli] = cellDecomposition[celli];
        }

        cellDist.correctBoundaryConditions();
        cellDist.write();

        Info<< nl << "Wrote decomposition to "
            << cellDist.objectRelPath()
            << " (volScalarField) for visualization."
            << endl;

        // Restore time
        runTime.setTime(oldTime, oldIndex);
    }
}


void writeMesh
(
    const fvMesh& mesh,
    const labelListList& cellProcAddressing
)
{
    const fileName outputDir
    (
        mesh.time().path()
      / mesh.time().timeName()
      / mesh.regionName()
      / polyMesh::meshSubDir
    );

    Info<< nl << "Writing merged mesh to "
        << mesh.time().relativePath(outputDir) << nl << endl;

    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
    topoSet::removeFiles(mesh);
}


void writeMaps
(
    const label masterInternalFaces,
    const labelUList& masterOwner,
    const polyMesh& procMesh,
    const labelUList& cellProcAddressing,
    const labelUList& faceProcAddressing,
    const labelUList& pointProcAddressing,
    const labelUList& boundProcAddressing
)
{
    const fileName outputDir
    (
        procMesh.time().caseName()
      / procMesh.facesInstance()
      / procMesh.regionName()
      / polyMesh::meshSubDir
    );

    IOobject ioAddr
    (
        "procAddressing",
        procMesh.facesInstance(),
        polyMesh::meshSubDir,
        procMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );


    Info<< "Writing addressing : " << outputDir << nl;

    // From processor point to reconstructed mesh point

    Info<< "    pointProcAddressing" << endl;
    ioAddr.rename("pointProcAddressing");
    IOListRef<label>(ioAddr, pointProcAddressing).write();

    // From processor face to reconstructed mesh face
    Info<< "    faceProcAddressing" << endl;
    ioAddr.rename("faceProcAddressing");
    labelIOList faceProcAddr(ioAddr, faceProcAddressing);

    // Now add turning index to faceProcAddressing.
    // See reconstructPar for meaning of turning index.
    forAll(faceProcAddr, procFacei)
    {
        const label masterFacei = faceProcAddr[procFacei];

        if
        (
           !procMesh.isInternalFace(procFacei)
         && masterFacei < masterInternalFaces
        )
        {
            // proc face is now external but used to be internal face.
            // Check if we have owner or neighbour.

            label procOwn = procMesh.faceOwner()[procFacei];
            label masterOwn = masterOwner[masterFacei];

            if (cellProcAddressing[procOwn] == masterOwn)
            {
                // No turning. Offset by 1.
                faceProcAddr[procFacei]++;
            }
            else
            {
                // Turned face.
                faceProcAddr[procFacei] = -1 - faceProcAddr[procFacei];
            }
        }
        else
        {
            // No turning. Offset by 1.
            faceProcAddr[procFacei]++;
        }
    }

    faceProcAddr.write();


    // From processor cell to reconstructed mesh cell
    Info<< "    cellProcAddressing" << endl;
    ioAddr.rename("cellProcAddressing");
    IOListRef<label>(ioAddr, cellProcAddressing).write();


    // From processor patch to reconstructed mesh patch
    Info<< "    boundaryProcAddressing" << endl;
    ioAddr.rename("boundaryProcAddressing");
    IOListRef<label>(ioAddr, boundProcAddressing).write();

    Info<< endl;
}


void printWarning()
{
    Info<<
"Merge individual processor meshes back into one master mesh.\n"
"Use if the original master mesh has been deleted or the processor meshes\n"
"have been modified (topology change).\n"
"This tool will write the resulting mesh to a new time step and construct\n"
"xxxxProcAddressing files in the processor meshes so reconstructPar can be\n"
"used to regenerate the fields on the master mesh.\n\n"
"Not well tested & use at your own risk!\n\n";
}


// Used to determine the correspondence between the edges.
// See faMeshReconstructor::calcAddressing
void determineFaEdgeMapping
(
    const uindirectPrimitivePatch& onePatch,// reconstructed faMesh patch
    const PtrList<faMesh>& procFaMeshes,    // individual faMeshes
    const labelListList& pointProcAddressing,   // procPolyMesh to reconstructed

    labelListList& faEdgeProcAddressing
)
{
    // Determines ordering from faMesh to onePatch (= recsontructed faMesh)

    // From two polyMesh points to patch-edge label
    EdgeMap<label> pointsToOnePatchEdge(onePatch.nEdges());
    {
        const edgeList& edges = onePatch.edges();
        const labelList& mp = onePatch.meshPoints();

        forAll(edges, edgei)
        {
            const edge meshE(mp, edges[edgei]);
            pointsToOnePatchEdge.insert(meshE, edgei);
        }
    }

    // Now we can create the mapping from patch-edge on processor
    // to patch-edge on onePatch

    faEdgeProcAddressing.resize_nocopy(procFaMeshes.size());
    forAll(procFaMeshes, proci)
    {
        const auto& procPatch = procFaMeshes[proci].patch();
        const auto& edges = procPatch.edges();
        const auto& mp = procPatch.meshPoints();
        const auto& ppAddressing = pointProcAddressing[proci];

        labelList& edgeProcAddr = faEdgeProcAddressing[proci];
        edgeProcAddr.resize_nocopy(edges.size());

        label edgei = 0;
        for
        (
            ;
            edgei < procPatch.nEdges(); //procPatch.nInternalEdges();
            edgei++
        )
        {
            const edge meshE(mp, edges[edgei]);
            const edge onePatchE(ppAddressing, meshE);

            edgeProcAddr[edgei] = pointsToOnePatchEdge[onePatchE];
        }
    }
}


void sortFaEdgeMapping
(
    const uindirectPrimitivePatch& onePatch,// reconstructed faMesh patch
    const PtrList<faMesh>& procFaMeshes,    // individual faMeshes
    const labelListList& pointProcAddressing,   // procPolyMesh to reconstructed

    labelListList& faEdgeProcAddressing,    // per proc the map to the master
    labelListList& singlePatchEdgeLabels    // per patch the map to the master
)
{
    // From faMeshReconstructor.C - edge shuffling on patches

    Map<label> remapGlobal(2*onePatch.nEdges());
    for (label edgei = 0; edgei < onePatch.nInternalEdges(); ++edgei)
    {
        remapGlobal.insert(edgei, remapGlobal.size());
    }


    const faMesh& procMesh0 = procFaMeshes[0];
    const label nNonProcessor = procMesh0.boundary().nNonProcessor();

    singlePatchEdgeLabels.resize_nocopy(nNonProcessor);

    forAll(singlePatchEdgeLabels, patchi)
    {
        if (isA<ignoreFaPatch>(procMesh0.boundary()[patchi]))
        {
            // These are not real edges
            continue;
        }

        forAll(procFaMeshes, proci)
        {
            const faMesh& procMesh = procFaMeshes[proci];
            const faPatch& fap = procMesh.boundary()[patchi];

            labelList patchEdgeLabels(fap.edgeLabels());

            // Renumber from local edges to global edges (natural order)
            for (label& edgeId : patchEdgeLabels)
            {
                edgeId = faEdgeProcAddressing[proci][edgeId];
            }

            // Combine from all processors
            singlePatchEdgeLabels[patchi].append(patchEdgeLabels);
        }

        // Sorted order will be the original non-decomposed order
        Foam::sort(singlePatchEdgeLabels[patchi]);

        for (const label sortedEdgei : singlePatchEdgeLabels[patchi])
        {
            remapGlobal.insert(sortedEdgei, remapGlobal.size());
        }
    }

    {
        // Use the map to rewrite the local faEdgeProcAddressing

        labelListList newEdgeProcAddr(faEdgeProcAddressing);

        forAll(procFaMeshes, proci)
        {
            const faMesh& procMesh = procFaMeshes[proci];

            label edgei = procMesh.nInternalEdges();

            for (const faPatch& fap : procMesh.boundary())
            {
                for (const label patchEdgei : fap.edgeLabels())
                {
                    const label globalEdgei =
                        faEdgeProcAddressing[proci][patchEdgei];

                    const auto fnd = remapGlobal.cfind(globalEdgei);
                    if (fnd.good())
                    {
                        newEdgeProcAddr[proci][edgei] = fnd.val();
                        ++edgei;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Failed to find edge " << globalEdgei
                            << " this indicates a programming error" << nl
                            << exit(FatalError);
                    }
                }
            }
        }
        faEdgeProcAddressing = std::move(newEdgeProcAddr);
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reconstruct a mesh using geometric/topological information only"
    );

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true); // constant(true), zero(true)

    argList::noParallel();

    argList::addVerboseOption();
    argList::addBoolOption
    (
        "addressing-only",
        "Create procAddressing only without overwriting the mesh"
    );
    argList::addOption
    (
        "mergeTol",
        "scalar",
        "The merge distance relative to the bounding box size (default 1e-7)"
    );
    argList::addBoolOption
    (
        "fullMatch",
        "Do (slower) geometric matching on all boundary faces"
    );
    argList::addBoolOption
    (
        "procMatch",
        "Do matching on processor faces only"
    );
    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress finiteArea mesh reconstruction",
        true  // Advanced option
    );

    #include "addAllRegionOptions.H"

    #include "setRootCase.H"
    #include "createTime.H"

    printWarning();

    const bool fullMatch = args.found("fullMatch");
    const bool procMatch = args.found("procMatch");
    const bool writeCellDist = args.found("cellDist");
    bool doFiniteArea = !args.found("no-finite-area");
    const bool writeAddrOnly = args.found("addressing-only");

    const scalar mergeTol =
        args.getOrDefault<scalar>("mergeTol", defaultMergeTol);

    if (fullMatch)
    {
        Info<< "Use geometric matching on all boundary faces." << nl << endl;
    }
    else if (procMatch)
    {
        Info<< "Use geometric matching on correct procBoundaries only." << nl
            << "This assumes a correct decomposition." << endl;
    }
    else
    {
        Info<< "Merge assuming correct, fully matched procBoundaries." << nl
            << endl;
    }

    if (fullMatch || procMatch)
    {
        const scalar writeTol =
            std::pow(scalar(10), -scalar(IOstream::defaultPrecision()));

        Info<< "Merge tolerance : " << mergeTol << nl
            << "Write tolerance : " << writeTol << endl;

        if
        (
            runTime.writeFormat() == IOstreamOption::ASCII
         && mergeTol < writeTol
        )
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << endl
                << "Your merging tolerance (" << mergeTol << ")"
                << " is finer than this." << endl
                << "Please change your writeFormat to binary"
                << " or increase the writePrecision" << endl
                << "or adjust the merge tolerance (-mergeTol)."
                << exit(FatalError);
        }
    }

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

    // Read all time databases
    PtrList<Time> databases(nProcs);

    Info<< "Found " << nProcs << " processor directories" << endl;
    forAll(databases, proci)
    {
        Info<< "    Reading database "
            << args.caseName()/("processor" + Foam::name(proci))
            << endl;

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
    Info<< endl;

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    // Loop over all times
    forAll(timeDirs, timei)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timei], timei);

        // Set time for all databases
        forAll(databases, proci)
        {
            databases[proci].setTime(timeDirs[timei], timei);
        }

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for any mesh changes
        label nMeshChanged = 0;
        boolList hasRegionMesh(regionNames.size(), false);
        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];

            IOobject facesIO
            (
                "faces",
                databases[0].timeName(),
                polyMesh::meshDir(regionName),
                databases[0],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );

            // Problem: faceCompactIOList recognises both 'faceList' and
            //          'faceCompactList' so we should be lenient when doing
            //          typeHeaderOk

            const bool ok = facesIO.typeHeaderOk<faceCompactIOList>(false);

            if (ok)
            {
                hasRegionMesh[regioni] = true;
                ++nMeshChanged;
            }
        }

        // Check for any mesh changes
        if (!nMeshChanged)
        {
            Info<< "No meshes" << nl << endl;
            continue;
        }

        Info<< endl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir = polyMesh::regionName(regionName);

            if (!hasRegionMesh[regioni])
            {
                Info<< "region=" << regionName << " (no mesh)" << nl << endl;
                continue;
            }

            if (regionNames.size() > 1)
            {
                Info<< "region=" << regionName << nl;
            }


            // Addressing from processor to reconstructed case
            labelListList cellProcAddressing(nProcs);
            labelListList faceProcAddressing(nProcs);
            labelListList pointProcAddressing(nProcs);
            labelListList boundProcAddressing(nProcs);


            // Internal faces on the final reconstructed mesh
            label masterInternalFaces;

            // Owner addressing on the final reconstructed mesh
            labelList masterOwner;

            if (procMatch)
            {
                // Read point on individual processors to determine
                // merge tolerance
                // (otherwise single cell domains might give problems)

                const boundBox bb = procBounds(databases, regionDir);
                const scalar mergeDist = mergeTol*bb.mag();

                Info<< "Overall mesh bounding box  : " << bb << nl
                    << "Relative tolerance         : " << mergeTol << nl
                    << "Absolute matching distance : " << mergeDist << nl
                    << endl;


                // Construct empty mesh.
                PtrList<fvMesh> masterMesh(nProcs);

                for (label proci=0; proci<nProcs; proci++)
                {
                    masterMesh.set
                    (
                        proci,
                        new fvMesh
                        (
                            IOobject
                            (
                                regionName,
                                runTime.timeName(),
                                runTime,
                                IOobject::NO_READ
                            ),
                            Zero
                        )
                    );

                    fvMesh meshToAdd
                    (
                        IOobject
                        (
                            regionName,
                            databases[proci].timeName(),
                            databases[proci]
                        )
                    );

                    // Initialize its addressing
                    cellProcAddressing[proci] = identity(meshToAdd.nCells());
                    faceProcAddressing[proci] = identity(meshToAdd.nFaces());
                    pointProcAddressing[proci] = identity(meshToAdd.nPoints());
                    boundProcAddressing[proci] =
                        identity(meshToAdd.boundaryMesh().size());

                    // Find geometrically shared points/faces.
                    autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                    (
                        fullMatch,
                        proci,
                        proci,
                        masterMesh[proci],
                        proci,
                        proci,
                        meshToAdd,
                        mergeDist
                    );

                    // Add elements to mesh
                    autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                    (
                        masterMesh[proci],
                        meshToAdd,
                        couples()
                    );

                    // Added processor
                    renumber(map().addedCellMap(), cellProcAddressing[proci]);
                    renumber(map().addedFaceMap(), faceProcAddressing[proci]);
                    renumber(map().addedPointMap(), pointProcAddressing[proci]);
                    renumber(map().addedPatchMap(), boundProcAddressing[proci]);
                }
                for (label step=2; step<nProcs*2; step*=2)
                {
                    for (label proci=0; proci<nProcs; proci+=step)
                    {
                        label next = proci + step/2;
                        if(next >= nProcs)
                        {
                            continue;
                        }

                        Info<< "Merging mesh " << proci << " with "
                            << next << endl;

                        // Find geometrically shared points/faces.
                        autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                        (
                            fullMatch,
                            proci,
                            next,
                            masterMesh[proci],
                            next,
                            proci+step,
                            masterMesh[next],
                            mergeDist
                        );

                        // Add elements to mesh
                        autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                        (
                            masterMesh[proci],
                            masterMesh[next],
                            couples()
                        );

                        // Processors that were already in masterMesh
                        for (label mergedI=proci; mergedI<next; mergedI++)
                        {
                            renumber
                            (
                                map().oldCellMap(),
                                cellProcAddressing[mergedI]
                            );

                            renumber
                            (
                                map().oldFaceMap(),
                                faceProcAddressing[mergedI]
                            );

                            renumber
                            (
                                map().oldPointMap(),
                                pointProcAddressing[mergedI]
                            );

                            // Note: boundary is special since can contain -1.
                            renumber
                            (
                                map().oldPatchMap(),
                                boundProcAddressing[mergedI]
                            );
                        }

                        // Added processor
                        for
                        (
                            label addedI=next;
                            addedI<min(proci+step, nProcs);
                            addedI++
                        )
                        {
                            renumber
                            (
                                map().addedCellMap(),
                                cellProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedFaceMap(),
                                faceProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedPointMap(),
                                pointProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedPatchMap(),
                                boundProcAddressing[addedI]
                            );
                        }

                        masterMesh.set(next, nullptr);
                    }
                }

                for (label proci=0; proci<nProcs; proci++)
                {
                    Info<< "Reading mesh to add from "
                        << databases[proci].caseName()
                        << " for time = " << databases[proci].timeName()
                        << nl << nl << endl;
                }

                // See if any points on the mastermesh have become connected
                // because of connections through processor meshes.
                mergeSharedPoints(mergeDist,masterMesh[0],pointProcAddressing);

                // Save some properties on the reconstructed mesh
                masterInternalFaces = masterMesh[0].nInternalFaces();
                masterOwner = masterMesh[0].faceOwner();


                if (writeAddrOnly)
                {
                    Info<< nl
                        << "Disabled writing of merged mesh (-addressing-only)"
                        << nl << nl;
                }
                else
                {
                    // Write reconstructed mesh
                    writeMesh(masterMesh[0], cellProcAddressing);
                    if (writeCellDist)
                    {
                        writeDistribution
                        (
                            runTime,
                            masterMesh[0],
                            cellProcAddressing
                        );
                    }
                }
            }
            else
            {
                // Load all meshes
                PtrList<fvMesh> fvMeshes(nProcs);
                for (label proci=0; proci<nProcs; proci++)
                {
                    fvMeshes.set
                    (
                        proci,
                        new fvMesh
                        (
                            IOobject
                            (
                                regionName,
                                databases[proci].timeName(),
                                databases[proci]
                            )
                        )
                    );
                }

                // Construct pointers to meshes
                UPtrList<polyMesh> meshes(fvMeshes.size());
                forAll(fvMeshes, proci)
                {
                    meshes.set(proci, &fvMeshes[proci]);
                }

                // Get pairs of patches to stitch. These pairs have to
                // - have ordered, opposite faces (so one to one correspondence)
                List<DynamicList<label>> localPatch;
                List<DynamicList<label>> remoteProc;
                List<DynamicList<label>> remotePatch;
                const label nGlobalPatches = polyMeshAdder::procPatchPairs
                (
                    meshes,
                    localPatch,
                    remoteProc,
                    remotePatch
                );

                // Collect matching boundary faces on patches-to-stitch
                labelListList localBoundaryFace;
                labelListList remoteFaceProc;
                labelListList remoteBoundaryFace;
                polyMeshAdder::patchFacePairs
                (
                    meshes,
                    localPatch,
                    remoteProc,
                    remotePatch,
                    localBoundaryFace,
                    remoteFaceProc,
                    remoteBoundaryFace
                );

                // All matched faces assumed to have vertex0 matched
                labelListList remoteFaceStart(meshes.size());
                forAll(meshes, proci)
                {
                    const labelList& procFaces = localBoundaryFace[proci];
                    remoteFaceStart[proci].setSize(procFaces.size(), 0);
                }



                labelListList patchMap(meshes.size());
                labelListList pointZoneMap(meshes.size());
                labelListList faceZoneMap(meshes.size());
                labelListList cellZoneMap(meshes.size());
                forAll(meshes, proci)
                {
                    const polyMesh& mesh = meshes[proci];
                    patchMap[proci] = identity(mesh.boundaryMesh().size());

                    // Remove excess patches
                    patchMap[proci].setSize(nGlobalPatches);

                    pointZoneMap[proci] = identity(mesh.pointZones().size());
                    faceZoneMap[proci] = identity(mesh.faceZones().size());
                    cellZoneMap[proci] = identity(mesh.cellZones().size());
                }


                refPtr<fvMesh> masterMeshPtr;
                {
                    // Do in-place addition on proc0.

                    const labelList oldFaceOwner(fvMeshes[0].faceOwner());

                    fvMeshAdder::add
                    (
                        0,              // index of mesh to modify (== mesh_)
                        fvMeshes,
                        oldFaceOwner,

                        // Coupling info
                        localBoundaryFace,
                        remoteFaceProc,
                        remoteBoundaryFace,

                        boundProcAddressing,
                        cellProcAddressing,
                        faceProcAddressing,
                        pointProcAddressing
                    );

                    // Remove zero-faces processor patches
                    const polyBoundaryMesh& pbm = fvMeshes[0].boundaryMesh();
                    labelList oldToNew(pbm.size(), -1);
                    label newi = 0;
                    // Non processor patches first
                    forAll(pbm, patchi)
                    {
                        const auto& pp = pbm[patchi];
                        if (!isA<processorPolyPatch>(pp) || pp.size())
                        {
                            oldToNew[patchi] = newi++;
                        }
                    }
                    const label nNonProcPatches = newi;

                    // Move all deletable patches to the end
                    forAll(oldToNew, patchi)
                    {
                        if (oldToNew[patchi] == -1)
                        {
                            oldToNew[patchi] = newi++;
                        }
                    }
                    fvMeshTools::reorderPatches
                    (
                        fvMeshes[0],
                        oldToNew,
                        nNonProcPatches,
                        false
                    );

                    masterMeshPtr.cref(fvMeshes[0]);
                }


                const fvMesh& masterMesh = masterMeshPtr();

                // Number of internal faces on the final reconstructed mesh
                masterInternalFaces = masterMesh.nInternalFaces();

                // Owner addressing on the final reconstructed mesh
                masterOwner = masterMesh.faceOwner();


                // Write reconstructed mesh
                // Override:
                //  - caseName
                //  - processorCase flag
                // so the resulting mesh goes to the correct location (even with
                // collated). The better way of solving this is to construct
                // (zero) mesh on the undecomposed runTime.

                if (writeAddrOnly)
                {
                    Info<< nl
                        << "Disabled writing of merged mesh (-addressing-only)"
                        << nl << nl;
                }
                else
                {
                    Time& masterTime = const_cast<Time&>(masterMesh.time());

                    const word oldCaseName = masterTime.caseName();
                    masterTime.caseName() = runTime.caseName();
                    const bool oldProcCase(masterTime.processorCase(false));

                    // Write reconstructed mesh
                    writeMesh(masterMesh, cellProcAddressing);
                    if (writeCellDist)
                    {
                        writeDistribution
                        (
                            runTime,
                            masterMesh,
                            cellProcAddressing
                        );
                    }
                    masterTime.caseName() = oldCaseName;
                    masterTime.processorCase(oldProcCase);
                }
            }


            // Write the addressing

            Info<< "Reconstructing addressing from processor meshes"
                << " to the newly reconstructed mesh" << nl << endl;


            // Read finite-area
            PtrList<faMesh> procFaMeshes(databases.size());
            PtrList<polyMesh> procMeshes(databases.size());

            forAll(databases, proci)
            {
                Info<< "Processor " << proci << nl
                    << "Read processor mesh: "
                    << (databases[proci].caseName()/regionDir) << endl;

                procMeshes.set
                (
                    proci,
                    new polyMesh
                    (
                        IOobject
                        (
                            regionName,
                            databases[proci].timeName(),
                            databases[proci]
                        )
                    )
                );
                const polyMesh& procMesh = procMeshes[proci];

                writeMaps
                (
                    masterInternalFaces,
                    masterOwner,
                    procMesh,
                    cellProcAddressing[proci],
                    faceProcAddressing[proci],
                    pointProcAddressing[proci],
                    boundProcAddressing[proci]
                );


                if (doFiniteArea)
                {
                    const word boundaryInst =
                        procMesh.time().findInstance
                        (
                            procMesh.meshDir(),
                            "boundary"
                        );

                    IOobject io
                    (
                        "faBoundary",
                        boundaryInst,
                        faMesh::meshDir(procMesh, word::null),
                        procMesh.time(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    );

                    if (io.typeHeaderOk<faBoundaryMesh>(true))
                    {
                        // Always based on the volume decomposition!
                        procFaMeshes.set(proci, new faMesh(procMesh));
                    }
                }
            }


            // Finite-area mapping
            doFiniteArea = false;
            forAll(procFaMeshes, proci)
            {
                if (procFaMeshes.set(proci))
                {
                    doFiniteArea = true;
                }
            }

            if (doFiniteArea)
            {
                // Addressing from processor to reconstructed case
                labelListList faFaceProcAddressing(nProcs);
                labelListList faEdgeProcAddressing(nProcs);
                labelListList faPointProcAddressing(nProcs);
                labelListList faBoundProcAddressing(nProcs);


                // boundProcAddressing
                // ~~~~~~~~~~~~~~~~~~~

                forAll(procFaMeshes, proci)
                {
                    const auto& procMesh = procFaMeshes[proci];
                    const auto& bm = procMesh.boundary();

                    faBoundProcAddressing[proci] = identity(bm.size());
                    // Mark processor patches
                    for
                    (
                        label patchi = bm.nNonProcessor();
                        patchi < bm.size();
                        ++patchi
                    )
                    {
                        faBoundProcAddressing[proci][patchi] = -1;
                    }
                }


                // Re-read reconstructed polyMesh. Note: could probably be avoided
                // by merging into loops above.
                const polyMesh masterMesh
                (
                    IOobject
                    (
                        regionName,
                        runTime.timeName(),
                        runTime
                    ),
                    true
                );


                // faceProcAddressing
                // ~~~~~~~~~~~~~~~~~~

                DynamicList<label> masterFaceLabels;
                label nPatchFaces = 0;
                forAll(procFaMeshes, proci)
                {
                    const auto& procMesh = procFaMeshes[proci];
                    const auto& procPolyFaces = procMesh.faceLabels();
                    const auto& fpa = faceProcAddressing[proci];

                    labelList& faceAddr = faFaceProcAddressing[proci];
                    faceAddr.resize_nocopy(procPolyFaces.size());

                    // Map to masterPolyMesh faces
                    forAll(procPolyFaces, i)
                    {
                        const label facei = procPolyFaces[i];
                        masterFaceLabels.append(fpa[facei]);
                        faceAddr[i] = nPatchFaces++;
                    }
                }


                // faMesh itself
                // ~~~~~~~~~~~~~

                // Set up to read-if-present
                IOobject io(masterMesh);
                io.readOpt(IOobject::READ_IF_PRESENT);

                // Construct without patches
                faMesh masterFaMesh
                (
                    masterMesh,
                    std::move(masterFaceLabels),
                    io
                );

                const uindirectPrimitivePatch& masterPatch =
                    masterFaMesh.patch();


                // pointProcAddressing
                // ~~~~~~~~~~~~~~~~~~~

                const auto& mpm = masterPatch.meshPointMap();
                forAll(procFaMeshes, proci)
                {
                    const auto& procPatch = procFaMeshes[proci].patch();
                    const auto& mp = procPatch.meshPoints();

                    labelList& pointAddr = faPointProcAddressing[proci];
                    pointAddr.resize_nocopy(mp.size());

                    forAll(mp, i)
                    {
                        pointAddr[i] = mpm[pointProcAddressing[proci][mp[i]]];
                    }
                }


                // edgeProcAddressing
                // ~~~~~~~~~~~~~~~~~~
                // 1. straight mapping from proc faMesh back to reconstructed
                //    faMesh
                determineFaEdgeMapping
                (
                    masterPatch,            // reconstructed faMesh patch
                    procFaMeshes,           // individual faMeshes
                    pointProcAddressing,    // procPolyMesh to reconstructed

                    faEdgeProcAddressing
                );


                // Per patch all edges on the masterPatch
                labelListList singlePatchEdgeLabels;

                // 2. edgeProcAddressing : fix ordering on patches
                sortFaEdgeMapping
                (
                    masterPatch,           // reconstructed faMesh patch
                    procFaMeshes,          // individual faMeshes
                    pointProcAddressing,   // procPolyMesh to reconstructed

                    faEdgeProcAddressing,  // per proc the map to the master
                    singlePatchEdgeLabels  // per patch the map to the master
                );


                const faMesh& procMesh0 = procFaMeshes[0];


                // Add faPatches
                // ~~~~~~~~~~~~~

                faPatchList completePatches(singlePatchEdgeLabels.size());
                label nPatches = 0;
                forAll(completePatches, patchi)
                {
                    const labelList& patchEdgeLabels =
                        singlePatchEdgeLabels[patchi];

                    const faPatch& fap = procMesh0.boundary()[patchi];

                    if (isA<ignoreFaPatch>(fap))
                    {
                        // These are not real edges
                        continue;
                    }

                    const label neiPolyPatchId = fap.ngbPolyPatchIndex();

                    completePatches.set
                    (
                        nPatches,
                        fap.clone
                        (
                            masterFaMesh.boundary(),
                            patchEdgeLabels,
                            nPatches,  // index
                            neiPolyPatchId
                        )
                    );

                    ++nPatches;
                }

                completePatches.resize(nPatches);

                // Serial mesh - no parallel communication

                const bool oldParRun = Pstream::parRun(false);

                masterFaMesh.addFaPatches(completePatches);

                Pstream::parRun(oldParRun);  // Restore parallel state


                // Write mesh & individual addressing
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                faMeshReconstructor::writeMesh
                (
                    masterMesh.time().timeName(),
                    masterFaMesh,
                    masterFaMesh.faceLabels()
                );

                forAll(procFaMeshes, proci)
                {
                    const faMesh& procMesh = procFaMeshes[proci];

                    faMeshReconstructor::writeAddressing
                    (
                        IOobject
                        (
                            "procAddressing",
                            masterMesh.time().timeName(),
                            faMesh::meshSubDir,
                            procMesh.thisDb(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            IOobject::NO_REGISTER
                        ),
                        faBoundProcAddressing[proci],
                        faFaceProcAddressing[proci],
                        faPointProcAddressing[proci],
                        faEdgeProcAddressing[proci]
                    );
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
