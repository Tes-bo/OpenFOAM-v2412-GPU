/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    surfaceMeshExtract

Group
    grpSurfaceUtilities

Description
    Extract patch or faceZone surfaces from a polyMesh.
    Depending on output surface format triangulates faces.

    Region numbers on faces not guaranteed to be the same as the patch indices.

    Optionally only extracts named patches.

    Optionally filters out points on faceZones, feature-edges and
    featurePoints and generates pointPatches
    for these - written to pointMesh/boundary.

    If run in parallel, processor patches get filtered out by default and
    the mesh is merged (based on topology).

Usage
    \b surfaceMeshExtract [OPTION] \<surfacefile\>

    Options:
      - \par -patches NAME | LIST
        Specify single patch or multiple patches (name or regex) to extract

      - \par -faceZones NAME | LIST
        Specify single zone or multiple face zones (name or regex) to extract

      - \par -exclude-patches NAME | LIST
        Exclude single or multiple patches (name or regex) from extraction

      - \par -excludeProcPatches
        Exclude processor patches (default if parallel)

      - \par -featureAngle \<angle\>
        Extract feature edges/points and put into separate point-patches

      - \par -extractZonePoints
        Extract all face zone points and put into separate point-patches

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "ListListOps.H"
#include "stringListOps.H"  // For stringListOps::findMatching()
#include "indirectPrimitivePatch.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "timeSelector.H"
#include "meshPointPatch.H"
#include "unitConversion.H"
#include "dummyTransform.H"
#include "syncTools.H"
#include "processorPointPatch.H"
#include "pointMeshTools.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ
(
    const fileName& path,
    const pointPatch& pp
)
{
    const meshPointPatch& ppp = refCast<const meshPointPatch>(pp);
    const polyMesh& mesh = ppp.boundaryMesh().mesh().mesh();

    // Count constraints
    label maxConstraint = 0;
    const auto& constraints = ppp.constraints();
    forAll(constraints, i)
    {
        maxConstraint = max(maxConstraint, constraints[i].first());
    }
    reduce(maxConstraint, maxOp<label>());

    const pointField localPoints(mesh.points(), ppp.meshPoints());

    if (maxConstraint == 3)
    {
        OBJstream os
        (
            path
          / ppp.name()+"_fixedPoints.obj"
        );
        os.write(localPoints);
        Info<< "Written pointPatch " << ppp.name() << " to " << os.name()
            << endl;
    }
    else if (maxConstraint == 2)
    {
        OBJstream os
        (
            path
          / ppp.name()+"_slidingPoints.obj"
        );
        forAll(localPoints, i)
        {
            os.write(localPoints[i], constraints[i].second());
        }
        Info<< "Written pointPatch " << ppp.name() << " to " << os.name()
            << " as coordinates and normals"
            << endl;
    }
    else if (maxConstraint == 1)
    {
        OBJstream os
        (
            path
          / ppp.name()+"_surfacePoints.obj"
        );
        os.write(localPoints);
        Info<< "Written pointPatch " << ppp.name() << " to " << os.name()
            << " as coordinates"
            << endl;
    }
}


labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const wordRes& allow,
    const wordRes& deny
)
{
    // Name-based selection
    labelList indices
    (
        stringListOps::findMatching
        (
            patches,
            allow,
            deny,
            nameOp<polyPatch>()
        )
    );


    // Remove undesirable patches

    label count = 0;
    for (const label patchi : indices)
    {
        const polyPatch& pp = patches[patchi];

        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (Pstream::parRun() && bool(isA<processorPolyPatch>(pp)))
        {
            break; // No processor patches for parallel output
        }

        indices[count] = patchi;
        ++count;
    }

    indices.resize(count);

    return indices;
}


label addMeshPointPatches
(
    const polyMesh& mesh,
    const pointMesh& pMesh,

    const uindirectPrimitivePatch& allBoundary,
    const labelUList& faceToZone,
    const surfZoneIdentifierList& surfZones,
    const labelList& faceZoneToCompactZone,

    const scalar edgeFeatureAngle,
    const scalar pointFeatureAngle,
    const bool verbose = true,
    const bool writePoints = false
)
{
    const auto& pointBm = pMesh.boundary();
    const auto& fzs = mesh.faceZones();
    const label nPointPatches = pointBm.size();
    const pointField& points = mesh.points();


    // Feature edge(points) internal to a zone
    labelListList zoneToMeshPoints;
    List<pointConstraintList> zoneToConstraints;

    // Feature edge(points) in between zones
    labelList twoZoneMeshPoints;
    pointConstraintList twoZoneConstraints;

    // Feature points on > 2 zones
    labelList multiZoneMeshPoints;
    pointConstraintList multiZoneConstraints;

    pointMeshTools::featurePointsEdges
    (
        mesh,
        allBoundary,
        // Per boundary face to zone
        faceToZone,
        // Number of zones
        surfZones.size(),
        edgeFeatureAngle,
        //const scalar pointFeatureAngle, //not yet done

        // Feature edge(points) internal to a zone
        zoneToMeshPoints,
        zoneToConstraints,

        // Feature edge(points) in between zones
        twoZoneMeshPoints,
        twoZoneConstraints,

        // Feature points on > 2 zones
        multiZoneMeshPoints,
        multiZoneConstraints
    );


    // Add per-zone patches
    if (faceZoneToCompactZone.size())
    {
        // Calculate point normals consistent across whole patch
        const pointField pointNormals
        (
            PatchTools::pointNormals
            (
                mesh,
                allBoundary
            )
        );

        forAll(faceZoneToCompactZone, zonei)
        {
            const label compacti = faceZoneToCompactZone[zonei];
            if (compacti != -1)
            {
                const word patchName(surfZones[compacti].name());

                if (pointBm.findPatchID(patchName) == -1)
                {
                    // Extract the points originating from the faceZone. Can
                    //  - re-call featurePointsEdges with 0 feature angle so
                    //    all points go into the feature edges
                    //  - mark using faceToZone the correct points
                    //  - or assume the whole faceZone was extracted:
                    const uindirectPrimitivePatch fzPatch
                    (
                        UIndirectList<face>
                        (
                            mesh.faces(),
                            fzs[zonei].addressing()
                        ),
                        points
                    );
                    const auto& mp = fzPatch.meshPoints();

                    const vector nullVector(vector::uniform(0));

                    // Extract pointNormal (or 0) on all patch/connected points
                    vectorField meshPointNormals(mesh.nPoints(), nullVector);
                    for (const label pointi : mp)
                    {
                        const label allPointi =
                            allBoundary.meshPointMap()[pointi];
                        meshPointNormals[pointi] = pointNormals[allPointi];
                    }
                    syncTools::syncPointList
                    (
                        mesh,
                        meshPointNormals,
                        maxMagSqrEqOp<vector>(),
                        nullVector
                    );

                    // Extract indices with non-zero pointNormal
                    DynamicList<label> meshPoints(mp.size());
                    forAll(meshPointNormals, pointi)
                    {
                        if (meshPointNormals[pointi] != nullVector)
                        {
                            meshPoints.append(pointi);
                        }
                    }

                    const_cast<pointBoundaryMesh&>(pointBm).push_back
                    (
                        new meshPointPatch
                        (
                            patchName,
                            meshPoints,
                            vectorField(meshPointNormals, meshPoints),
                            pointBm.size(),
                            pointBm,
                            meshPointPatch::typeName
                        )
                    );

                    if (verbose)
                    {
                        const auto& ppp = pointBm.last();
                        Info<< "Added zone pointPatch " << ppp.name()
                            << " with "
                            << returnReduce(meshPoints.size(), sumOp<label>())
                            << " points" << endl;
                    }
                    if (writePoints)
                    {
                        writeOBJ(mesh.path(), pointBm.last());
                    }
                }
            }
        }
    }


    // Add per-patch feature-edges
    forAll(zoneToMeshPoints, zonei)
    {
        const label nPoints =
            returnReduce(zoneToMeshPoints[zonei].size(), sumOp<label>());

        const word patchName(surfZones[zonei].name() + "Edges");

        if (nPoints && (pointBm.findPatchID(patchName) == -1))
        {
            const_cast<pointBoundaryMesh&>(pointBm).push_back
            (
                new meshPointPatch
                (
                    patchName,
                    zoneToMeshPoints[zonei],
                    zoneToConstraints[zonei],
                    pointBm.size(),
                    pointBm,
                    meshPointPatch::typeName
                )
            );

            if (verbose)
            {
                const auto& ppp = pointBm.last();
                Info<< "Added feature-edges pointPatch " << ppp.name()
                    << " with " << nPoints << " points" << endl;
            }
            if (writePoints)
            {
                writeOBJ(mesh.path(), pointBm.last());
            }
        }
    }


    // Add inter-patch points

    const word allEdgePatchName("boundaryEdges");
    const label nPatchEdgePoints =
        returnReduce(twoZoneMeshPoints.size(), sumOp<label>());
    if (nPatchEdgePoints && (pointBm.findPatchID(allEdgePatchName) == -1))
    {
        const_cast<pointBoundaryMesh&>(pointBm).push_back
        (
            new meshPointPatch
            (
                allEdgePatchName,
                twoZoneMeshPoints,
                twoZoneConstraints,
                pointBm.size(),
                pointBm,
                meshPointPatch::typeName
            )
        );

        if (verbose)
        {
            const auto& ppp = pointBm.last();
            Info<< "Added inter-patch pointPatch " << ppp.name()
                << " with " << nPatchEdgePoints << " points" << endl;
        }
        if (writePoints)
        {
            writeOBJ(mesh.path(), pointBm.last());
        }
    }


    const word allPointPatchName("boundaryPoints");
    const label nMultiPoints =
        returnReduce(multiZoneMeshPoints.size(), sumOp<label>());
    if (nMultiPoints && (pointBm.findPatchID(allPointPatchName) == -1))
    {
        const_cast<pointBoundaryMesh&>(pointBm).push_back
        (
            new meshPointPatch
            (
                allPointPatchName,
                multiZoneMeshPoints,
                multiZoneConstraints,
                pointBm.size(),
                pointBm,
                meshPointPatch::typeName
            )
        );

        if (verbose)
        {
            const auto& ppp = pointBm.last();
            Info<< "Added multi-patch pointPatch " << ppp.name()
                << " with " << nMultiPoints << " points" << endl;
        }
        if (writePoints)
        {
            writeOBJ(mesh.path(), pointBm.last());
        }
    }


    // Shuffle into order
    labelList oldToNew(pointBm.size());
    label newPatchi = 0;
    forAll(pointBm, patchi)
    {
        if (!isA<processorPointPatch>(pointBm[patchi]))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }
    forAll(pointBm, patchi)
    {
        if (isA<processorPointPatch>(pointBm[patchi]))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }
    const_cast<pointBoundaryMesh&>(pointBm).reorder(oldToNew, true);

    return pointBm.size() - nPointPatches;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract patch or faceZone surfaces from a polyMesh."
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addArgument("output", "The output surface file");

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "excludeProcPatches",
        "Exclude processor patches"
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to extract\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'"
    );
    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to extract.\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Specify single patch or multiple patches to exclude from -patches."
        " Eg, 'outlet' or '( inlet \".*Wall\" )'",
        true  // mark as an advanced option
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2306});
    argList::addOption
    (
        "featureAngle",
        "angle",
        "Auto-extract feature edges/points and put into separate point-patches"
    );
    argList::addBoolOption
    (
        "extractZonePoints",
        "Extract point-patches for selected faceZones"
    );
    argList::addBoolOption
    (
        "writeOBJ",
        "Write added pointPatch points to .obj files"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const auto userOutFileName = args.get<fileName>(1);

    if (!userOutFileName.has_ext())
    {
        FatalErrorInFunction
            << "Missing extension on output name " << userOutFileName
            << exit(FatalError);
    }

    Info<< "Extracting surface from boundaryMesh ..." << nl << nl;

    const bool includeProcPatches =
        (!UPstream::parRun() && !args.found("excludeProcPatches"));

    if (includeProcPatches)
    {
        Info<< "Including all processor patches." << nl << endl;
    }
    else if (UPstream::parRun())
    {
        Info<< "Excluding all processor patches." << nl << endl;
    }

    wordRes includePatches, excludePatches;
    if (args.readListIfPresent<wordRe>("patches", includePatches))
    {
        Info<< "Including patches " << flatOutput(includePatches)
            << nl << endl;
    }
    if (args.readListIfPresent<wordRe>("exclude-patches", excludePatches))
    {
        Info<< "Excluding patches " << flatOutput(excludePatches)
            << nl << endl;
    }

    // Non-mandatory
    const wordRes selectedFaceZones(args.getList<wordRe>("faceZones", false));
    if (selectedFaceZones.size())
    {
        Info<< "Including faceZones " << flatOutput(selectedFaceZones)
            << nl << endl;
    }

    scalar featureAngle = 180.0;
    const bool specifiedFeature = args.readIfPresent
    (
        "featureAngle",
        featureAngle
    );

    const bool extractZonePoints = args.found("extractZonePoints");
    const bool writeOBJ = args.found("writeOBJ");

    Info<< "Reading mesh from time " << runTime.value() << endl;

    #include "createNamedPolyMesh.H"
    if (specifiedFeature)
    {
        Info<< "Detecting all sharp (>" << featureAngle
            << " degrees) patch edges." << nl << endl;

        if (extractZonePoints)
        {
            Info<< "Extracting all faceZone points as pointPatches."
                << nl << endl;
        }

        //#include "createNamedPointMesh.H"
        // Do not read constant/pointMesh - construct from polyMesh only
        Info<< "Create pointMesh for time = "
             << runTime.timeName() << Foam::nl << Foam::endl;
        (void)pointMesh::New(mesh);
    }


    // User specified times
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        Info<< "Time [" << timeIndex << "] = " << runTime.timeName();

        fileName outFileName;
        if (timeDirs.size() == 1)
        {
            outFileName = userOutFileName;
        }
        else
        {
            polyMesh::readUpdateState meshState = mesh.readUpdate();
            if (timeIndex && meshState == polyMesh::UNCHANGED)
            {
                Info<<"  ... no mesh change." << nl;
                continue;
            }

            // The filename based on the original, but with additional
            // time information. The extension was previously checked that
            // it exists
            const auto dot = userOutFileName.rfind('.');

            outFileName =
                userOutFileName.substr(0, dot) + "_"
              + Foam::name(runTime.value()) + "."
              + userOutFileName.ext();
        }

        Info<< nl;

        // Create local surface from:
        // - explicitly named patches only (-patches (at your option)
        // - all patches (default in sequential mode)
        // - all non-processor patches (default in parallel mode)
        // - all non-processor patches (sequential mode, -excludeProcPatches
        //   (at your option)

        // Construct table of patches to include.
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        const labelList patchIds =
        (
            (includePatches.size() || excludePatches.size())
          ? getSelectedPatches(bMesh, includePatches, excludePatches)
          : includeProcPatches
          ? identity(bMesh.size())
          : identity(bMesh.nNonProcessor())
        );

        labelList faceZoneIds;

        const faceZoneMesh& fzm = mesh.faceZones();

        if (selectedFaceZones.size())
        {
            faceZoneIds = fzm.indices(selectedFaceZones);

            Info<< "Additionally extracting faceZones "
                << fzm.names(selectedFaceZones) << nl;
        }


        // From (name of) patch to compact 'zone' index
        HashTable<label> compactZoneID(1024);
        // Mesh face and compact zone indx
        DynamicList<label> faceLabels;
        DynamicList<label> compactZones;
        // Per compact 'zone' index the name and location
        surfZoneIdentifierList surfZones;

        // Per local patch to compact 'zone' index (or -1)
        labelList patchToCompactZone(bMesh.size(), -1);
        // Per local faceZone to compact 'zone' index (or -1)
        labelList faceZoneToCompactZone(bMesh.size(), -1);
        {
            // Collect sizes. Hash on names to handle local-only patches (e.g.
            //  processor patches)
            HashTable<label> patchSize(1024);
            label nFaces = 0;
            for (const label patchi : patchIds)
            {
                const polyPatch& pp = bMesh[patchi];
                patchSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }

            HashTable<label> zoneSize(1024);
            for (const label zonei : faceZoneIds)
            {
                const faceZone& pp = fzm[zonei];
                zoneSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }


            Pstream::mapCombineGather(patchSize, plusEqOp<label>());
            Pstream::mapCombineGather(zoneSize, plusEqOp<label>());

            if (Pstream::master())
            {
                // Allocate compact numbering for all patches/faceZones
                forAllConstIters(patchSize, iter)
                {
                    compactZoneID.insert(iter.key(), compactZoneID.size());
                }

                forAllConstIters(zoneSize, iter)
                {
                    compactZoneID.insert(iter.key(), compactZoneID.size());
                }
            }
            Pstream::broadcast(compactZoneID);


            // Zones
            surfZones.resize_nocopy(compactZoneID.size());
            forAllConstIters(compactZoneID, iter)
            {
                surfZones[*iter] = surfZoneIdentifier(iter.key(), *iter);
                Info<< "surfZone " << *iter
                    <<  " : "      << surfZones[*iter].name()
                    << endl;
            }


            // Rework HashTable into labelList just for speed of conversion
            forAllConstIters(compactZoneID, iter)
            {
                label patchi = bMesh.findPatchID(iter.key());
                if (patchi != -1)
                {
                    patchToCompactZone[patchi] = iter.val();
                }
                else
                {
                    label zoneI = fzm.findZoneID(iter.key());
                    faceZoneToCompactZone[zoneI] = iter.val();
                }
            }


            faceLabels.setCapacity(nFaces);
            compactZones.setCapacity(nFaces);

            // Collect faces on patches
            for (const label patchi : patchIds)
            {
                const polyPatch& pp = bMesh[patchi];
                forAll(pp, i)
                {
                    faceLabels.append(pp.start()+i);
                    compactZones.append(patchToCompactZone[pp.index()]);
                }
            }
            // Collect faces on faceZones
            for (const label zonei : faceZoneIds)
            {
                const faceZone& pp = fzm[zonei];
                forAll(pp, i)
                {
                    faceLabels.append(pp[i]);
                    compactZones.append(faceZoneToCompactZone[pp.index()]);
                }
            }
        }


        // Addressing engine for all faces
        const uindirectPrimitivePatch allBoundary
        (
            UIndirectList<face>(mesh.faces(), faceLabels),
            mesh.points()
        );


        // Find correspondence to master points
        labelList pointToGlobal;
        labelList uniqueMeshPoints;
        autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
        (
            allBoundary.meshPoints(),
            allBoundary.meshPointMap(),
            pointToGlobal,
            uniqueMeshPoints
        );

        // Gather all unique points on master
        List<pointField> gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = pointField
        (
            mesh.points(),
            uniqueMeshPoints
        );
        Pstream::gatherList(gatheredPoints);

        // Gather all faces
        List<faceList> gatheredFaces(Pstream::nProcs());
        gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
        for (face& f : gatheredFaces[Pstream::myProcNo()])
        {
            inplaceRenumber(pointToGlobal, f);
        }
        Pstream::gatherList(gatheredFaces);

        // Gather all ZoneIDs
        List<labelList> gatheredZones(Pstream::nProcs());
        gatheredZones[Pstream::myProcNo()] = compactZones;
        Pstream::gatherList(gatheredZones);

        // On master combine all points, faces, zones
        if (Pstream::master())
        {
            pointField allPoints = ListListOps::combine<pointField>
            (
                gatheredPoints,
                accessOp<pointField>()
            );
            gatheredPoints.clear();

            faceList allFaces = ListListOps::combine<faceList>
            (
                gatheredFaces,
                accessOp<faceList>()
            );
            gatheredFaces.clear();

            labelList allZones = ListListOps::combine<labelList>
            (
                gatheredZones,
                accessOp<labelList>()
            );
            gatheredZones.clear();


            UnsortedMeshedSurface<face> unsortedFace
            (
                std::move(allPoints),
                std::move(allFaces),
                std::move(allZones),
                surfZones
            );


            MeshedSurface<face> sortedFace(unsortedFace);

            fileName globalCasePath
            (
                outFileName.isAbsolute()
              ? outFileName
              : (
                    runTime.processorCase()
                  ? runTime.globalPath()/outFileName
                  : runTime.path()/outFileName
                )
            );

            Info<< "Writing merged surface to " << globalCasePath << endl;

            sortedFace.write(globalCasePath);
        }


        if (specifiedFeature)
        {
            // Add edge patches
            const auto& pMesh = pointMesh::New(mesh);

            const label nAdded = addMeshPointPatches
            (
                mesh,
                pMesh,

                allBoundary,    // all patches together
                compactZones,   // originating compactZone
                surfZones,      // per compactZone the index
                (
                    extractZonePoints
                  ? faceZoneToCompactZone   // per faceZone the compactZone
                  : labelList::null()
                ),

                featureAngle,
                featureAngle,

                true,
                writeOBJ
            );

            if (nAdded)
            {
                pMesh.boundary().write();
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
