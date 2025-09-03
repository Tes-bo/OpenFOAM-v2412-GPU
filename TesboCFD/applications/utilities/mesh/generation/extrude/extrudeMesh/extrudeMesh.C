/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    extrudeMesh

Group
    grpMeshGenerationUtilities

Description
    Extrude mesh from existing patch (by default outwards facing normals;
    optional flips faces) or from patch read from file.

    Note: Merges close points so be careful.

    Type of extrusion prescribed by run-time selectable model.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "edgeCollapser.H"
#include "perfectInterface.H"
#include "addPatchCellLayer.H"
#include "fvMesh.H"
#include "MeshedSurfaces.H"
#include "globalIndex.H"
#include "cellSet.H"
#include "fvMeshTools.H"

#include "extrudedMesh.H"
#include "extrudeModel.H"

#include "wedge.H"
#include "wedgePolyPatch.H"
#include "planeExtrusion.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "processorMeshes.H"
#include "hexRef8Data.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum ExtrudeMode
{
    MESH,
    PATCH,
    SURFACE
};

static const Enum<ExtrudeMode> ExtrudeModeNames
{
    { ExtrudeMode::MESH, "mesh" },
    { ExtrudeMode::PATCH, "patch" },
    { ExtrudeMode::SURFACE, "surface" },
};


label findPatchID(const polyBoundaryMesh& patches, const word& name)
{
    const label patchID = patches.findPatchID(name);

    if (patchID == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << name
            << " in the source mesh.\n"
            << "Valid patch names are " << patches.names()
            << exit(FatalError);
    }
    return patchID;
}


labelList patchFaces(const polyBoundaryMesh& patches, const wordRes& names)
{
    const labelList patchIDs(patches.indices(names));

    label n = 0;

    for (label patchi : patchIDs)
    {
        n += patches[patchi].size();
    }
    labelList faceLabels(n);
    n = 0;
    for (label patchi : patchIDs)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, j)
        {
            faceLabels[n++] = pp.start()+j;
        }
    }

    return faceLabels;
}


void zoneFaces
(
    const faceZoneMesh& fzs,
    const wordRes& names,
    labelList& faceLabels,
    bitSet& faceFlip
)
{
    const labelList zoneIDs(fzs.indices(names));

    label n = 0;

    for (label zonei : zoneIDs)
    {
        n += fzs[zonei].size();
    }
    faceLabels.setSize(n);
    faceFlip.setSize(n);
    n = 0;
    for (label zonei : zoneIDs)
    {
        const auto& pp = fzs[zonei];
        const boolList& ppFlip = pp.flipMap();
        forAll(pp, i)
        {
            faceLabels[n] = pp[i];
            faceFlip[n] = ppFlip[i];
            n++;
        }
    }
}


void updateFaceLabels(const mapPolyMesh& map, labelList& faceLabels)
{
    const labelList& reverseMap = map.reverseFaceMap();

    labelList newFaceLabels(faceLabels.size());
    label newI = 0;

    forAll(faceLabels, i)
    {
        label oldFacei = faceLabels[i];

        if (reverseMap[oldFacei] >= 0)
        {
            newFaceLabels[newI++] = reverseMap[oldFacei];
        }
    }
    newFaceLabels.setSize(newI);
    faceLabels.transfer(newFaceLabels);
}


void updateCellSet(const mapPolyMesh& map, labelHashSet& cellLabels)
{
    const labelList& reverseMap = map.reverseCellMap();

    labelHashSet newCellLabels(2*cellLabels.size());

    forAll(cellLabels, i)
    {
        label oldCelli = cellLabels[i];

        if (reverseMap[oldCelli] >= 0)
        {
            newCellLabels.insert(reverseMap[oldCelli]);
        }
    }
    cellLabels.transfer(newCellLabels);
}


template<class PatchType>
void changeFrontBackPatches
(
    polyMesh& mesh,
    const word& frontPatchName,
    const word& backPatchName
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label frontPatchi = findPatchID(patches, frontPatchName);
    label backPatchi = findPatchID(patches, backPatchName);

    DynamicList<polyPatch*> newPatches(patches.size());

    forAll(patches, patchi)
    {
        const polyPatch& pp(patches[patchi]);

        if (patchi == frontPatchi || patchi == backPatchi)
        {
            newPatches.append
            (
                new PatchType
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    pp.index(),
                    mesh.boundaryMesh(),
                    PatchType::typeName
                )
            );
        }
        else
        {
            newPatches.append(pp.clone(mesh.boundaryMesh()).ptr());
        }
    }

    // Edit patches
    mesh.removeBoundary();
    mesh.addPatches(newPatches, true);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extrude mesh from existing patch."
    );

    #include "addRegionOption.H"
    argList::addOption("dict", "file", "Alternative extrudeMeshDict");
    #include "setRootCase.H"
    #include "createTimeExtruded.H"

    // Specified region or default region
    #include "getRegionOption.H"

    {
        Info<< "Create mesh";
        if (!polyMesh::regionName(regionName).empty())
        {
            Info<< ' ' << regionName;
        }
        Info<< " for time = " << runTimeExtruded.timeName() << nl << nl;
    }

    const IOdictionary dict
    (
        IOobject::selectIO
        (
            IOobject
            (
                "extrudeMeshDict",
                runTimeExtruded.system(),
                runTimeExtruded,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            ),
            args.getOrDefault<fileName>("dict", "")
        )
    );

    // Point generator
    autoPtr<extrudeModel> model(extrudeModel::New(dict));

    // Whether to flip normals
    const bool flipNormals(dict.get<bool>("flipNormals"));

    // What to extrude
    const ExtrudeMode mode = ExtrudeModeNames.get("constructFrom", dict);

    // Any merging of small edges
    const scalar mergeTol(dict.getOrDefault<scalar>("mergeTol", 1e-4));

    Info<< "Extruding from " << ExtrudeModeNames[mode]
        << " using model " << model().type() << endl;
    if (flipNormals)
    {
        Info<< "Flipping normals before extruding" << endl;
    }
    if (mergeTol > 0)
    {
        Info<< "Collapsing edges < " << mergeTol << " of bounding box" << endl;
    }
    else
    {
        Info<< "Not collapsing any edges after extrusion" << endl;
    }
    Info<< endl;


    // Generated mesh (one of either)
    autoPtr<fvMesh> meshFromMesh;
    autoPtr<polyMesh> meshFromSurface;

    // Faces on front and back for stitching (in case of mergeFaces)
    word frontPatchName;
    labelList frontPatchFaces;
    word backPatchName;
    labelList backPatchFaces;

    // Optional added cells (get written to cellSet)
    labelHashSet addedCellsSet;

    // Optional refinement data
    autoPtr<hexRef8Data> refDataPtr;

    if (mode == PATCH || mode == MESH)
    {
        if (flipNormals && mode == MESH)
        {
            FatalErrorInFunction
                << "Flipping normals not supported for extrusions from mesh."
                << exit(FatalError);
        }

        fileName sourceCasePath(dict.get<fileName>("sourceCase").expand());
        fileName sourceRootDir = sourceCasePath.path();
        fileName sourceCaseDir = sourceCasePath.name();
        if (Pstream::parRun())
        {
            sourceCaseDir =
                sourceCaseDir/("processor" + Foam::name(Pstream::myProcNo()));
        }
        wordRes sourcePatches;
        wordRes sourceFaceZones;
        if
        (
            dict.readIfPresent
            (
                "sourcePatches",
                sourcePatches,
                keyType::LITERAL
            )
        )
        {
            if (sourcePatches.size() == 1)
            {
                frontPatchName = sourcePatches[0];
            }

            Info<< "Extruding patches " << sourcePatches
                << " on mesh " << sourceCasePath << nl
                << endl;
        }
        else
        {
            dict.readEntry("sourceFaceZones", sourceFaceZones);
            Info<< "Extruding faceZones " << sourceFaceZones
                << " on mesh " << sourceCasePath << nl
                << endl;
        }


        Time runTime
        (
            Time::controlDictName,
            sourceRootDir,
            sourceCaseDir
        );

        #include "createNamedMesh.H"

        const polyBoundaryMesh& patches = mesh.boundaryMesh();


        // Extrusion engine. Either adding to existing mesh or
        // creating separate mesh.
        addPatchCellLayer layerExtrude(mesh, (mode == MESH));

        labelList meshFaces;
        bitSet faceFlips;
        if (sourceFaceZones.size())
        {
            zoneFaces(mesh.faceZones(), sourceFaceZones, meshFaces, faceFlips);
        }
        else
        {
            meshFaces = patchFaces(patches, sourcePatches);
            faceFlips.setSize(meshFaces.size());
        }

        if (mode == PATCH && flipNormals)
        {
            // Cheat. Flip patch faces in mesh. This invalidates the
            // mesh (open cells) but does produce the correct extrusion.
            polyTopoChange meshMod(mesh);
            forAll(meshFaces, i)
            {
                label meshFacei = meshFaces[i];

                label patchi = patches.whichPatch(meshFacei);
                label own = mesh.faceOwner()[meshFacei];
                label nei = -1;
                if (patchi == -1)
                {
                    nei = mesh.faceNeighbour()[meshFacei];
                }

                label zoneI = mesh.faceZones().whichZone(meshFacei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    label index = mesh.faceZones()[zoneI].whichFace(meshFacei);
                    zoneFlip = mesh.faceZones()[zoneI].flipMap()[index];
                }

                meshMod.modifyFace
                (
                    mesh.faces()[meshFacei].reverseFace(),  // modified face
                    meshFacei,                      // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    true,                           // face flip
                    patchi,                         // patch for face
                    zoneI,                          // zone for face
                    zoneFlip                        // face flip in zone
                );
            }

            // Change the mesh. No inflation.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map());

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }


        indirectPrimitivePatch extrudePatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                meshFaces
            ),
            mesh.points()
        );

        // Determine extrudePatch normal
        pointField extrudePatchPointNormals
        (
            PatchTools::pointNormals(mesh, extrudePatch, faceFlips)
        );


        // Precalculate mesh edges for pp.edges.
        const labelList meshEdges
        (
            extrudePatch.meshEdges
            (
                mesh.edges(),
                mesh.pointEdges()
            )
        );

        // Global face indices engine
        const globalIndex globalFaces(mesh.nFaces());

        // Determine extrudePatch.edgeFaces in global numbering (so across
        // coupled patches)
        labelListList edgeGlobalFaces
        (
            addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                extrudePatch
            )
        );


        // Determine what patches boundary edges need to get extruded into.
        // This might actually cause edge-connected processors to become
        // face-connected so might need to introduce new processor boundaries.
        // Calculates:
        //  - per pp.edge the patch to extrude into
        //  - any additional processor boundaries (patchToNbrProc = map from
        //    new patchID to neighbour processor)
        //  - number of new patches (nPatches)

        labelList edgePatchID;
        labelList edgeZoneID;
        boolList edgeFlip;
        labelList inflateFaceID;
        label nPatches;
        Map<label> nbrProcToPatch;
        Map<label> patchToNbrProc;
        addPatchCellLayer::calcExtrudeInfo
        (
            true,   // for internal edges get zone info from any face

            mesh,
            globalFaces,
            edgeGlobalFaces,
            extrudePatch,

            edgePatchID,
            nPatches,
            nbrProcToPatch,
            patchToNbrProc,

            edgeZoneID,
            edgeFlip,
            inflateFaceID
        );


        // Add any patches.

        const label nAdded = returnReduce
        (
            nPatches - mesh.boundaryMesh().size(),
            sumOp<label>()
        );

        Info<< "Adding overall " << nAdded << " processor patches." << endl;

        if (nAdded > 0)
        {
            DynamicList<polyPatch*> newPatches(nPatches);
            forAll(mesh.boundaryMesh(), patchi)
            {
                newPatches.append
                (
                    mesh.boundaryMesh()[patchi].clone
                    (
                        mesh.boundaryMesh()
                    ).ptr()
                );
            }
            for
            (
                label patchi = mesh.boundaryMesh().size();
                patchi < nPatches;
                patchi++
            )
            {
                label nbrProci = patchToNbrProc[patchi];

                Pout<< "Adding patch " << patchi
                    << " between " << Pstream::myProcNo()
                    << " and " << nbrProci
                    << endl;

                newPatches.append
                (
                    new processorPolyPatch
                    (
                        0,                  // size
                        mesh.nFaces(),      // start
                        patchi,             // index
                        mesh.boundaryMesh(),// polyBoundaryMesh
                        Pstream::myProcNo(),// myProcNo
                        nbrProci            // neighbProcNo
                    )
                );
            }

            // Add patches. Do no parallel updates.
            mesh.removeFvBoundary();
            mesh.addFvPatches(newPatches, true);
        }



        // Only used for addPatchCellLayer into new mesh
        labelList exposedPatchID;
        if (mode == PATCH)
        {
            dict.readEntry("exposedPatchName", backPatchName);
            exposedPatchID.setSize
            (
                extrudePatch.size(),
                findPatchID(patches, backPatchName)
            );
        }

        // Determine points and extrusion
        pointField layer0Points(extrudePatch.nPoints());
        pointField displacement(extrudePatch.nPoints());
        forAll(displacement, pointi)
        {
            const vector& patchNormal = extrudePatchPointNormals[pointi];

            // layer0 point
            layer0Points[pointi] = model()
            (
                extrudePatch.localPoints()[pointi],
                patchNormal,
                0
            );
            // layerN point
            point extrudePt = model()
            (
                extrudePatch.localPoints()[pointi],
                patchNormal,
                model().nLayers()
            );
            displacement[pointi] = extrudePt - layer0Points[pointi];
        }


        // Check if wedge (has layer0 different from original patch points)
        // If so move the mesh to starting position.
        if (gMax(mag(layer0Points-extrudePatch.localPoints())) > SMALL)
        {
            Info<< "Moving mesh to layer0 points since differ from original"
                << " points - this can happen for wedge extrusions." << nl
                << endl;

            pointField newPoints(mesh.points());
            forAll(extrudePatch.meshPoints(), i)
            {
                newPoints[extrudePatch.meshPoints()[i]] = layer0Points[i];
            }
            mesh.movePoints(newPoints);
        }


        // Layers per face
        labelList nFaceLayers(extrudePatch.size(), model().nLayers());

        // Layers per point
        labelList nPointLayers(extrudePatch.nPoints(), model().nLayers());

        // Displacement for first layer
        vectorField firstLayerDisp(displacement*model().sumThickness(1));

        // Expansion ratio not used.
        scalarField ratio(extrudePatch.nPoints(), 1.0);


        // Load any refinement data
        if (mode == MESH)
        {
            // Tricky: register hexRef8 onto the database
            // since the mesh does not yet exist. It should not be registered
            // onto the input mesh.
            refDataPtr.reset
            (
                new hexRef8Data
                (
                    IOobject
                    (
                        "dummy",
                        mesh.facesInstance(),
                        polyMesh::meshSubDir,
                        runTimeExtruded,        //mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    )
                )
            );
        }


        // Topo change container. Either copy an existing mesh or start
        // with empty storage (number of patches only needed for checking)
        autoPtr<polyTopoChange> meshMod
        (
            (
                mode == MESH
              ? new polyTopoChange(mesh)
              : new polyTopoChange(patches.size())
            )
        );

        layerExtrude.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            ratio,              // expansion ratio
            extrudePatch,       // patch faces to extrude
            faceFlips,          // side to extrude (for internal/coupled faces)

            edgePatchID,        // if boundary edge: patch for extruded face
            edgeZoneID,         // optional zone for extruded face
            edgeFlip,
            inflateFaceID,      // mesh face that zone/patch info is from

            exposedPatchID,     // if new mesh: patches for exposed faces
            nFaceLayers,
            nPointLayers,
            firstLayerDisp,
            meshMod()
        );

        // Reset points according to extrusion model
        forAll(layerExtrude.addedPoints(), pointi)
        {
            const labelList& pPoints = layerExtrude.addedPoints()[pointi];
            forAll(pPoints, pPointi)
            {
                label meshPointi = pPoints[pPointi];

                point modelPt
                (
                    model()
                    (
                        extrudePatch.localPoints()[pointi],
                        extrudePatchPointNormals[pointi],
                        pPointi+1       // layer
                    )
                );

                const_cast<DynamicList<point>&>
                (
                    meshMod().points()
                )[meshPointi] = modelPt;
            }
        }

        // Store faces on front and exposed patch (if mode=patch there are
        // only added faces so cannot used map to old faces)
        const labelListList& layerFaces = layerExtrude.layerFaces();
        backPatchFaces.setSize(layerFaces.size());
        frontPatchFaces.setSize(layerFaces.size());
        forAll(backPatchFaces, patchFacei)
        {
            backPatchFaces[patchFacei]  = layerFaces[patchFacei].first();
            frontPatchFaces[patchFacei] = layerFaces[patchFacei].last();
        }


        // Create dummy fvSchemes, fvSolution
        fvMeshTools::createDummyFvMeshFiles
        (
            mesh,
            polyMesh::regionName(regionName),
            true
        );

        // Create actual mesh from polyTopoChange container
        autoPtr<mapPolyMesh> map = meshMod().makeMesh
        (
            meshFromMesh,
            IOobject
            (
                regionName,
                runTimeExtruded.constant(),
                runTimeExtruded,
                IOobject::READ_IF_PRESENT,  // Read fv* if present
                IOobject::AUTO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh
        );

        layerExtrude.updateMesh
        (
            map(),
            identity(extrudePatch.size()),
            identity(extrudePatch.nPoints())
        );

        // Calculate face labels for front and back.
        frontPatchFaces = renumber
        (
            map().reverseFaceMap(),
            frontPatchFaces
        );
        backPatchFaces = renumber
        (
            map().reverseFaceMap(),
            backPatchFaces
        );

        // Update
        if (refDataPtr)
        {
            refDataPtr->updateMesh(map());
        }

        // Store added cells
        if (mode == MESH)
        {
            const labelListList addedCells
            (
                layerExtrude.addedCells
                (
                    *meshFromMesh,
                    layerExtrude.layerFaces()
                )
            );
            forAll(addedCells, facei)
            {
                const labelList& aCells = addedCells[facei];
                addedCellsSet.insert(aCells);
            }
        }
    }
    else
    {
        // Read from surface
        fileName surfName(dict.get<fileName>("surface").expand());

        Info<< "Extruding surfaceMesh read from file " << surfName << nl
            << endl;

        MeshedSurface<face> fMesh(surfName);

        if (flipNormals)
        {
            Info<< "Flipping faces." << nl << endl;
            faceList& faces = const_cast<faceList&>(fMesh.surfFaces());
            forAll(faces, i)
            {
                faces[i] = fMesh[i].reverseFace();
            }
        }

        Info<< "Extruding surface with :" << nl
                << "    points     : " << fMesh.points().size() << nl
                << "    faces      : " << fMesh.size() << nl
                << "    normals[0] : " << fMesh.faceNormals()[0]
                << nl
                << endl;

        meshFromSurface.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                fMesh,
                model()
            )
        );


        // Get the faces on front and back
        frontPatchName = "originalPatch";
        frontPatchFaces = patchFaces
        (
            meshFromSurface().boundaryMesh(),
            wordRes(1, frontPatchName)
        );
        backPatchName = "otherSide";
        backPatchFaces = patchFaces
        (
            meshFromSurface().boundaryMesh(),
            wordRes(1, backPatchName)
        );
    }


    polyMesh& mesh =
    (
        meshFromMesh
      ? *meshFromMesh
      : *meshFromSurface
    );


    const boundBox& bb = mesh.bounds();
    const scalar mergeDim = mergeTol * bb.minDim();

    Info<< "Mesh bounding box : " << bb << nl
        << "        with span : " << bb.span() << nl
        << "Merge distance    : " << mergeDim << nl
        << endl;


    // Collapse edges
    // ~~~~~~~~~~~~~~

    if (mergeDim > 0)
    {
        Info<< "Collapsing edges < " << mergeDim << " ..." << nl << endl;

        // Edge collapsing engine
        edgeCollapser collapser(mesh);

        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        bitSet collapseEdge(mesh.nEdges());
        Map<point> collapsePointToLocation(mesh.nPoints());

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            scalar d = e.mag(points);

            if (d < mergeDim)
            {
                Info<< "Merging edge " << e << " since length " << d
                    << " << " << mergeDim << nl;

                collapseEdge.set(edgeI);
                collapsePointToLocation.set(e[1], points[e[0]]);
            }
        }

        List<pointEdgeCollapse> allPointInfo;
        const globalIndex globalPoints(mesh.nPoints());
        labelList pointPriority(mesh.nPoints(), Zero);

        collapser.consistentCollapse
        (
            globalPoints,
            pointPriority,
            collapsePointToLocation,
            collapseEdge,
            allPointInfo
        );

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Put all modifications into meshMod
        bool anyChange = collapser.setRefinement(allPointInfo, meshMod);

        if (returnReduceOr(anyChange))
        {
            // Construct new mesh from polyTopoChange.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map());

            // Update stored data
            updateFaceLabels(map(), frontPatchFaces);
            updateFaceLabels(map(), backPatchFaces);
            updateCellSet(map(), addedCellsSet);

            if (refDataPtr)
            {
                refDataPtr->updateMesh(map());
            }

            // Move mesh (if inflation used)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }
    }


    // Change the front and back patch types as required
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    word frontBackType;

    if (isType<extrudeModels::wedge>(model()))
    {
        changeFrontBackPatches<wedgePolyPatch>
        (
            mesh,
            frontPatchName,
            backPatchName
        );
    }
    else if (isType<extrudeModels::plane>(model()))
    {
        changeFrontBackPatches<emptyPolyPatch>
        (
            mesh,
            frontPatchName,
            backPatchName
        );
    }


    // Merging front and back patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (dict.get<bool>("mergeFaces"))
    {
        if (mode == MESH)
        {
            FatalErrorInFunction
                << "Cannot stitch front and back of extrusion since"
                << " in 'mesh' mode (extrusion appended to mesh)."
                << exit(FatalError);
        }

        Info<< "Assuming full 360 degree axisymmetric case;"
            << " stitching faces on patches "
            << frontPatchName << " and "
            << backPatchName << " together ..." << nl << endl;

        if (frontPatchFaces.size() != backPatchFaces.size())
        {
            FatalErrorInFunction
                << "Differing number of faces on front ("
                << frontPatchFaces.size() << ") and back ("
                << backPatchFaces.size() << ")"
                << exit(FatalError);
        }


        polyTopoChanger stitcher(mesh);
        stitcher.setSize(1);

        const word cutZoneName("originalCutFaceZone");

        List<faceZone*> fz
        (
            1,
            new faceZone
            (
                cutZoneName,
                frontPatchFaces,
                false, // none are flipped
                0,
                mesh.faceZones()
            )
        );

        mesh.addZones(List<pointZone*>(0), fz, List<cellZone*>(0));

        // Add the perfect interface mesh modifier
        perfectInterface perfectStitcher
        (
            "couple",
            0,
            stitcher,
            cutZoneName,
            word::null,         // dummy patch name
            word::null          // dummy patch name
        );

        // Topo change container
        polyTopoChange meshMod(mesh);

        perfectStitcher.setRefinement
        (
            indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh.faces(),
                    frontPatchFaces
                ),
                mesh.points()
            ),
            indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh.faces(),
                    backPatchFaces
                ),
                mesh.points()
            ),
            meshMod
        );

        // Construct new mesh from polyTopoChange.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map());

        // Update local data
        updateCellSet(map(), addedCellsSet);

        if (refDataPtr)
        {
            refDataPtr->updateMesh(map());
        }

        // Move mesh (if inflation used)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
    }

    mesh.setInstance(runTimeExtruded.constant());
    Info<< "Writing mesh to " << mesh.objectPath() << nl << endl;

    if (!mesh.write())
    {
        FatalErrorInFunction
            << exit(FatalError);
    }
    // Remove any left-over files
    processorMeshes::removeFiles(mesh);

    // Need writing cellSet
    if (returnReduceOr(addedCellsSet.size()))
    {
        cellSet addedCells(mesh, "addedCells", addedCellsSet);
        Info<< "Writing added cells to cellSet " << addedCells.name()
            << nl << endl;
        if (!addedCells.write())
        {
            FatalErrorInFunction
                << addedCells.name()
                << exit(FatalError);
        }
    }

    if (refDataPtr)
    {
        refDataPtr->write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
