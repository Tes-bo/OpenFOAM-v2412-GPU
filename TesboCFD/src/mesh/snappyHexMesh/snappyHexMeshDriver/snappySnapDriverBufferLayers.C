/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

Description
    All to do with snapping to the surface with buffer layer

\*----------------------------------------------------------------------------*/

#include "snappySnapDriver.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OBJstream.H"
#include "mapPolyMesh.H"
#include "snapParameters.H"
#include "unitConversion.H"
#include "PatchTools.H"
#include "profiling.H"
#include "addPatchCellLayer.H"
#include "snappyLayerDriver.H"
#include "weightedPosition.H"

#include "localPointRegion.H"
#include "pointConstraints.H"
#include "displacementMotionSolver.H"
#include "meshPointPatch.H"
#include "processorPointPatch.H"
#include "dummyTransform.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "tetDecomposer.H"
#include "tetMatcher.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avg
(
    const polyMesh& mesh,
    const bitSet& isMasterPoint,
    const indirectPrimitivePatch& pp,
    const pointField& localPoints
)
{
    const labelListList& pointEdges = pp.pointEdges();
    const labelList& meshPoints = pp.meshPoints();
    const edgeList& edges = pp.edges();

    Field<weightedPosition> wps
    (
        pointEdges.size(),
        pTraits<weightedPosition>::zero
    );

    // Calculate sum of all contributions (so not positions)
    forAll(pointEdges, verti)
    {
        weightedPosition& wp = wps[verti];
        for (const label edgei : pointEdges[verti])
        {
            const label otherVerti = edges[edgei].otherVertex(verti);

            if (isMasterPoint[meshPoints[otherVerti]])
            {
                wp.first() += 1.0;
                wp.second() += localPoints[otherVerti];
            }
        }
    }

    weightedPosition::syncPoints(mesh, meshPoints, wps);

    tmp<pointField> tavg(new pointField(wps.size()));
    pointField& avg = tavg.ref();

    forAll(wps, verti)
    {
        const weightedPosition& wp = wps[verti];

        if (mag(wp.first()) < VSMALL)
        {
            // Set to zero?
            avg[verti] = Zero;
        }
        else
        {
            avg[verti] = wp.second()/wp.first();
        }
    }

    return tavg;
}


Foam::tmp<Foam::pointField>
Foam::snappySnapDriver::smoothLambdaMuPatchDisplacement
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const List<labelPair>& baffles
)
{
    const bitSet isMasterPoint(syncTools::getMasterPoints(mesh));

    pointField newLocalPoints(pp.localPoints());

    const label iters = 90;
    const scalar lambda = 0.33;
    const scalar mu = 0.34;

    for (label iter = 0; iter < iters; iter++)
    {
        // Lambda
        newLocalPoints =
            (1 - lambda)*newLocalPoints
          + lambda*avg(mesh, isMasterPoint, pp, newLocalPoints);

        // Mu
        newLocalPoints =
            (1 + mu)*newLocalPoints
          - mu*avg(mesh, isMasterPoint, pp, newLocalPoints);
    }
    return newLocalPoints-pp.localPoints();
}


Foam::tmp<Foam::scalarField>
Foam::snappySnapDriver::wantedThickness
(
    const indirectPrimitivePatch& pp,
    const scalar cellSizeFraction
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& owner = mesh.faceOwner();

    // Undistorted edge length
    const scalar edge0Len =
        meshRefiner_.meshCutter().level0EdgeLength();


    tmp<scalarField> tthickness(tmp<scalarField>::New(pp.nPoints()));
    scalarField& thickness = tthickness.ref();

    labelList maxPointLevel(pp.nPoints(), labelMin);

    forAll(pp, i)
    {
        label ownLevel = cellLevel[owner[pp.addressing()[i]]];

        const face& f = pp.localFaces()[i];

        forAll(f, fp)
        {
            maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxPointLevel,
        maxEqOp<label>(),
        labelMin            // null value
    );

    forAll(thickness, pointi)
    {
        const scalar edgeLen = edge0Len/(1<<maxPointLevel[pointi]);
        thickness[pointi] = cellSizeFraction*edgeLen;
    }
    return tthickness;
}


//const Foam::pointMesh& Foam::snappySnapDriver::makePointMesh
//(
//    const indirectPrimitivePatch& pp,
//    const pointConstraintList& pointConstraints,
//    const word& allEdgePatchName,
//    const word& allPointPatchName
//) const
//{
//    fvMesh& mesh = meshRefiner_.mesh();
//
//    if (pointConstraints.size() != pp.nPoints())
//    {
//        FatalErrorInFunction<< "pointConstraints:" << pointConstraints.size()
//            << " pp:" << pp.nPoints() << exit(FatalError);
//    }
//
//
//    // Expand pointConstraints to all meshPoints
//    pointConstraintList meshPointConstraints(mesh.nPoints());
//    UIndirectList<pointConstraint>(meshPointConstraints, pp.meshPoints()) =
//        pointConstraints;
//
//    const auto combineConstraints = [&]
//    (
//        pointConstraint& x,
//        const pointConstraint& y
//    )
//    {
//        x.combine(y);
//    };
//
//
//    syncTools::syncPointList
//    (
//        mesh,
//        meshPointConstraints,
//        combineConstraints,
//        pointConstraint(),
//        dummyTransform()
//    );
//
//
//    // Sort point constraints:
//    // - 0 : does not happen?
//    // - 1 : attract to surface
//    // - 2 : attract to feature edge
//    // - 3 : attract to feature point
//
//    DynamicList<label> featEdgeMeshPoints(pointConstraints.size());
//    DynamicList<label> featPointMeshPoints(pointConstraints.size());
//    forAll(meshPointConstraints, pointi)
//    {
//        if (meshPointConstraints[pointi].first() == 2)
//        {
//            featEdgeMeshPoints.append(pointi);
//        }
//        else if (meshPointConstraints[pointi].first() == 3)
//        {
//            featPointMeshPoints.append(pointi);
//        }
//    }
//
//    // Lookup / construct (from polyPatches) the pointMesh
//    auto& pMesh = pointMesh::New(meshRefiner_.mesh());
//
//    pointBoundaryMesh& pointBm =
//        const_cast<pointBoundaryMesh&>(pMesh.boundary());
//
//    // Check if already has constraint patches
//    label edgePatchi = pointBm.findPatchID(allEdgePatchName);
//    if (edgePatchi != -1)
//    {
//        // Delete patch. TBD: clear patchGroup
//        pointBm.set(edgePatchi, nullptr);
//    }
//    label pointPatchi = pointBm.findPatchID(allPointPatchName);
//    if (pointPatchi != -1)
//    {
//        // Delete patch. TBD: clear patchGroup
//        pointBm.set(pointPatchi, nullptr);
//    }
//
//    // Add additional point patches in order:
//    // - polyPatch based
//    // - featEdge-based constraints
//    // - featPoint-based constraints (so can override edge-based constraints)
//    // - processor boundaries (or should this be all constraint patches, e.g.
//    //   symmetry plane. Note: hopefully this is already handled in the
//    //   feature extraction ...)
//
//    if (returnReduce(featEdgeMeshPoints.size(), sumOp<label>()))
//    {
//        if (edgePatchi != -1)
//        {
//            // Override patch
//            const_cast<pointBoundaryMesh&>(pointBm).set
//            (
//                edgePatchi,
//                new meshPointPatch
//                (
//                    allEdgePatchName,
//                    featEdgeMeshPoints,
//                    List<pointConstraint>
//                    (
//                        meshPointConstraints,
//                        featEdgeMeshPoints
//                    ),
//                    edgePatchi,
//                    pointBm,
//                    meshPointPatch::typeName
//                )
//            );
//        }
//        else
//        {
//            // Append
//            const_cast<pointBoundaryMesh&>(pointBm).push_back
//            (
//                new meshPointPatch
//                (
//                    allEdgePatchName,
//                    featEdgeMeshPoints,
//                    List<pointConstraint>
//                    (
//                        meshPointConstraints,
//                        featEdgeMeshPoints
//                    ),
//                    pointBm.size(),
//                    pointBm,
//                    meshPointPatch::typeName
//                )
//            );
//        }
//    }
//    if (returnReduce(featPointMeshPoints.size(), sumOp<label>()))
//    {
//        if (pointPatchi != -1)
//        {
//            // Override patch
//            const_cast<pointBoundaryMesh&>(pointBm).set
//            (
//                pointPatchi,
//                new meshPointPatch
//                (
//                    allPointPatchName,
//                    featPointMeshPoints,
//                    List<pointConstraint>
//                    (
//                        meshPointConstraints,
//                        featPointMeshPoints
//                    ),
//                    pointPatchi,
//                    pointBm,
//                    meshPointPatch::typeName
//                )
//            );
//        }
//        else
//        {
//            // Append
//            const_cast<pointBoundaryMesh&>(pointBm).push_back
//            (
//                new meshPointPatch
//                (
//                    allPointPatchName,
//                    featPointMeshPoints,
//                    List<pointConstraint>
//                    (
//                        meshPointConstraints,
//                        featPointMeshPoints
//                    ),
//                    pointBm.size(),
//                    pointBm,
//                    meshPointPatch::typeName
//                )
//            );
//        }
//    }
//
//    // Shuffle into order
//    labelList oldToNew(pointBm.size());
//    label newPatchi = 0;
//    forAll(pointBm, patchi)
//    {
//        if (!isA<processorPointPatch>(pointBm[patchi]))
//        {
//            oldToNew[patchi] = newPatchi++;
//        }
//    }
//    forAll(pointBm, patchi)
//    {
//        if (isA<processorPointPatch>(pointBm[patchi]))
//        {
//            oldToNew[patchi] = newPatchi++;
//        }
//    }
//    pointBm.reorder(oldToNew, true);
//
//    return pMesh;
//}


Foam::autoPtr<Foam::displacementMotionSolver>
Foam::snappySnapDriver::makeMotionSolver
(
    const pointMesh& pMesh,
    const dictionary& snapDict,
    const labelList& adaptPatchIDs
//    const pointConstraintList& pointConstraints,
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    tmp<pointVectorField> tallDisp
    (
        meshRefinement::makeDisplacementField
        (
            pMesh,
            adaptPatchIDs
        )
    );

    // Make sure the pointDisplacement is not registered (since
    // displacementMotionSolver itself holds it)
    tallDisp.ref().checkOut();

    autoPtr<displacementMotionSolver> motionPtr
    (
        displacementMotionSolver::New
        (
            snapDict.get<word>("solver"),
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "motionSolverDict",
                    pMesh.thisDb().time().constant(),
                    pMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                snapDict
            ),
            tallDisp(),
            pointIOField
            (
                IOobject
                (
                    "points0",
                    pMesh.thisDb().time().constant(),
                    polyMesh::meshSubDir,
                    pMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh.points()
            )
        )
    );
    return motionPtr;
}


void Foam::snappySnapDriver::setDisplacement
(
    const indirectPrimitivePatch& pp,
    const pointField& patchDisp,    // displacement w.r.t. current mesh
    const labelList& adaptPatchIDs,
    const pointField& points0,
    pointVectorField& fld           // displacement w.r.t. points0
)
{
    const pointMesh& pMesh = fld.mesh();
    const pointField& points = pMesh().points();
    const labelList& meshPoints = pp.meshPoints();

    if
    (
        (points0.size() != points.size())
     || (points0.size() != fld.size())
     || (points0.size() != pMesh.size())
     || (meshPoints.size() != patchDisp.size())
    )
    {
        FatalErrorInFunction
            << "Sizing :"
            << " points0.size():" << points0.size()
            << " points.size():" << points.size()
            << " fld.size():" << fld.size()
            << " patchDisp.size():" << patchDisp.size()
            << " meshPoints.size():" << meshPoints.size()
            << " mesh.nPoints():" << pMesh.size()
            << exit(FatalError);
    }


    // Problem is that the patchDisp might not be consistent (parallel etc)
    // across shared points so expand to mesh points

    pointField meshDisp(pMesh.size(), Zero);

    forAll(meshPoints, patchPointi)
    {
        const label meshPointi = meshPoints[patchPointi];

        meshDisp[meshPointi] =
            patchDisp[patchPointi]+points[meshPointi]-points0[meshPointi];
    }

    // Assign to bc
    pointVectorField::Boundary& bfld = fld.boundaryFieldRef();
    for (const label patchi : adaptPatchIDs)
    {
        bfld[patchi] == pointField(meshDisp, bfld[patchi].patch().meshPoints());
    }

    // Apply multi-patch constraints. Problem: most patches originate from
    // meshing so are type 'wall'. The pointConstraints are only on any points
    // remaining from the starting mesh.
    const pointConstraints& pcs = pointConstraints::New(pMesh);
    pcs.constrainDisplacement(fld, true);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::snappySnapDriver::addBufferLayers
(
    const indirectPrimitivePatch& pp,
    const pointField& thickness,
    // Layer mesh modifier
    addPatchCellLayer& addLayer
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Introduce single layer of cells. Straight from snappyLayerDriver

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges
    // between coupled
    // faces (all edges between same (global)face indices get extruded).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
        )
    );



    // Use global edge - face connectivity to
    // - disable any non-manifold extrusions
    // - boundary edges can be extruded - is handled by addPatchCellLayer
    // - what about setting numLayers to 0 to disable layers? Should that
    //   affect buffer layer addition? Guess not since buffer layer is
    //   to aid snapping ...
    labelList nFaceLayers(pp.size(), 1);
    labelList nPointLayers(pp.nPoints(), 1);
    forAll(edgeGlobalFaces, edgei)
    {
        if (edgeGlobalFaces[edgei].size() > 2)
        {
            const edge& e = pp.edges()[edgei];
            const edge meshE(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);

            nPointLayers[e[0]] = 0;
            nPointLayers[e[1]] = 0;
        }
        //else if (edgeGlobalFaces[edgei].size() == 1)
        //{
        //    const edge& e = pp.edges()[edgei];
        //    const edge meshE(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);
        //
        //    nPointLayers[e[0]] = 0;
        //    nPointLayers[e[1]] = 0;
        //}
        else if (edgeGlobalFaces[edgei].size() == 0)
        {
            const edge& e = pp.edges()[edgei];
            const edge meshE(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);
            FatalErrorInFunction << "Edge:" << meshE
                << " At:" << meshE.line(mesh.points())
                << " has no faces!" << exit(FatalError);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nPointLayers,
        minEqOp<label>(),
        labelMax            // null value
    );

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        const UIndirectList<label> pLayers(nPointLayers, f);
        if (!pLayers.found(label(1)))
        {
            nFaceLayers[facei] = 0;
        }
    }



    // Determine patches for extruded boundary edges. Adds any
    // additional processor patches (since extruding coupled edge can cause
    // additional connectivity)

    labelList edgePatchID;
    labelList edgeZoneID;
    boolList edgeFlip;
    labelList inflateFaceID;
    snappyLayerDriver::determineSidePatches
    (
        meshRefiner_,
        globalFaces,
        edgeGlobalFaces,
        pp,

        edgePatchID,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );


    // Mesh topo change engine. Insert current mesh.
    polyTopoChange meshMod(mesh);


    // Add topo regardless of whether extrudeStatus is extruderemove.
    // Not add layer if patchDisp is zero.

    addLayer.setRefinement
    (
        globalFaces,
        edgeGlobalFaces,

        scalarField(pp.nPoints(), 1),  // expansion ratio
        pp,
        bitSet(pp.size()),          // no flip

        edgePatchID,    // boundary patch for extruded boundary edges
        edgeZoneID,     // zone for extruded edges
        edgeFlip,
        inflateFaceID,


        labelList(0),   // exposed patchIDs, not used for adding layers
        nFaceLayers,
        nPointLayers,
        thickness,  //patchAttraction,

        meshMod
    );

    // Apply the stored topo changes to the current mesh.
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
    mapPolyMesh& map = *mapPtr;

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map.hasMotionPoints())
    {
        mesh.movePoints(map.preMotionPoints());
    }
    else
    {
        // Hack to remove meshPhi - mapped incorrectly. TBD.
        mesh.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh.setInstance(meshRefiner_.timeName());

    // Update numbering on layer
    addLayer.updateMesh
    (
        map,
        identity(pp.size()),
        identity(pp.nPoints())
    );

    // Update intersections
    const labelListList addedCells(addLayer.addedCells());
    bitSet isChangedFace(mesh.nFaces());
    for (const labelList& faceToCells : addedCells)
    {
        for (const label celli : faceToCells)
        {
            isChangedFace.set(mesh.cells()[celli]);
        }
    }
    meshRefiner_.updateMesh(map, isChangedFace.toc());

    if (debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing mesh-with-layer to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            meshRefiner_.timeName()
        );
    }
    return mapPtr;
}


void Foam::snappySnapDriver::doSnapBufferLayers
(
    const dictionary& snapDict,
    const dictionary& motionDict,
    const meshRefinement::FaceMergeType mergeType,
    const scalar featureCos,
    const scalar planarAngle,
    const snapParameters& snapParams
)
{
    addProfiling(snap, "snappyHexMesh::snap");
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;


    // Name of pointPatch for any points attracted to feature edges
    //const word allEdgePatchName("boundaryEdges");
    // Name of pointPatch for any points attracted to feature points
    //const word allPointPatchName("boundaryPoints");


    // faceZone handling
    // ~~~~~~~~~~~~~~~~~
    //
    // We convert all faceZones into baffles during snapping so we can use
    // a standard mesh motion (except for the mesh checking which for baffles
    // created from internal faces should check across the baffles). The state
    // is stored in two variables:
    //      baffles : pairs of boundary faces
    //      duplicateFace : from mesh face to its baffle colleague (or -1 for
    //                      normal faces)
    // There are three types of faceZones according to the faceType property:
    //
    // internal
    // --------
    // - baffles: need to be checked across
    // - duplicateFace: from face to duplicate face. Contains
    //   all faces on faceZone to prevents merging patch faces.
    //
    // baffle
    // ------
    // - baffles: no need to be checked across
    // - duplicateFace: contains all faces on faceZone to prevent
    //   merging patch faces.
    //
    // boundary
    // --------
    // - baffles: no need to be checked across. Also points get duplicated
    //            so will no longer be baffles
    // - duplicateFace: contains no faces on faceZone since both sides can
    //   merge faces independently.



    // faceZones of type internal
    const labelList internalFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::INTERNAL
            )
        )
    );


    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    {
        List<labelPair> baffles;
        labelList originatingFaceZone;
        meshRefiner_.createZoneBaffles
        (
            identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone
        );
    }


    // Duplicate points on faceZones of type boundary
    meshRefiner_.dupNonManifoldBoundaryPoints();

    // Extract baffles across internal faceZones (for checking mesh quality
    // across)
    labelPairList internalBaffles
    (
        meshRefiner_.subsetBaffles
        (
            mesh,
            internalFaceZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        )
    );


    bool doFeatures = true; //false;
    label nFeatIter = 1;
    if (snapParams.nFeatureSnap() > 0)
    {
        doFeatures = true;

        if (!dryRun_)
        {
            nFeatIter = snapParams.nFeatureSnap();
        }

        Info<< "Snapping to features in " << nFeatIter
            << " iterations ..." << endl;
    }

    // Get the labels of added patches.
    const labelList adaptPatchIDs(meshRefiner_.meshedPatches());

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
        )
    );


    if (debug)
    {
        const fileName dir(mesh.time().path());
        mkDir(dir);
        OBJstream str
        (
            dir
          / "pp_initial_" + meshRefiner_.timeName() + ".obj"
        );
        str.write
        (
            ppPtr().localFaces(),
            ppPtr().localPoints(),
            false
        );
    }


    // Maximum distance to attract to nearest feature on surface
    scalarField snapDist(calcSnapDistance(mesh, snapParams, ppPtr()));


    // Smooth patch points. Equivalent of preSmoothPatch but using mesh
    // motion solver.
    vectorField patchDisp
    (
        smoothLambdaMuPatchDisplacement
        (
            mesh,
            ppPtr(),
            internalBaffles
        )
    );
    pointField ppLocalPoints(ppPtr().localPoints()+patchDisp);

    const bool smoothInternal = true;
    if (smoothInternal)
    {
        // Create pointMesh with the correct patches
        //const pointMesh& pMesh = makePointMesh
        //(
        //    ppPtr(),
        //    patchConstraints,
        //    allEdgePatchName,
        //    allPointPatchName
        //);
        const pointMesh& pMesh = pointMesh::New(mesh);

        autoPtr<displacementMotionSolver> motionPtr
        (
            makeMotionSolver
            (
                pMesh,
                snapDict,
                adaptPatchIDs
                //patchConstraints
            )
        );

        // Insert as bc to motionSmoother. Note that this takes displacement
        // relative to points0
        setDisplacement
        (
            ppPtr(),
            patchDisp,
            adaptPatchIDs,
            motionPtr().points0(),
            motionPtr().pointDisplacement()
        );
        // Solve internal displacement
        tmp<pointField> tnewPoints(motionPtr->newPoints());

        // Move points
        mesh.movePoints(tnewPoints);

        // Update pp for new mesh points. Ok as long as we also update geometry
        // and wanted displacement (usually zero if mesh motion has succeeded)
        ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
        patchDisp -= (ppLocalPoints-ppPtr().localPoints());

        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing smoothed mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
        }
    }


    // Construct iterative mesh mover.
    Info<< "Constructing mesh displacer ..." << endl;
    Info<< "Using mesh parameters " << motionDict << nl << endl;

    autoPtr<displacementMotionSolver> motionPtr
    (
        makeMotionSolver
        (
            pointMesh::New(mesh),
            snapDict,
            adaptPatchIDs
            //patchConstraints
        )
    );
    autoPtr<motionSmoother> meshMoverPtr
    (
        new motionSmoother
        (
            mesh,
            ppPtr(),
            adaptPatchIDs,
            motionPtr().pointDisplacement(),
            motionDict,
            dryRun_
        )
    );


    // Check initial mesh
    Info<< "Checking initial mesh ..." << endl;
    labelHashSet wrongFaces(mesh.nFaces()/100);
    motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces, dryRun_);
    const label nInitErrors = returnReduce
    (
        wrongFaces.size(),
        sumOp<label>()
    );

    Info<< "Detected " << nInitErrors << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;


    Info<< "Checked initial mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;



    //- Only if in feature attraction mode:
    // Point on nearest feature
    vectorField patchFeaturePoint(ppPtr().nPoints(), Zero);
    // Constraints at feature
    List<pointConstraint> patchConstraints(ppPtr().nPoints());

    for (label iter = 0; iter < nFeatIter; iter++)
    {
        Info<< nl
            << "Morph iteration " << iter << nl
            << "-----------------" << endl;


        // Calculate displacement at every patch point if we need it:
        // - if automatic near-surface detection
        // - if face splitting active
        pointField nearestPoint(ppPtr().nPoints(), vector::max);
        vectorField nearestNormal(ppPtr().nPoints(), Zero);


        const bool strictRegionSnap
        (
            iter < nFeatIter/2
          ? snapParams.strictRegionSnap()
          : Switch(true)
        );

        vectorField disp = calcNearestSurface
        (
            strictRegionSnap,               // attract points to region only
            meshRefiner_,
            globalToMasterPatch_,           // for if strictRegionSnap
            globalToSlavePatch_,            // for if strictRegionSnap
            ppPtr(),
            ppLocalPoints,
            snapDist,                       // max snap distance

            nearestPoint,
            nearestNormal
        );

        // Override displacement at thin gaps
        if (snapParams.detectNearSurfacesSnap())
        {
            detectNearSurfaces
            (
                Foam::cos(degToRad(planarAngle)),// planar cos for gaps
                ppPtr(),
                ppLocalPoints,
                nearestPoint,   // surfacepoint from nearest test
                nearestNormal,  // surfacenormal from nearest test

                disp
            );
        }

        // Override displacement with feature edge attempt
        if (doFeatures)
        {
            //- Any faces to split
            DynamicList<label> splitFaces;
            //- Indices in face to split across
            DynamicList<labelPair> splits;
            //- Patches for split face
            DynamicList<labelPair> splitPatches;

            // Offset to project to nearest feature. Use in combination with
            // patchConstraints.
            vectorField patchAttraction;

            disp = calcNearestSurfaceFeature
            (
                snapParams,
                false,              // no alignMeshEdges
                true,               // check >=3 patch points

                iter,
                featureCos,
                scalar(iter+1)/nFeatIter,

                snapDist,
                disp,               // nearest surface
                nearestNormal,
                ppPtr(),
                ppLocalPoints,

                patchAttraction,    // offset wrt ppLocalPoints to nearest
                                    // feature edge/point
                patchConstraints,   // feature type + constraint

                splitFaces,
                splits,
                splitPatches
            );


            // Freeze points on exposed points/faces
            freezeExposedPoints
            (
                meshRefiner_,
                "frozenFaces",      // faceZone name
                "frozenPoints",     // pointZone name
                ppPtr(),
                disp
            );

            patchFeaturePoint = ppLocalPoints+patchAttraction;

            if (debug)
            {
                OBJstream str
                (
                    mesh.time().path()
                  / "calcNearestSurfaceFeature"
                  + meshRefiner_.timeName()
                  + ".obj"
                );
                forAll(ppLocalPoints, pointi)
                {
                    const point& pt = ppLocalPoints[pointi];
                    str.write(linePointRef(pt, pt+disp[pointi]));
                }
            }



            // Split any faces:
            // - does not move/add any points
            // - tries to move new faces to correct patch

            if (returnReduce(splitFaces.size(), sumOp<label>()))
            {
                polyTopoChange meshMod(mesh);

                // Insert the mesh changes
                meshRefiner_.doSplitFaces
                (
                    splitFaces,
                    splits,
                    splitPatches,
                    meshMod
                );

                // Save old meshPoints before changing mesh
                const Map<label> oldMeshPointMap(ppPtr->meshPointMap());
                Pout<< "old pp points:" << ppPtr->nPoints()
                    << " oldMeshPointMap:" << oldMeshPointMap.size()
                    << endl;

                // Remove any unnecessary fields
                meshMoverPtr.clear();
                motionPtr.clear();
                ppPtr.clear();
                mesh.clearOut();
                mesh.moving(false);

                // Change the mesh (no inflation)
                autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
                mapPolyMesh& map = *mapPtr;

                // Update fields
                mesh.updateMesh(map);

                // Move mesh (since morphing might not do this)
                if (map.hasMotionPoints())
                {
                    mesh.movePoints(map.preMotionPoints());
                }
                else
                {
                    mesh.clearOut();
                }

                // Reset the instance for if in overwrite mode
                mesh.setInstance(meshRefiner_.timeName());
                meshRefiner_.setInstance(mesh.facesInstance());

                // Update intersections on split faces
                {
                    DynamicList<label> changedFaces(splitFaces.size());
                    Map<label> splitFacesMap(splitFaces.size());
                    forAll(splitFaces, i)
                    {
                        splitFacesMap.insert(splitFaces[i], i);
                    }
                    forAll(map.faceMap(), facei)
                    {
                        if (splitFacesMap.find(map.faceMap()[facei]))
                        {
                            changedFaces.append(facei);
                        }
                    }
                    // Update intersections on changed faces
                    meshRefiner_.updateMesh
                    (
                        map,
                        meshRefiner_.growFaceCellFace(changedFaces)
                    );
                }

                if (debug&meshRefinement::MESH)
                {
                    const_cast<Time&>(mesh.time())++;
                    Info<< "Writing split-faces mesh to time "
                        << meshRefiner_.timeName() << endl;
                    meshRefiner_.write
                    (
                        meshRefinement::debugType(debug),
                        meshRefinement::writeType
                        (
                            meshRefinement::writeLevel()
                          | meshRefinement::WRITEMESH
                        ),
                        mesh.time().path()/meshRefiner_.timeName()
                    );
                }



                // Update local mesh data
                // ~~~~~~~~~~~~~~~~~~~~~~

                Info<< "Updating for face-splitting" << endl;

                // baffles
                forAll(internalBaffles, i)
                {
                    labelPair& baffle = internalBaffles[i];
                    baffle.first() = map.reverseFaceMap()[baffle.first()];
                    baffle.second() = map.reverseFaceMap()[baffle.second()];
                }

                // re-do patch (since faces might have been split)
                ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
                motionPtr = makeMotionSolver
                (
                    pointMesh::New(mesh),
                    snapDict,
                    adaptPatchIDs
                    //patchConstraints
                );
                meshMoverPtr.reset
                (
                    new motionSmoother
                    (
                        mesh,
                        ppPtr(),
                        adaptPatchIDs,
                        motionPtr().pointDisplacement(),
                        motionDict,
                        dryRun_
                    )
                );

                const auto& mp = ppPtr->meshPoints();
                // pointMap (new-to-old) for pp points. Note: no points changed
                // but local point ordering might have changed since faces
                // split.
                labelList ppMap(mp.size());
                forAll(mp, i)
                {
                    ppMap[i] = oldMeshPointMap[mp[i]];
                }
                // patchDisp
                meshRefinement::updateList(ppMap, vector::zero, patchDisp);
                // snapDist
                meshRefinement::updateList(ppMap, scalar(0), snapDist);
                // patchFeaturePoint
                meshRefinement::updateList
                (
                    ppMap,
                    vector::zero,
                    patchFeaturePoint
                );
                // patchConstraints
                meshRefinement::updateList
                (
                    ppMap,
                    pointConstraint(),
                    patchConstraints
                );
//                // nearestPoint
//                meshRefinement::updateList(ppMap, vector::zero, nearestPoint);
//                // nearestNormal
//                meshRefinement::updateList(ppMap, vector::zero, nearestNormal);
                // disp
                meshRefinement::updateList(ppMap, vector::zero, disp);
                // ppLocalPoints
                meshRefinement::updateList(ppMap, vector::zero, ppLocalPoints);

                Info<< "DONE Updating for face-splitting" << endl;
            }
        }


        // Attract/slide the mesh a bit
        {
            //// Create pointMesh with the correct patches
            //// - points on single polyPatches stay as is
            //// - points on two polyPatches go to allEdgePatchName
            //// - points on >two polyPatches go to allPointPatchName
            //const pointMesh& pMesh = makePointMesh
            //(
            //    ppPtr(),
            //    patchConstraints,
            //    allEdgePatchName,
            //    allPointPatchName
            //);
            //const pointMesh& pMesh = pointMesh::New(mesh);
            //
            //autoPtr<displacementMotionSolver> motionPtr
            //(
            //    makeMotionSolver
            //    (
            //        pMesh,
            //        snapDict,
            //        adaptPatchIDs
            //        //patchConstraints
            //    )
            //);

            // Insert as bc to motionSmoother. Note that this takes displacement
            // relative to points0
            setDisplacement
            (
                ppPtr(),
                disp,
                adaptPatchIDs,
                motionPtr().points0(),
                motionPtr().pointDisplacement()
            );

            // Solve internal displacement
            tmp<pointField> tnewPoints(motionPtr->newPoints());

            // Move points
            if (false)
            {
                // 1. Directly move points
                mesh.movePoints(tnewPoints);
                // Optional? Only geometry used is ppLocalPoints which we keep
                // and update 'by hand'.
                //ppPtr().movePoints(tnewPoints);
            }
            else
            {
                // 2. Use motionSmoother
                // Set initial distribution of displacement field (on patches)
                // from patchDisp and make displacement consistent with b.c.
                // on displacement pointVectorField.
                const vectorField newDisp(tnewPoints()-meshMoverPtr().oldPoints());

                meshMoverPtr().displacement().vectorField::operator=(newDisp);
                meshMoverPtr().setDisplacement(disp);

                // Apply internal displacement to mesh.
                const bool meshOk = scaleMesh
                (
                    snapParams,
                    nInitErrors,
                    internalBaffles,
                    meshMoverPtr()
                );

                if (!meshOk)
                {
                    WarningInFunction
                        << "Did not successfully snap mesh."
                        << " Continuing to snap to resolve easy" << nl
                        << "    surfaces but the"
                        << " resulting mesh will not satisfy your quality"
                        << " constraints" << nl << endl;
                }

                // Use current mesh as base mesh
                meshMoverPtr().correct();
            }

            // Update pp for new mesh points. Ok as long as we also update
            // geometry and wanted displacement (usually zero if mesh motion
            // has succeeded)
            ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
            patchDisp -= (ppLocalPoints-ppPtr().localPoints());

            if (debug & meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing partially moved mesh to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    meshRefiner_.timeName()
                );
            }
        }

        // Split problematic cells
        {
            Info<< nl << "Checking moved mesh ..." << endl;
            faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/1000);
            motionSmoother::checkMesh
            (
                false,
                mesh,
                motionDict,
                identity(mesh.nFaces()),
                internalBaffles,
                wrongFaces,
                false           // dryRun_
            );
            const label nWrong(returnReduce(wrongFaces.size(), sumOp<label>()));
            Info<< "Detected " << nWrong
                << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;

            //if (nWrong)
            if (false)
            {
                bitSet decomposeCell(mesh.nCells());
                for (const label facei : wrongFaces)
                {
                    const label own = mesh.faceOwner()[facei];
                    if (!tetMatcher::test(mesh, own))
                    {
                        decomposeCell.set(own);
                    }
                    if (mesh.isInternalFace(facei))
                    {
                        const label nei = mesh.faceNeighbour()[facei];
                        if (!tetMatcher::test(mesh, nei))
                        {
                            decomposeCell.set(nei);
                        }
                    }
                }
                Pout<< "spliyyinG :" << decomposeCell.count() << " cells"
                    << endl;

                tetDecomposer tetDecomp(mesh);
                polyTopoChange meshMod(mesh);
                tetDecomp.setRefinement
                (
                    tetDecomposer::FACE_CENTRE_TRIS,
                    decomposeCell,
                    meshMod
                );

                // Save old meshPoints before changing mesh
                const Map<label> oldMeshPointMap(ppPtr->meshPointMap());

                Pout<< "old pp points:" << ppPtr->nPoints()
                    << " oldMeshPointMap:" << oldMeshPointMap.size()
                    << endl;

                // Remove any unnecessary fields
                meshMoverPtr.clear();
                motionPtr.clear();
                ppPtr.clear();
                mesh.clearOut();
                mesh.moving(false);

                // Change the mesh (no inflation)
                autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
                mapPolyMesh& map = *mapPtr;

                // Update fields
                mesh.updateMesh(map);

                // Move mesh (since morphing does not do this)
                if (map.hasMotionPoints())
                {
                    mesh.movePoints(map.preMotionPoints());
                }
                else
                {
                    mesh.clearOut();
                }

                // Reset the instance for if in overwrite mode
                mesh.setInstance(meshRefiner_.timeName());
                meshRefiner_.setInstance(mesh.facesInstance());

                //- Update numbering on tet-decomposition engine
                tetDecomp.updateMesh(map);

                bitSet isChangedFace(mesh.nFaces());
                forAll(map.cellMap(), celli)
                {
                    if (decomposeCell[map.cellMap()[celli]])
                    {
                        isChangedFace.set(mesh.cells()[celli]);
                    }
                }
                syncTools::syncFaceList
                (
                    mesh,
                    isChangedFace,
                    orEqOp<unsigned int>()
                );

                Pout<< "isChangedFace :" << decomposeCell.count() << " faces"
                    << endl;


                // Update intersection info
                meshRefiner_.updateMesh(map, isChangedFace.toc());


                if (debug&meshRefinement::MESH)
                {
                    const_cast<Time&>(mesh.time())++;
                    Info<< "Writing tet-decomp mesh to time "
                        << meshRefiner_.timeName() << endl;
                    meshRefiner_.write
                    (
                        meshRefinement::debugType(debug),
                        meshRefinement::writeType
                        (
                            meshRefinement::writeLevel()
                          | meshRefinement::WRITEMESH
                        ),
                        mesh.time().path()/meshRefiner_.timeName()
                    );
                }


                // Update local mesh data
                // ~~~~~~~~~~~~~~~~~~~~~~

                Info<< "Updating for tet-decomp" << endl;

                // baffles
                forAll(internalBaffles, i)
                {
                    labelPair& baffle = internalBaffles[i];
                    baffle.first() = map.reverseFaceMap()[baffle.first()];
                    baffle.second() = map.reverseFaceMap()[baffle.second()];
                }

                // re-do patch (since faces might have been split)
                ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
    Pout<< "new pp points:" << ppPtr->nPoints()
        << " new meshPointMap:" << ppPtr->meshPointMap().size()
        << endl;
                const auto& mp = ppPtr->meshPoints();
                // pointMap (new-to-old) for pp points. Might have new
                // face-centre points - these get mapped from any point on
                // the originating face.
                labelList ppMap(mp.size(), -1);
                forAll(mp, i)
                {
                    const label oldMeshPointi = map.pointMap()[mp[i]];
                    const auto mpFnd = oldMeshPointMap.find(oldMeshPointi);
                    if (mpFnd)
                    {
                        ppMap[i] = mpFnd();
                    }
                }

                // patchDisp
                meshRefinement::updateList(ppMap, vector::zero, patchDisp);
                // snapDist
                meshRefinement::updateList(ppMap, scalar(0), snapDist);
                // patchFeaturePoint
                meshRefinement::updateList
                (
                    ppMap,
                    vector::zero,
                    patchFeaturePoint
                );
                // patchConstraints
                meshRefinement::updateList
                (
                    ppMap,
                    pointConstraint(),
                    patchConstraints
                );
//                // nearestPoint
//                meshRefinement::updateList(ppMap, vector::zero, nearestPoint);
//                // nearestNormal
//                meshRefinement::updateList(ppMap, vector::zero, nearestNormal);
//                // disp
//                meshRefinement::updateList(ppMap, vector::zero, disp);

                // Maximum distance to attract to nearest feature on surface
                snapDist = calcSnapDistance(mesh, snapParams, ppPtr());
                ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());

                Info<< "DONE Updating for tet-decomp" << endl;
            }
        }
    }

//XXXXXX
/*
    {
        // Introduce single layer of cells. Straight from snappyLayerDriver

        // Global face indices engine
        const globalIndex globalFaces(mesh.nFaces());

        // Determine extrudePatch.edgeFaces in global numbering (so across
        // coupled patches). This is used only to string up edges
        // between coupled
        // faces (all edges between same (global)face indices get extruded).
        labelListList edgeGlobalFaces
        (
            addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                ppPtr()
            )
        );

        // Determine patches for extruded boundary edges. Calculates if any
        // additional processor patches need to be constructed.

        labelList edgePatchID;
        labelList edgeZoneID;
        boolList edgeFlip;
        labelList inflateFaceID;
        snappyLayerDriver::determineSidePatches
        (
            meshRefiner_,
            globalFaces,
            edgeGlobalFaces,
            ppPtr(),

            edgePatchID,
            edgeZoneID,
            edgeFlip,
            inflateFaceID
        );


        // Mesh topo change engine. Insert current mesh.
        polyTopoChange meshMod(mesh);

        // Layer mesh modifier
        addPatchCellLayer addLayer(mesh);

        // Extrude very thin layer of cells
        pointField extrusion(PatchTools::pointNormals(mesh, ppPtr()));
        const tmp<scalarField> thickness
        (
            wantedThickness
            (
                ppPtr(),
                1e-3        // cellSizeFraction
            )
        );
        extrusion *= thickness;

        // Add topo regardless of whether extrudeStatus is extruderemove.
        // Not add layer if patchDisp is zero.

        const label ppNFaces = ppPtr().size();
        const label ppNPoints = ppPtr().nPoints();

        addLayer.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            scalarField(ppNPoints, 1),  // expansion ratio
            ppPtr(),
            bitSet(ppPtr().size()),     // no flip

            edgePatchID,    // boundary patch for extruded boundary edges
            edgeZoneID,     // zone for extruded edges
            edgeFlip,
            inflateFaceID,


            labelList(0),   // exposed patchIDs, not used for adding layers
            labelList(ppNFaces, 1),
            labelList(ppNPoints, 1),
            extrusion,  //patchAttraction,

            meshMod
        );

        // Save old meshPoints before changing mesh
        const Map<label> oldMeshPointMap(ppPtr->meshPointMap());
        Pout<< "old pp points:" << ppPtr->nPoints()
            << " oldMeshPointMap:" << oldMeshPointMap.size()
            << endl;

        // Remove any unnecessary fields
        meshMoverPtr.clear();
        motionPtr.clear();
        ppPtr.clear();
        mesh.clearOut();
        mesh.moving(false);

        // Apply the stored topo changes to the current mesh.
        autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map.hasMotionPoints())
        {
            mesh.movePoints(map.preMotionPoints());
        }
        else
        {
            // Hack to remove meshPhi - mapped incorrectly. TBD.
            mesh.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh.setInstance(meshRefiner_.timeName());

        // Re-do the patch
        ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
    Pout<< "new pp points:" << ppPtr->nPoints()
        << " new meshPointMap:" << ppPtr->meshPointMap().size()
        << endl;
        // Map the old-patch-point data to the new patch points
        // pointMap (new-to-old) for new pp points
        labelList ppMap(ppPtr->nPoints(), -1);
        {
            const labelListList& oldAddedPoints = addLayer.addedPoints();
            forAll(oldAddedPoints, oldPatchPointi)
            {
                const label oldPointi = oldAddedPoints[oldPatchPointi].last();
                const label newPointi = map.reversePointMap()[oldPointi];
                const label newPatchPointi = ppPtr().meshPointMap()[newPointi];

                ppMap[newPatchPointi] = oldPatchPointi;
            }
        }

        // Update attraction
        meshRefinement::updateList(ppMap, vector::zero, patchDisp);
        // snapDist
        meshRefinement::updateList(ppMap, scalar(0), snapDist);
        // patchFeaturePoint
        meshRefinement::updateList
        (
            ppMap,
            vector::zero,
            patchFeaturePoint
        );
        // patchConstraints
        meshRefinement::updateList
        (
            ppMap,
            pointConstraint(),
            patchConstraints
        );
        // ppLocalPoints
        //meshRefinement::updateList(ppMap, vector::zero, ppLocalPoints);
        ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());


        // Update numbering on layer
        addLayer.updateMesh
        (
            map,
            identity(ppNFaces),
            identity(ppNPoints)
        );

        // Update intersections
        const labelListList addedCells(addLayer.addedCells());
        bitSet isChangedFace(mesh.nFaces());
        for (const labelList& faceToCells : addedCells)
        {
            for (const label celli : faceToCells)
            {
                isChangedFace.set(mesh.cells()[celli]);
            }
        }
        meshRefiner_.updateMesh(map, isChangedFace.toc());


        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing mesh-with-layer to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
        }
        if (debug)
        {
            OBJstream str
            (
                mesh.time().path()
              / "new_projection"
              + meshRefiner_.timeName()
              + ".obj"
            );
            forAll(ppLocalPoints, pointi)
            {
                const point& pt = ppLocalPoints[pointi];
                str.write(linePointRef(pt, patchFeaturePoint[pointi]));
            }
            Pout<< "** writing mapped attraction to " << str.name() << endl;
        }
    }
*/
//XXXX
    {
        // Layer mesh modifier
        addPatchCellLayer addLayer(mesh);

        // Save old mesh points (to construct the new-to-old local patch points)
        const Map<label> oldMeshPointMap(ppPtr->meshPointMap());

        meshMoverPtr.clear();
        motionPtr.clear();

        const pointField thickness
        (
            wantedThickness(ppPtr(), 1e-3)             //cellSizeFraction
          * PatchTools::pointNormals(mesh, ppPtr())
        );
        autoPtr<mapPolyMesh> mapPtr = addBufferLayers
        (
            ppPtr(),
            thickness,
            addLayer
        );

        // Re-do the patch
        ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
        
        // Map the old-patch-point data to the new patch points
        // pointMap (new-to-old) for new pp points
        labelList ppMap(ppPtr->nPoints(), -1);
        {
            const labelListList& addedPoints = addLayer.addedPoints();
            forAll(addedPoints, oldPatchPointi)
            {
                const label newPointi = addedPoints[oldPatchPointi].last();
                const label newPatchPointi = ppPtr().meshPointMap()[newPointi];
                ppMap[newPatchPointi] = oldPatchPointi;
            }
        }

        // Update attraction
        meshRefinement::updateList(ppMap, vector::zero, patchDisp);
        // snapDist
        meshRefinement::updateList(ppMap, scalar(0), snapDist);
        // patchFeaturePoint
        meshRefinement::updateList
        (
            ppMap,
            vector::zero,
            patchFeaturePoint
        );
        // patchConstraints
        meshRefinement::updateList
        (
            ppMap,
            pointConstraint(),
            patchConstraints
        );
        // ppLocalPoints
        ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
    }
//XXXXXX
    const bool snapToGeometry = true;
    if (snapToGeometry)
    {
        // Create pointMesh with the correct patches
        // - points on single polyPatches stay as is
        // - points on two polyPatches go to allEdgePatchName
        // - points on >two polyPatches go to allPointPatchName
        //const pointMesh& pMesh = makePointMesh
        //(
        //    ppPtr(),
        //    patchConstraints,
        //    allEdgePatchName,
        //    allPointPatchName
        //);
        const pointMesh& pMesh = pointMesh::New(mesh);

        autoPtr<displacementMotionSolver> motionPtr
        (
            makeMotionSolver
            (
                pMesh,
                snapDict,
                adaptPatchIDs
                //patchConstraints
            )
        );

        // Insert as bc to motionSmoother. Note that this takes displacement
        // relative to points0
        setDisplacement
        (
            ppPtr(),
            patchFeaturePoint-ppLocalPoints,
            adaptPatchIDs,
            motionPtr().points0(),
            motionPtr().pointDisplacement()
        );


        if (debug)
        {
            OBJstream str
            (
                mesh.time().path()
              / "buffer_layer_new_projection"
              + meshRefiner_.timeName()
              + ".obj"
            );
            forAll(ppLocalPoints, pointi)
            {
                const point& pt = ppLocalPoints[pointi];
                str.write(linePointRef(pt, patchFeaturePoint[pointi]));
            }
            Pout<< "** writing mapped attraction to " << str.name() << endl;
        }


        // Solve internal displacement
        tmp<pointField> tnewPoints(motionPtr->newPoints());

        // Move points
        mesh.movePoints(tnewPoints);

        // Update pp for new mesh points. Ok as long as we also update geometry
        // and wanted displacement (usually zero if mesh motion has succeeded)
        ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
        patchDisp -= (ppLocalPoints-ppPtr().localPoints());

        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing smoothed LAYER mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
        }
    }


    // Merge any introduced baffles (from faceZones of faceType 'internal')
    {
        autoPtr<mapPolyMesh> mapPtr = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false   // baffle zones
        );

        if (mapPtr.valid())
        {
            if (debug & meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing baffle-merged mesh to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    meshRefiner_.timeName()
                );
            }
        }
    }

    // Repatch faces according to nearest. Do not repatch baffle faces.
    {
        labelList duplicateFace(getInternalOrBaffleDuplicateFace());

        repatchToSurface(snapParams, adaptPatchIDs, duplicateFace);
    }

    if (debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
    }
}


// ************************************************************************* //
