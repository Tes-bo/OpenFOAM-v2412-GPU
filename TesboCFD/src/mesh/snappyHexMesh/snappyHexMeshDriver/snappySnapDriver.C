/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Description
    All to do with snapping to the surface

\*----------------------------------------------------------------------------*/

#include "snappySnapDriver.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "mapPolyMesh.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "mergePoints.H"
#include "snapParameters.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "unitConversion.H"
#include "localPointRegion.H"
#include "PatchTools.H"
#include "refinementFeatures.H"
#include "weightedPosition.H"
#include "profiling.H"
#include "addPatchCellLayer.H"
#include "displacementMotionSolver.H"
#include "snappyLayerDriver.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappySnapDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate geometrically collocated points, Requires bitSet to be
// sized and initialised!
Foam::label Foam::snappySnapDriver::getCollocatedPoints
(
    const scalar tol,
    const pointField& points,
    bitSet& isCollocatedPoint
)
{
    labelList pointMap;
    label nUnique = Foam::mergePoints
    (
        points,                         // points
        tol,                            // mergeTol
        false,                          // verbose
        pointMap
    );
    bool hasMerged = (nUnique < points.size());

    if (!returnReduceOr(hasMerged))
    {
        return 0;
    }

    // Determine which merged points are referenced more than once
    label nCollocated = 0;

    // Per old point the newPoint. Or -1 (not set yet) or -2 (already seen
    // twice)
    labelList firstOldPoint(nUnique, -1);
    forAll(pointMap, oldPointi)
    {
        label newPointi = pointMap[oldPointi];

        if (firstOldPoint[newPointi] == -1)
        {
            // First use of oldPointi. Store.
            firstOldPoint[newPointi] = oldPointi;
        }
        else if (firstOldPoint[newPointi] == -2)
        {
            // Third or more reference of oldPointi -> non-manifold
            isCollocatedPoint.set(oldPointi);
            nCollocated++;
        }
        else
        {
            // Second reference of oldPointi -> non-manifold
            isCollocatedPoint.set(firstOldPoint[newPointi]);
            nCollocated++;

            isCollocatedPoint.set(oldPointi);
            nCollocated++;

            // Mark with special value to save checking next time round
            firstOldPoint[newPointi] = -2;
        }
    }
    return returnReduce(nCollocated, sumOp<label>());
}


Foam::tmp<Foam::pointField> Foam::snappySnapDriver::smoothInternalDisplacement
(
    const meshRefinement& meshRefiner,
    const motionSmoother& meshMover
)
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const polyMesh& mesh = meshMover.mesh();

    // Get neighbour refinement
    const hexRef8& cutter = meshRefiner.meshCutter();
    const labelList& cellLevel = cutter.cellLevel();


    // Get the faces on the boundary
    bitSet isFront(mesh.nFaces(), pp.addressing());

    // Walk out from the surface a bit. Poor man's FaceCellWave.
    // Commented out for now - not sure if needed and if so how much
    //for (label iter = 0; iter < 2; iter++)
    //{
    //    bitSet newIsFront(mesh.nFaces());
    //
    //    forAll(isFront, facei)
    //    {
    //        if (isFront.test(facei))
    //        {
    //            label own = mesh.faceOwner()[facei];
    //            const cell& ownFaces = mesh.cells()[own];
    //            newIsFront.set(ownFaces);
    //
    //            if (mesh.isInternalFace(facei))
    //            {
    //                label nei = mesh.faceNeighbour()[facei];
    //                const cell& neiFaces = mesh.cells()[nei];
    //                newIsFront.set(neiFaces);
    //            }
    //        }
    //    }
    //
    //    syncTools::syncFaceList
    //    (
    //        mesh,
    //        newIsFront,
    //        orEqOp<unsigned int>()
    //    );
    //
    //    isFront = newIsFront;
    //}

    // Mark all points on faces
    //  - not on the boundary
    //  - inbetween differing refinement levels
    bitSet isMovingPoint(mesh.nPoints());

    label nInterface = 0;

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        label ownLevel = cellLevel[mesh.faceOwner()[facei]];
        label neiLevel = cellLevel[mesh.faceNeighbour()[facei]];

        if (!isFront.test(facei) && ownLevel != neiLevel)
        {
            const face& f = mesh.faces()[facei];
            isMovingPoint.set(f);

            ++nInterface;
        }
    }

    labelList neiCellLevel;
    syncTools::swapBoundaryCellList(mesh, cellLevel, neiCellLevel);

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        label ownLevel = cellLevel[mesh.faceOwner()[facei]];
        label neiLevel = neiCellLevel[facei-mesh.nInternalFaces()];

        if (!isFront.test(facei) && ownLevel != neiLevel)
        {
            const face& f = mesh.faces()[facei];
            isMovingPoint.set(f);

            ++nInterface;
        }
    }

    if (debug)
    {
        Info<< "Found " << returnReduce(nInterface, sumOp<label>())
            << " faces out of " << mesh.globalData().nTotalFaces()
            << " inbetween refinement regions." << endl;
    }

    // Make sure that points that are coupled to a moving point are marked
    // as well
    syncTools::syncPointList(mesh, isMovingPoint, maxEqOp<unsigned int>(), 0);

    // Unmark any point on the boundary. If we're doing zero iterations of
    // face-cell wave we might have coupled points not being unmarked.
    isMovingPoint.unset(pp.meshPoints());

    // Make sure that points that are coupled to meshPoints but not on a patch
    // are unmarked as well
    syncTools::syncPointList(mesh, isMovingPoint, minEqOp<unsigned int>(), 1);


    // Calculate average of connected cells
    Field<weightedPosition> sumLocation
    (
        mesh.nPoints(),
        pTraits<weightedPosition>::zero
    );

    forAll(isMovingPoint, pointi)
    {
        if (isMovingPoint.test(pointi))
        {
            const labelList& pCells = mesh.pointCells(pointi);

            sumLocation[pointi].first() = pCells.size();
            for (const label celli : pCells)
            {
                sumLocation[pointi].second() += mesh.cellCentres()[celli];
            }
        }
    }

    // Add coupled contributions
    weightedPosition::syncPoints(mesh, sumLocation);

    auto tdisplacement = tmp<pointField>::New(mesh.nPoints(), Zero);
    auto& displacement = tdisplacement.ref();

    label nAdapted = 0;

    forAll(displacement, pointi)
    {
        const weightedPosition& wp = sumLocation[pointi];
        if (mag(wp.first()) > VSMALL)
        {
            displacement[pointi] =
                wp.second()/wp.first()
              - mesh.points()[pointi];
            nAdapted++;
        }
    }

    Info<< "Smoothing " << returnReduce(nAdapted, sumOp<label>())
        << " points inbetween refinement regions."
        << endl;

    return tdisplacement;
}


// Calculate displacement as average of patch points.
Foam::tmp<Foam::pointField> Foam::snappySnapDriver::smoothPatchDisplacement
(
    const motionSmoother& meshMover,
    const List<labelPair>& baffles
)
{
    const indirectPrimitivePatch& pp = meshMover.patch();

    // Calculate geometrically non-manifold points on the patch to be moved.
    bitSet nonManifoldPoint(pp.nPoints());
    label nNonManifoldPoints = getCollocatedPoints
    (
        SMALL,
        pp.localPoints(),
        nonManifoldPoint
    );
    Info<< "Found " << nNonManifoldPoints << " non-manifold point(s)."
        << endl;


    // Average points
    // ~~~~~~~~~~~~~~

    // We determine three points:
    // - average of (centres of) connected patch faces
    // - average of (centres of) connected internal mesh faces
    // - as fallback: centre of any connected cell
    // so we can do something moderately sensible for non/manifold points.

    // Note: the averages are calculated properly parallel. This is
    // necessary to get the points shared by processors correct.


    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = pp.points();
    const polyMesh& mesh = meshMover.mesh();

    // Get labels of faces to count (master of coupled faces and baffle pairs)
    bitSet isMasterFace(syncTools::getMasterFaces(mesh));

    {
        forAll(baffles, i)
        {
            label f0 = baffles[i].first();
            label f1 = baffles[i].second();

            if (isMasterFace.test(f0))
            {
                // Make f1 a slave
                isMasterFace.unset(f1);
            }
            else if (isMasterFace.test(f1))
            {
                isMasterFace.unset(f0);
            }
            else
            {
                FatalErrorInFunction
                    << "Both sides of baffle consisting of faces " << f0
                    << " and " << f1 << " are already slave faces."
                    << abort(FatalError);
            }
        }
    }


    // Get average position of boundary face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Field<weightedPosition> avgBoundary
    (
        pointFaces.size(),
        pTraits<weightedPosition>::zero
    );
    {
        forAll(pointFaces, patchPointi)
        {
            const labelList& pFaces = pointFaces[patchPointi];

            forAll(pFaces, pfi)
            {
                label facei = pFaces[pfi];

                if (isMasterFace.test(pp.addressing()[facei]))
                {
                    avgBoundary[patchPointi].first() += 1.0;
                    avgBoundary[patchPointi].second() +=
                        pp[facei].centre(points);
                }
            }
        }

        // Add coupled contributions
        weightedPosition::syncPoints(mesh, pp.meshPoints(), avgBoundary);

        // Normalise
        forAll(avgBoundary, i)
        {
            // Note: what if there is no master boundary face?
            if (mag(avgBoundary[i].first()) > VSMALL)
            {
                avgBoundary[i].second() /= avgBoundary[i].first();
            }
        }
    }


    // Get average position of internal face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Field<weightedPosition> avgInternal;
    {
        Field<weightedPosition> globalSum
        (
            mesh.nPoints(),
            pTraits<weightedPosition>::zero
        );

        // Note: no use of pointFaces
        const faceList& faces = mesh.faces();

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            const face& f = faces[facei];
            const point& fc = mesh.faceCentres()[facei];

            forAll(f, fp)
            {
                weightedPosition& wp = globalSum[f[fp]];
                wp.first() += 1.0;
                wp.second() += fc;
            }
        }

        // Count coupled faces as internal ones (but only once)
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchi)
        {
            if
            (
                patches[patchi].coupled()
             && refCast<const coupledPolyPatch>(patches[patchi]).owner()
            )
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchi]);

                const vectorField::subField faceCentres = pp.faceCentres();

                forAll(pp, i)
                {
                    const face& f = pp[i];
                    const point& fc = faceCentres[i];

                    forAll(f, fp)
                    {
                        weightedPosition& wp = globalSum[f[fp]];
                        wp.first() += 1.0;
                        wp.second() += fc;
                    }
                }
            }
        }

        // Add coupled contributions
        weightedPosition::syncPoints(mesh, globalSum);

        avgInternal.setSize(meshPoints.size());

        forAll(avgInternal, patchPointi)
        {
            label meshPointi = meshPoints[patchPointi];
            const weightedPosition& wp = globalSum[meshPointi];

            avgInternal[patchPointi].first() = wp.first();
            if (mag(wp.first()) < VSMALL)
            {
                // Set to zero?
                avgInternal[patchPointi].second() = wp.second();
            }
            else
            {
                avgInternal[patchPointi].second() = wp.second()/wp.first();
            }
        }
    }


    // Precalculate any cell using mesh point (replacement of pointCells()[])
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList anyCell(mesh.nPoints(), -1);
    forAll(mesh.faceOwner(), facei)
    {
        label own = mesh.faceOwner()[facei];
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }


    // Displacement to calculate.
    auto tpatchDisp = tmp<pointField>::New(meshPoints.size(), Zero);
    auto& patchDisp = tpatchDisp.ref();

    forAll(pointFaces, i)
    {
        label meshPointi = meshPoints[i];
        const point& currentPos = pp.points()[meshPointi];

        // Now we have the two average points and their counts:
        //  avgBoundary and avgInternal
        // Do some blending between the two.
        // Note: the following section has some reasoning behind it but the
        // blending factors can be experimented with.

        const weightedPosition& internal = avgInternal[i];
        const weightedPosition& boundary = avgBoundary[i];

        point newPos;

        if (!nonManifoldPoint.test(i))
        {
            // Points that are manifold. Weight the internal and boundary
            // by their number of faces and blend with
            scalar internalBlend = 0.1;
            scalar blend = 0.1;

            point avgPos =
                (
                   internalBlend*internal.first()*internal.second()
                  +(1-internalBlend)*boundary.first()*boundary.second()
                )
              / (
                    internalBlend*internal.first()
                   +(1-internalBlend)*boundary.first()
                );

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else if (internal.first() == 0)
        {
            // Non-manifold without internal faces. Use any connected cell
            // as internal point instead. Use precalculated any cell to avoid
            // e.g. pointCells()[meshPointi][0]

            const point& cc = mesh.cellCentres()[anyCell[meshPointi]];

            scalar cellCBlend = 0.8;
            scalar blend = 0.1;

            point avgPos = (1-cellCBlend)*boundary.second() + cellCBlend*cc;

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else
        {
            // Non-manifold point with internal faces connected to them
            scalar internalBlend = 0.9;
            scalar blend = 0.1;

            point avgPos =
                internalBlend*internal.second()
              + (1-internalBlend)*boundary.second();

            newPos = (1-blend)*avgPos + blend*currentPos;
        }

        patchDisp[i] = newPos - currentPos;
    }

    return tpatchDisp;
}
//XXXXXXX
//Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avg
//(
//    const indirectPrimitivePatch& pp,
//    const pointField& localPoints
//)
//{
//    const labelListList& pointEdges = pp.pointEdges();
//    const edgeList& edges = pp.edges();
//
//    auto tavg = tmp<pointField>::New(pointEdges.size(), Zero);
//    auto& avg = tavg.ref();
//
//    forAll(pointEdges, verti)
//    {
//        vector& avgPos = avg[verti];
//
//        const labelList& pEdges = pointEdges[verti];
//
//        forAll(pEdges, myEdgei)
//        {
//            const edge& e = edges[pEdges[myEdgei]];
//
//            label otherVerti = e.otherVertex(verti);
//
//            avgPos += localPoints[otherVerti];
//        }
//
//        avgPos /= pEdges.size();
//    }
//    return tavg;
//}
//Foam::tmp<Foam::pointField>
//Foam::snappySnapDriver::smoothLambdaMuPatchDisplacement
//(
//    const motionSmoother& meshMover,
//    const List<labelPair>& baffles
//)
//{
//    const indirectPrimitivePatch& pp = meshMover.patch();
//    pointField newLocalPoints(pp.localPoints());
//
//    const label iters = 90;
//    const scalar lambda = 0.33;
//    const scalar mu = 0.34;
//
//    for (label iter = 0; iter < iters; iter++)
//    {
//        // Lambda
//        newLocalPoints =
//            (1 - lambda)*newLocalPoints
//          + lambda*avg(pp, newLocalPoints);
//
//        // Mu
//        newLocalPoints =
//            (1 + mu)*newLocalPoints
//          - mu*avg(pp, newLocalPoints);
//    }
//    return newLocalPoints-pp.localPoints();
//}
//XXXXXXX


Foam::tmp<Foam::scalarField> Foam::snappySnapDriver::edgePatchDist
(
    const pointMesh& pMesh,
    const indirectPrimitivePatch& pp
)
{
    const polyMesh& mesh = pMesh();

    // Set initial changed points to all the patch points
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppi)
    {
        wallInfo[ppi] = pointEdgePoint(pp.localPoints()[ppi], 0.0);
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        mesh,
        pp.meshPoints(),
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // max iterations
    );

    // Copy edge values into scalarField
    auto tedgeDist = tmp<scalarField>::New(mesh.nEdges());
    auto& edgeDist = tedgeDist.ref();

    forAll(allEdgeInfo, edgei)
    {
        edgeDist[edgei] = Foam::sqrt(allEdgeInfo[edgei].distSqr());
    }

    return tedgeDist;
}


void Foam::snappySnapDriver::dumpMove
(
    const fileName& fName,
    const pointField& meshPts,
    const pointField& surfPts
)
{
    // Dump direction of growth into file
    Info<< "Dumping move direction to " << fName << endl;

    OFstream nearestStream(fName);

    label verti = 0;

    forAll(meshPts, pti)
    {
        meshTools::writeOBJ(nearestStream, meshPts[pti]);
        verti++;

        meshTools::writeOBJ(nearestStream, surfPts[pti]);
        verti++;

        nearestStream<< "l " << verti-1 << ' ' << verti << nl;
    }
}


// Check whether all displacement vectors point outwards of patch. Return true
// if so.
bool Foam::snappySnapDriver::outwardsDisplacement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp
)
{
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        vector disp(patchDisp[pointi]);

        scalar magDisp = mag(disp);

        if (magDisp > SMALL)
        {
            disp /= magDisp;

            bool outwards = meshTools::visNormal(disp, faceNormals, pFaces);

            if (!outwards)
            {
                Warning<< "Displacement " << patchDisp[pointi]
                    << " at mesh point " << pp.meshPoints()[pointi]
                    << " coord " << pp.points()[pp.meshPoints()[pointi]]
                    << " points through the surrounding patch faces" << endl;
                return false;
            }
        }
        else
        {
            //? Displacement small but in wrong direction. Would probably be ok.
        }
    }
    return true;
}


void Foam::snappySnapDriver::freezeExposedPoints
(
    const meshRefinement& meshRefiner,
    const word& fzName,     // faceZone name
    const word& pzName,     // pointZone name
    const indirectPrimitivePatch& outside,
    vectorField& outsideDisp
)
{
    const fvMesh& mesh = meshRefiner.mesh();
    const pointZoneMesh& pointZones = mesh.pointZones();

    bitSet isFrozenPoint(mesh.nPoints());

    // Add frozen points
    const label pointZonei = pointZones.findZoneID(pzName);
    if (pointZonei != -1)
    {
        isFrozenPoint.set(pointZones[pointZonei]);
    }

    // Add (inside) points of frozen faces
    const faceZoneMesh& faceZones = mesh.faceZones();
    const label faceZonei = faceZones.findZoneID(fzName);
    if (faceZonei != -1)
    {
        const uindirectPrimitivePatch pp
        (
            UIndirectList<face>(mesh.faces(), faceZones[faceZonei]),
            mesh.points()
        );

        // Count number of faces per edge
        const labelList nEdgeFaces(meshRefiner.countEdgeFaces(pp));

        // Freeze all internal points
        forAll(nEdgeFaces, edgei)
        {
            if (nEdgeFaces[edgei] != 1)
            {
                const edge& e = pp.edges()[edgei];
                isFrozenPoint.set(pp.meshPoints()[e[0]]);
                isFrozenPoint.set(pp.meshPoints()[e[1]]);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        isFrozenPoint,
        orEqOp<unsigned int>(),
        0u
    );

    for (const label pointi : isFrozenPoint)
    {
        const auto iter = outside.meshPointMap().find(pointi);
        if (iter.good())
        {
            outsideDisp[iter.val()] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappySnapDriver::snappySnapDriver
(
    meshRefinement& meshRefiner,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const bool dryRun
)
:
    meshRefiner_(meshRefiner),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch),
    dryRun_(dryRun)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::snappySnapDriver::calcSnapDistance
(
    const fvMesh& mesh,
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp
)
{
    const edgeList& edges = pp.edges();
    const labelListList& pointEdges = pp.pointEdges();
    const pointField& localPoints = pp.localPoints();

    scalarField maxEdgeLen(localPoints.size(), -GREAT);

    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];

        forAll(pEdges, pEdgei)
        {
            const edge& e = edges[pEdges[pEdgei]];

            scalar len = e.mag(localPoints);

            maxEdgeLen[pointi] = max(maxEdgeLen[pointi], len);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxEdgeLen,
        maxEqOp<scalar>(),  // combine op
        -GREAT              // null value
    );

    return scalarField(snapParams.snapTol()*maxEdgeLen);
}


void Foam::snappySnapDriver::preSmoothPatch
(
    const meshRefinement& meshRefiner,
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    addProfiling(smooth, "snappyHexMesh::snap::smoothing");
    const fvMesh& mesh = meshRefiner.mesh();

    labelList checkFaces;

    if (snapParams.nSmoothInternal() > 0)
    {
        Info<< "Smoothing patch and internal points ..." << endl;
    }
    else
    {
        Info<< "Smoothing patch points ..." << endl;
    }

    vectorField& pointDisp = meshMover.pointDisplacement().primitiveFieldRef();

    for
    (
        label smoothIter = 0;
        smoothIter < snapParams.nSmoothPatch();
        smoothIter++
    )
    {
        Info<< "Smoothing iteration " << smoothIter << endl;
        checkFaces.setSize(mesh.nFaces());
        forAll(checkFaces, facei)
        {
            checkFaces[facei] = facei;
        }

        // If enabled smooth the internal points
        if (snapParams.nSmoothInternal() > smoothIter)
        {
            // Override values on internal points on refinement interfaces
            pointDisp = smoothInternalDisplacement(meshRefiner, meshMover);
        }

        // Smooth the patch points
        pointField patchDisp(smoothPatchDisplacement(meshMover, baffles));
        //pointField patchDisp
        //(
        //  smoothLambdaMuPatchDisplacement(meshMover, baffles)
        //);

        // Take over patch displacement as boundary condition on
        // pointDisplacement
        meshMover.setDisplacement(patchDisp);

        // Start off from current mesh.points()
        meshMover.correct();

        scalar oldErrorReduction = -1;

        for (label snapIter = 0; snapIter < 2*snapParams.nSnap(); snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if (snapIter == snapParams.nSnap())
            {
                Info<< "Displacement scaling for error reduction set to 0."
                    << endl;
                oldErrorReduction = meshMover.setErrorReduction(0.0);
            }

            // Try to adapt mesh to obtain displacement by smoothly
            // decreasing displacement at error locations.
            if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
        }

        if (oldErrorReduction >= 0)
        {
            meshMover.setErrorReduction(oldErrorReduction);
        }
        Info<< endl;
    }


    // The current mesh is the starting mesh to smooth from.
    meshMover.correct();

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing patch smoothed mesh to time "
            << meshRefiner.timeName() << '.' << endl;
        meshRefiner.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner.timeName()
        );
        Info<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }

    Info<< "Patch points smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// Get (pp-local) indices of points that are both on zone and on patched surface
void Foam::snappySnapDriver::getZoneSurfacePoints
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp,
    const word& zoneName,

    bitSet& pointOnZone
)
{
    label zonei = mesh.faceZones().findZoneID(zoneName);

    if (zonei == -1)
    {
        FatalErrorInFunction
            << "Cannot find zone " << zoneName
            << exit(FatalError);
    }

    const faceZone& fZone = mesh.faceZones()[zonei];


    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    forAll(fZone, i)
    {
        const face& f = mesh.faces()[fZone[i]];

        forAll(f, fp)
        {
            label meshPointi = f[fp];

            const auto iter = pp.meshPointMap().cfind(meshPointi);

            if (iter.good())
            {
                const label pointi = iter.val();
                pointOnZone[pointi] = true;
            }
        }
    }
}


Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avgCellCentres
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp
)
{
    const labelListList& pointFaces = pp.pointFaces();

    Field<weightedPosition> avgBoundary
    (
        pointFaces.size(),
        pTraits<weightedPosition>::zero
    );

    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        avgBoundary[pointi].first() = pFaces.size();
        forAll(pFaces, pfi)
        {
            label facei = pFaces[pfi];
            label own = mesh.faceOwner()[pp.addressing()[facei]];
            avgBoundary[pointi].second() += mesh.cellCentres()[own];
        }
    }

    // Add coupled contributions
    weightedPosition::syncPoints(mesh, pp.meshPoints(), avgBoundary);

    auto tavgBoundary = tmp<pointField>::New(avgBoundary.size());
    weightedPosition::getPoints(avgBoundary, tavgBoundary.ref());

    return tavgBoundary;
}


//Foam::tmp<Foam::scalarField> Foam::snappySnapDriver::calcEdgeLen
//(
//    const indirectPrimitivePatch& pp
//) const
//{
//    // Get local edge length based on refinement level
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // (Ripped from snappyLayerDriver)
//
//    auto tedgeLen = tmp<scalarField>::New(pp.nPoints());
//    auto& edgeLen = tedgeLen.ref();
//    {
//        const fvMesh& mesh = meshRefiner_.mesh();
//        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
//        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
//
//        labelList maxPointLevel(pp.nPoints(), labelMin);
//
//        forAll(pp, i)
//        {
//            label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];
//            const face& f = pp.localFaces()[i];
//            forAll(f, fp)
//            {
//                maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
//            }
//        }
//
//        syncTools::syncPointList
//        (
//            mesh,
//            pp.meshPoints(),
//            maxPointLevel,
//            maxEqOp<label>(),
//            labelMin            // null value
//        );
//
//
//        forAll(maxPointLevel, pointi)
//        {
//            // Find undistorted edge size for this level.
//            edgeLen[pointi] = edge0Len/(1<<maxPointLevel[pointi]);
//        }
//    }
//    return tedgeLen;
//}


void Foam::snappySnapDriver::detectNearSurfaces
(
    const scalar planarCos,
    const indirectPrimitivePatch& pp,
    const pointField& localPoints,
    const pointField& nearestPoint,
    const vectorField& nearestNormal,

    vectorField& disp
) const
{
    Info<< "Detecting near surfaces ..." << endl;

    const labelList& meshPoints = pp.meshPoints();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const fvMesh& mesh = meshRefiner_.mesh();

    //// Get local edge length based on refinement level
    //const scalarField edgeLen(calcEdgeLen(pp));
    //
    //// Generate rays for every surface point
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //{
    //    const vector n = normalised(vector::one);
    //
    //    pointField start(14*pp.nPoints());
    //    pointField end(start.size());
    //
    //    label rayi = 0;
    //    forAll(localPoints, pointi)
    //    {
    //        const point& pt = localPoints[pointi];
    //
    //        // Along coordinate axes
    //
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.y() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.y() += edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.z() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.z() += edgeLen[pointi];
    //        }
    //
    //        // At 45 degrees
    //
    //        const vector vec(edgeLen[pointi]*n);
    //
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //    }
    //
    //    labelList surface1;
    //    List<pointIndexHit> hit1;
    //    labelList region1;
    //    vectorField normal1;
    //
    //    labelList surface2;
    //    List<pointIndexHit> hit2;
    //    labelList region2;
    //    vectorField normal2;
    //    surfaces.findNearestIntersection
    //    (
    //        unzonedSurfaces,    // surfacesToTest,
    //        start,
    //        end,
    //
    //        surface1,
    //        hit1,
    //        region1,
    //        normal1,
    //
    //        surface2,
    //        hit2,
    //        region2,
    //        normal2
    //    );
    //
    //    // All intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "surfaceHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with rays to " << str.name()
    //            << endl;
    //
    //        forAll(hit1, i)
    //        {
    //            if (hit1[i].hit())
    //            {
    //                str.writeLine(start[i], hit1[i].point());
    //            }
    //            if (hit2[i].hit())
    //            {
    //                str.writeLine(start[i], hit2[i].point());
    //            }
    //        }
    //    }
    //
    //    // Co-planar intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "coplanarHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with co-planar surfaces to "
    //            << str.name() << endl;
    //
    //        forAll(localPoints, pointi)
    //        {
    //            bool hasNormal = false;
    //            point surfPointA;
    //            vector surfNormalA;
    //            point surfPointB;
    //            vector surfNormalB;
    //
    //            bool isCoplanar = false;
    //
    //            label rayi = 14*pointi;
    //            for (label i = 0; i < 14; i++)
    //            {
    //                if (hit1[rayi].hit())
    //                {
    //                    const point& pt = hit1[rayi].point();
    //                    const vector& n = normal1[rayi];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //                if (hit2[rayi].hit())
    //                {
    //                    const point& pt = hit2[rayi].point();
    //                    const vector& n = normal2[rayi];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //
    //                rayi++;
    //            }
    //
    //            if (isCoplanar)
    //            {
    //                str.writeLine(surfPointA, surfPointB);
    //            }
    //        }
    //    }
    //}


    const pointField avgCc(avgCellCentres(mesh, pp));

    // Construct rays through localPoints to beyond cell centre
    pointField start(pp.nPoints());
    pointField end(pp.nPoints());
    forAll(localPoints, pointi)
    {
        const point& pt = localPoints[pointi];
        const vector d = 2*(avgCc[pointi]-pt);
        start[pointi] = pt - d;
        end[pointi] = pt + d;
    }


    autoPtr<OBJstream> gapStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        gapStr.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "detectNearSurfaces_" + meshRefiner_.timeName() + ".obj"
            )
        );
    }


    const bitSet isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            meshPoints
        )
    );

    label nOverride = 0;

    // 1. All points to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        const labelList unzonedSurfaces =
            surfaceZonesInfo::getUnnamedSurfaces
            (
                meshRefiner_.surfaces().surfZones()
            );

        // Do intersection test
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;
        surfaces.findNearestIntersection
        (
            unzonedSurfaces,
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );


        forAll(localPoints, pointi)
        {
            // Current location
            const point& pt = localPoints[pointi];

            bool override = false;

            //if (hit1[pointi].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointi],
            //            nearestNormal[pointi],
            //            hit1[pointi].point(),
            //            normal1[pointi]
            //        )
            //    )
            //    {
            //        disp[pointi] = hit1[pointi].point()-pt;
            //        override = true;
            //    }
            //}
            //if (hit2[pointi].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointi],
            //            nearestNormal[pointi],
            //            hit2[pointi].point(),
            //            normal2[pointi]
            //        )
            //    )
            //    {
            //        disp[pointi] = hit2[pointi].point()-pt;
            //        override = true;
            //    }
            //}

            if (hit1[pointi].hit() && hit2[pointi].hit())
            {
                if
                (
                    meshRefiner_.isGap
                    (
                        planarCos,
                        hit1[pointi].point(),
                        normal1[pointi],
                        hit2[pointi].point(),
                        normal2[pointi]
                    )
                )
                {
                    // TBD: check if the attraction (to nearest) would attract
                    // good enough and not override attraction

                    if (gapStr)
                    {
                        gapStr().writeLine(pt, hit2[pointi].point());
                    }

                    // Choose hit2 : nearest to end point (so inside the domain)
                    disp[pointi] = hit2[pointi].point()-pt;
                    override = true;
                }
            }

            if (override && isPatchMasterPoint[pointi])
            {
                nOverride++;
            }
        }
    }


    // 2. All points on zones to their respective surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Surfaces with zone information
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

        const labelList zonedSurfaces = surfaceZonesInfo::getNamedSurfaces
        (
            surfZones
        );

        forAll(zonedSurfaces, i)
        {
            label zoneSurfi = zonedSurfaces[i];
            const labelList surfacesToTest(1, zoneSurfi);

            const wordList& faceZoneNames =
                surfZones[zoneSurfi].faceZoneNames();
            forAll(faceZoneNames, namei)
            {
                const word& faceZoneName = faceZoneNames[namei];

                // Get indices of points both on faceZone and on pp.
                bitSet pointOnZone(pp.nPoints());
                getZoneSurfacePoints
                (
                    mesh,
                    pp,
                    faceZoneName,
                    pointOnZone
                );
                const labelList zonePointIndices(pointOnZone.toc());

                // Do intersection test
                labelList surface1;
                List<pointIndexHit> hit1;
                labelList region1;
                vectorField normal1;

                labelList surface2;
                List<pointIndexHit> hit2;
                labelList region2;
                vectorField normal2;
                surfaces.findNearestIntersection
                (
                    surfacesToTest,
                    pointField(start, zonePointIndices),
                    pointField(end, zonePointIndices),

                    surface1,
                    hit1,
                    region1,
                    normal1,

                    surface2,
                    hit2,
                    region2,
                    normal2
                );


                forAll(hit1, i)
                {
                    label pointi = zonePointIndices[i];

                    // Current location
                    const point& pt = localPoints[pointi];

                    bool override = false;

                    //if (hit1[i].hit())
                    //{
                    //    if
                    //    (
                    //        meshRefiner_.isGap
                    //        (
                    //            planarCos,
                    //            nearestPoint[pointi],
                    //            nearestNormal[pointi],
                    //            hit1[i].point(),
                    //            normal1[i]
                    //        )
                    //    )
                    //    {
                    //        disp[pointi] = hit1[i].point()-pt;
                    //        override = true;
                    //    }
                    //}
                    //if (hit2[i].hit())
                    //{
                    //    if
                    //    (
                    //        meshRefiner_.isGap
                    //        (
                    //            planarCos,
                    //            nearestPoint[pointi],
                    //            nearestNormal[pointi],
                    //            hit2[i].point(),
                    //            normal2[i]
                    //        )
                    //    )
                    //    {
                    //        disp[pointi] = hit2[i].point()-pt;
                    //        override = true;
                    //    }
                    //}

                    if (hit1[i].hit() && hit2[i].hit())
                    {
                        if
                        (
                            meshRefiner_.isGap
                            (
                                planarCos,
                                hit1[i].point(),
                                normal1[i],
                                hit2[i].point(),
                                normal2[i]
                            )
                        )
                        {
                            if (gapStr)
                            {
                                gapStr().writeLine(pt, hit2[i].point());
                            }

                            disp[pointi] = hit2[i].point()-pt;
                            override = true;
                        }
                    }

                    if (override && isPatchMasterPoint[pointi])
                    {
                        nOverride++;
                    }
                }
            }
        }
    }

    Info<< "Overriding nearest with intersection of close gaps at "
        << returnReduce(nOverride, sumOp<label>())
        << " out of " << returnReduce(pp.nPoints(), sumOp<label>())
        << " points." << endl;
}


void Foam::snappySnapDriver::calcNearestSurface
(
    const refinementSurfaces& surfaces,

    const labelList& surfacesToTest,
    const labelListList& regionsToTest,

    const pointField& localPoints,
    const labelList& zonePointIndices,

    scalarField& minSnapDist,
    labelList& snapSurf,
    vectorField& patchDisp,

    // Optional: nearest point, normal
    pointField& nearestPoint,
    vectorField& nearestNormal
)
{
    // Find nearest for points both on faceZone and pp.
    List<pointIndexHit> hitInfo;
    labelList hitSurface;

    if (nearestNormal.size() == localPoints.size())
    {
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            surfacesToTest,
            regionsToTest,

            pointField(localPoints, zonePointIndices),
            sqr(scalarField(minSnapDist, zonePointIndices)),

            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                label pointi = zonePointIndices[i];
                nearestPoint[pointi] = hitInfo[i].point();
                nearestNormal[pointi] = hitNormal[i];
            }
        }
    }
    else
    {
        surfaces.findNearest
        (
            surfacesToTest,
            regionsToTest,

            pointField(localPoints, zonePointIndices),
            sqr(scalarField(minSnapDist, zonePointIndices)),

            hitSurface,
            hitInfo
        );
    }

    forAll(hitInfo, i)
    {
        if (hitInfo[i].hit())
        {
            label pointi = zonePointIndices[i];

            patchDisp[pointi] = hitInfo[i].point() - localPoints[pointi];
            minSnapDist[pointi] = mag(patchDisp[pointi]);
            snapSurf[pointi] = hitSurface[i];
        }
    }
}


Foam::vectorField Foam::snappySnapDriver::calcNearestSurface
(
    const bool strictRegionSnap,
    const meshRefinement& meshRefiner,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const indirectPrimitivePatch& pp,
    const pointField& localPoints,
    const scalarField& snapDist,
    pointField& nearestPoint,
    vectorField& nearestNormal
)
{
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;
    if (strictRegionSnap)
    {
        Info<< "    non-zone points : attract to local region on surface only"
            << nl
            << "    zone points     : attract to local region on surface only"
            << nl
            << endl;
    }
    else
    {
        Info<< "    non-zone points :"
            << " attract to nearest of all non-zone surfaces"
            << nl
            << "    zone points     : attract to zone surface only" << nl
            << endl;
    }


    const refinementSurfaces& surfaces = meshRefiner.surfaces();
    const fvMesh& mesh = meshRefiner.mesh();

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), Zero);

    if (returnReduceOr(localPoints.size()))
    {
        // Current surface snapped to. Used to check whether points have been
        // snapped at all
        labelList snapSurf(localPoints.size(), -1);

        // Current best snap distance (since point might be on multiple
        // regions)
        scalarField minSnapDist(snapDist);


        if (strictRegionSnap)
        {
            // Attract patch points to same region only

            forAll(surfaces.surfaces(), surfi)
            {
                label geomi = surfaces.surfaces()[surfi];
                label nRegions = surfaces.geometry()[geomi].regions().size();

                const labelList surfacesToTest(1, surfi);

                for (label regioni = 0; regioni < nRegions; regioni++)
                {
                    label globali = surfaces.globalRegion(surfi, regioni);
                    label masterPatchi = globalToMasterPatch[globali];

                    // Get indices of points both on patch and on pp
                    labelList zonePointIndices
                    (
                        getFacePoints
                        (
                            pp,
                            mesh.boundaryMesh()[masterPatchi]
                        )
                    );

                    calcNearestSurface
                    (
                        surfaces,

                        surfacesToTest,
                        labelListList(1, labelList(1, regioni)), //regionsToTest

                        localPoints,
                        zonePointIndices,

                        minSnapDist,
                        snapSurf,
                        patchDisp,

                        // Optional: nearest point, normal
                        nearestPoint,
                        nearestNormal
                    );

                    if (globalToSlavePatch[globali] != masterPatchi)
                    {
                        label slavePatchi = globalToSlavePatch[globali];

                        // Get indices of points both on patch and on pp
                        labelList zonePointIndices
                        (
                            getFacePoints
                            (
                                pp,
                                mesh.boundaryMesh()[slavePatchi]
                            )
                        );

                        calcNearestSurface
                        (
                            surfaces,

                            surfacesToTest,
                            labelListList(1, labelList(1, regioni)),

                            localPoints,
                            zonePointIndices,

                            minSnapDist,
                            snapSurf,
                            patchDisp,

                            // Optional: nearest point, normal
                            nearestPoint,
                            nearestNormal
                        );
                    }
                }
            }
        }
        else
        {
            // Divide surfaces into zoned and unzoned
            const labelList unzonedSurfaces =
                surfaceZonesInfo::getUnnamedSurfaces
                (
                    meshRefiner.surfaces().surfZones()
                );


            // 1. All points to non-interface surfaces
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            List<pointIndexHit> hitInfo;
            labelList hitSurface;

            if (nearestNormal.size() == localPoints.size())
            {
                labelList hitRegion;
                vectorField hitNormal;
                surfaces.findNearestRegion
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),
                    hitSurface,
                    hitInfo,
                    hitRegion,
                    hitNormal
                );

                forAll(hitInfo, pointi)
                {
                    if (hitInfo[pointi].hit())
                    {
                        nearestPoint[pointi] = hitInfo[pointi].point();
                        nearestNormal[pointi] = hitNormal[pointi];
                    }
                }
            }
            else
            {
                surfaces.findNearest
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),        // sqr of attract distance
                    hitSurface,
                    hitInfo
                );
            }

            forAll(hitInfo, pointi)
            {
                if (hitInfo[pointi].hit())
                {
                    patchDisp[pointi] =
                        hitInfo[pointi].point()
                      - localPoints[pointi];

                    snapSurf[pointi] = hitSurface[pointi];
                }
            }


            const labelList zonedSurfaces = surfaceZonesInfo::getNamedSurfaces
            (
                meshRefiner.surfaces().surfZones()
            );


            // 2. All points on zones to their respective surface
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // (ignoring faceZone subdivision)

            // Surfaces with zone information
            const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

            forAll(zonedSurfaces, i)
            {
                label surfi = zonedSurfaces[i];
                const labelList surfacesToTest(1, surfi);
                const label geomi = surfaces.surfaces()[surfi];
                const label nRegions =
                    surfaces.geometry()[geomi].regions().size();

                const wordList& faceZoneNames =
                    surfZones[surfi].faceZoneNames();

                // Get indices of points both on any faceZone and on pp.
                bitSet pointOnZone(pp.nPoints());
                forAll(faceZoneNames, locali)
                {
                    getZoneSurfacePoints
                    (
                        mesh,
                        pp,
                        faceZoneNames[locali],
                        pointOnZone
                    );
                }
                const labelList zonePointIndices(pointOnZone.toc());

                calcNearestSurface
                (
                    surfaces,

                    surfacesToTest,
                    labelListList(1, identity(nRegions)),

                    localPoints,
                    zonePointIndices,

                    minSnapDist,
                    snapSurf,
                    patchDisp,

                    // Optional: nearest point, normal
                    nearestPoint,
                    nearestNormal
                );
            }
        }


        // Check if all points are being snapped
        forAll(snapSurf, pointi)
        {
            if (snapSurf[pointi] == -1)
            {
                static label nWarn = 0;

                if (nWarn < 100)
                {
                    WarningInFunction
                        << "For point:" << pointi
                        << " coordinate:" << localPoints[pointi]
                        << " did not find any surface within:"
                        << minSnapDist[pointi] << " metre." << endl;
                    nWarn++;
                    if (nWarn == 100)
                    {
                        WarningInFunction
                            << "Reached warning limit " << nWarn
                            << ". Suppressing further warnings." << endl;
                    }
                }
            }
        }

        {
            const bitSet isPatchMasterPoint
            (
                meshRefinement::getMasterPoints
                (
                    mesh,
                    pp.meshPoints()
                )
            );

            scalarField magDisp(mag(patchDisp));

            Info<< "Wanted displacement : average:"
                <<  meshRefinement::gAverage(isPatchMasterPoint, magDisp)
                << " min:" << gMin(magDisp)
                << " max:" << gMax(magDisp) << endl;
        }
    }

    Info<< "Calculated surface displacement in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement. Can not happen for triSurfaceMesh but
    // can happen for some analytical shapes?
    forAll(patchDisp, patchPointi)
    {
        scalar magDisp = mag(patchDisp[patchPointi]);

        if (magDisp > snapDist[patchPointi])
        {
            patchDisp[patchPointi] *= snapDist[patchPointi] / magDisp;

            Pout<< "Limiting displacement for " << patchPointi
                << " from " << magDisp << " to " << snapDist[patchPointi]
                << endl;
        }
    }

    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cannot use VGREAT)
    );

    return patchDisp;
}


void Foam::snappySnapDriver::smoothDisplacement
(
    const snapParameters& snapParams,
    motionSmoother& meshMover
) const
{
    if (dryRun_)
    {
        return;
    }

    const fvMesh& mesh = meshRefiner_.mesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    Info<< "Smoothing displacement ..." << endl;

    // Set edge diffusivity as inverse of distance to patch
    scalarField edgeGamma(1.0/(edgePatchDist(meshMover.pMesh(), pp) + SMALL));
    //scalarField edgeGamma(mesh.nEdges(), 1.0);
    //scalarField edgeGamma(wallGamma(mesh, pp, 10, 1));

    // Get displacement field
    pointVectorField& disp = meshMover.displacement();

    for (label iter = 0; iter < snapParams.nSmoothDispl(); iter++)
    {
        if ((iter % 10) == 0)
        {
            Info<< "Iteration " << iter << endl;
        }
        pointVectorField oldDisp(disp);
        meshMover.smooth(oldDisp, edgeGamma, disp);
    }
    Info<< "Displacement smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing smoothed mesh to time " << meshRefiner_.timeName()
            << endl;

        // Moving mesh creates meshPhi. Can be cleared out by a mesh.clearOut
        // but this will also delete all pointMesh but not pointFields which
        // gives an illegal situation.

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
        Info<< "Writing displacement field ..." << endl;
        disp.write();
        tmp<pointScalarField> magDisp(mag(disp));
        magDisp().write();

        Info<< "Writing actual patch displacement ..." << endl;
        vectorField actualPatchDisp(disp, pp.meshPoints());
        dumpMove
        (
            mesh.time().path()
          / "actualPatchDisplacement_" + meshRefiner_.timeName() + ".obj",
            pp.localPoints(),
            pp.localPoints() + actualPatchDisp
        );
    }
}


bool Foam::snappySnapDriver::scaleMesh
(
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    addProfiling(scale, "snappyHexMesh::snap::scale");
    const fvMesh& mesh = meshRefiner_.mesh();

    // Relax displacement until correct mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList checkFaces(identity(mesh.nFaces()));

    scalar oldErrorReduction = -1;

    bool meshOk = false;

    Info<< "Moving mesh ..." << endl;
    for (label iter = 0; iter < 2*snapParams.nSnap(); iter++)
    {
        Info<< nl << "Iteration " << iter << endl;

        if (iter == snapParams.nSnap())
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        meshOk = meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors);

        if (meshOk)
        {
            Info<< "Successfully moved mesh" << endl;
            break;
        }
        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing scaled mesh to time " << meshRefiner_.timeName()
                << endl;
            mesh.write();

            Info<< "Writing displacement field ..." << endl;
            meshMover.displacement().write();
            tmp<pointScalarField> magDisp(mag(meshMover.displacement()));
            magDisp().write();
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }
    Info<< "Moved mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    return meshOk;
}


// After snapping: correct patching according to nearest surface.
// Code is very similar to calcNearestSurface.
// - calculate face-wise snap distance as max of point-wise
// - calculate face-wise nearest surface point
// - repatch face according to patch for surface point.
Foam::autoPtr<Foam::mapPolyMesh> Foam::snappySnapDriver::repatchToSurface
(
    const snapParameters& snapParams,
    const labelList& adaptPatchIDs,
    const labelList& preserveFaces
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    Info<< "Repatching faces according to nearest surface ..." << endl;

    // Get the labels of added patches.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones());
    labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones());


    // Faces that do not move
    bitSet isZonedFace(mesh.nFaces());
    {
        // 1. Preserve faces in preserveFaces list
        forAll(preserveFaces, facei)
        {
            if (preserveFaces[facei] != -1)
            {
                isZonedFace.set(facei);
            }
        }

        // 2. All faces on zoned surfaces
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();
        const faceZoneMesh& fZones = mesh.faceZones();

        forAll(zonedSurfaces, i)
        {
            const label zoneSurfi = zonedSurfaces[i];
            const wordList& fZoneNames = surfZones[zoneSurfi].faceZoneNames();
            forAll(fZoneNames, i)
            {
                const faceZone& fZone = fZones[fZoneNames[i]];
                isZonedFace.set(fZone);
            }
        }
    }


    // Determine per pp face which patch it should be in
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Patch that face should be in
    labelList closestPatch(pp.size(), -1);
    {
        // face snap distance as max of point snap distance
        scalarField faceSnapDist(pp.size(), -GREAT);
        {
            // Distance to attract to nearest feature on surface
            const scalarField snapDist
            (
                calcSnapDistance
                (
                    mesh,
                    snapParams,
                    pp
                )
            );

            const faceList& localFaces = pp.localFaces();

            forAll(localFaces, facei)
            {
                const face& f = localFaces[facei];

                forAll(f, fp)
                {
                    faceSnapDist[facei] = max
                    (
                        faceSnapDist[facei],
                        snapDist[f[fp]]
                    );
                }
            }
        }

        pointField localFaceCentres(mesh.faceCentres(), pp.addressing());

        // Get nearest surface and region
        labelList hitSurface;
        labelList hitRegion;
        surfaces.findNearestRegion
        (
            unzonedSurfaces,
            localFaceCentres,
            sqr(faceSnapDist),    // sqr of attract distance
            hitSurface,
            hitRegion
        );

        // Get patch
        forAll(pp, i)
        {
            label facei = pp.addressing()[i];

            if (hitSurface[i] != -1 && !isZonedFace.test(facei))
            {
                closestPatch[i] = globalToMasterPatch_
                [
                    surfaces.globalRegion
                    (
                        hitSurface[i],
                        hitRegion[i]
                    )
                ];
            }
        }
    }


    // Change those faces for which there is a different closest patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList ownPatch(mesh.nFaces(), -1);
    labelList neiPatch(mesh.nFaces(), -1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, i)
        {
            ownPatch[pp.start()+i] = patchi;
            neiPatch[pp.start()+i] = patchi;
        }
    }

    label nChanged = 0;
    forAll(closestPatch, i)
    {
        label facei = pp.addressing()[i];

        if (closestPatch[i] != -1 && closestPatch[i] != ownPatch[facei])
        {
            ownPatch[facei] = closestPatch[i];
            neiPatch[facei] = closestPatch[i];
            nChanged++;
        }
    }

    Info<< "Repatched " << returnReduce(nChanged, sumOp<label>())
        << " faces in = " << mesh.time().cpuTimeIncrement() << " s\n" << nl
        << endl;

    return meshRefiner_.createBaffles(ownPatch, neiPatch);
}


void Foam::snappySnapDriver::detectWarpedFaces
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const faceList& localFaces = pp.localFaces();
    const pointField& localPoints = pp.localPoints();
    const labelList& bFaces = pp.addressing();

    splitFaces.clear();
    splitFaces.setCapacity(bFaces.size());
    splits.clear();
    splits.setCapacity(bFaces.size());

    // Determine parallel consistent normals on points
    const vectorField pointNormals(PatchTools::pointNormals(mesh, pp));

    face f0(4);
    face f1(4);

    forAll(localFaces, facei)
    {
        const face& f = localFaces[facei];

        if (f.size() >= 4)
        {
            // See if splitting face across diagonal would make two faces with
            // biggish normal angle

            labelPair minDiag(-1, -1);
            scalar minCos(GREAT);

            for (label startFp = 0; startFp < f.size()-2; startFp++)
            {
                label minFp = f.rcIndex(startFp);

                for
                (
                    label endFp = f.fcIndex(f.fcIndex(startFp));
                    endFp < f.size() && endFp != minFp;
                    endFp++
                )
                {
                    // Form two faces
                    f0.setSize(endFp-startFp+1);
                    label i0 = 0;
                    for (label fp = startFp; fp <= endFp; fp++)
                    {
                        f0[i0++] = f[fp];
                    }
                    f1.setSize(f.size()+2-f0.size());
                    label i1 = 0;
                    for (label fp = endFp; fp != startFp; fp = f.fcIndex(fp))
                    {
                        f1[i1++] = f[fp];
                    }
                    f1[i1++] = f[startFp];

                    //Info<< "Splitting face:" << f << " into f0:" << f0
                    //    << " f1:" << f1 << endl;

                    const vector n0 = f0.areaNormal(localPoints);
                    const scalar n0Mag = mag(n0);

                    const vector n1 = f1.areaNormal(localPoints);
                    const scalar n1Mag = mag(n1);

                    if (n0Mag > ROOTVSMALL && n1Mag > ROOTVSMALL)
                    {
                        scalar cosAngle = (n0/n0Mag) & (n1/n1Mag);
                        if (cosAngle < minCos)
                        {
                            minCos = cosAngle;
                            minDiag = labelPair(startFp, endFp);
                        }
                    }
                }
            }


            if (minCos < featureCos)
            {
                splitFaces.append(bFaces[facei]);
                splits.append(minDiag);
            }
        }
    }
}


Foam::labelList Foam::snappySnapDriver::getInternalOrBaffleDuplicateFace() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner_.getZones(fzTypes);
    }

    List<labelPair> baffles
    (
        meshRefiner_.subsetBaffles
        (
            mesh,
            internalOrBaffleFaceZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        )
    );

    labelList faceToDuplicate(mesh.nFaces(), -1);
    forAll(baffles, i)
    {
        const labelPair& p = baffles[i];
        faceToDuplicate[p[0]] = p[1];
        faceToDuplicate[p[1]] = p[0];
    }

    return faceToDuplicate;
}


void Foam::snappySnapDriver::doSnap
(
    const dictionary& snapDict,
    const dictionary& motionDict,
    const meshRefinement::FaceMergeType mergeType,
    const scalar featureCos,
    const scalar planarAngle,
    const snapParameters& snapParams
)
{
    if (meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER2)
    {
        // Buffer-layer replacement for this routine

        doSnapBufferLayers
        (
            snapDict,
            motionDict,
            mergeType,
            featureCos,
            planarAngle,
            snapParams
        );
        return;
    }

    addProfiling(snap, "snappyHexMesh::snap");
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;

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



    // Get labels of patches where optional buffer layers are added
    DynamicList<label> bufPatchIDs;
    if (meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER)
    {
        bufPatchIDs.setCapacity(globalToMasterPatch_.size());

        const auto& addLayers =
            meshRefiner_.surfaces().addBufferLayers();

        // Normal patches
        forAll(globalToMasterPatch_, globalRegioni)
        {
            if (addLayers[globalRegioni])
            {
                const label masterP =
                    globalToMasterPatch_[globalRegioni];
                const label slaveP =
                    globalToSlavePatch_[globalRegioni];

                bufPatchIDs.append(masterP);
                if (slaveP != masterP)
                {
                    bufPatchIDs.append(slaveP);
                }
            }
        }

        // Temporary patches from faceZones
        for (const auto& fz : mesh.faceZones())
        {
            label mpI, spI;
            surfaceZonesInfo::faceZoneType type;
            if (meshRefiner_.getFaceZoneInfo(fz.name(), mpI, spI, type))
            {
                bufPatchIDs.appendUniq(mpI);
                bufPatchIDs.appendUniq(spI);
            }
        }
    }



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
    List<labelPair> baffles;
    {
        labelList originatingFaceZone;
        meshRefiner_.createZoneBaffles
        (
            identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone
        );
    }

    // Duplicate points on faceZones of type boundary. Renumber baffles
    // (probably not necessary - faceIDs should not change)
    {
        autoPtr<mapPolyMesh> map = meshRefiner_.dupNonManifoldBoundaryPoints();
        if (map)
        {
            const labelList& reverseFaceMap = map->reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }
    }


    bool doFeatures = false;
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


    bool meshOk = false;


    // Get the labels of added patches.
    const labelList adaptPatchIDs(meshRefiner_.meshedPatches());



    {
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                adaptPatchIDs
            )
        );

        // Distance to attract to nearest feature on surface
        scalarField snapDist(calcSnapDistance(mesh, snapParams, ppPtr()));


        // Construct iterative mesh mover.
        Info<< "Constructing mesh displacer ..." << endl;
        Info<< "Using mesh parameters " << motionDict << nl << endl;

        autoPtr<motionSmoother> meshMoverPtr
        (
            new motionSmoother
            (
                mesh,
                ppPtr(),
                adaptPatchIDs,
                meshRefinement::makeDisplacementField
                (
                    pointMesh::New(mesh),
                    adaptPatchIDs
                ),
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



        // Pre-smooth patch vertices (so before determining nearest)
        preSmoothPatch
        (
            meshRefiner_,
            snapParams,
            nInitErrors,
            internalBaffles,
            meshMoverPtr()
        );

        // Reset moving flag in case we do any topo changes
        mesh.moving(false);


        // Optionally add buffer layers
        if (meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER)
        {
            ////- merge zone baffles (since current buffer layer insertion
            ////- does not handle non-manifold edges ... TBD
            //autoPtr<mapPolyMesh> mapPtr = meshRefiner_.mergeZoneBaffles
            //(
            //    true,   // internal zones
            //    false   // baffle zones
            //);
            //
            //if (mapPtr)
            //{
            //    if (debug & meshRefinement::MESH)
            //    {
            //        const_cast<Time&>(mesh.time())++;
            //        Info<< "Writing baffle-merged mesh to time "
            //            << meshRefiner_.timeName() << endl;
            //        meshRefiner_.write
            //        (
            //            meshRefinement::debugType(debug),
            //            meshRefinement::writeType
            //            (
            //                meshRefinement::writeLevel()
            //              | meshRefinement::WRITEMESH
            //            ),
            //            meshRefiner_.timeName()
            //        );
            //    }
            //}

            ////- use bufferLayer insertion on internalFaceZones
            //Info<< "Adding buffer layers ..." << endl;
            //
            //{
            //    // Remove references to pp
            //    meshMoverPtr.clear();
            //
            //    const labelList meshFaces
            //    (
            //        mesh.faceZones().selection
            //        (
            //            internalFaceZones
            //        ).sortedToc()
            //    );
            //    ppPtr.reset
            //    (
            //        new indirectPrimitivePatch
            //        (
            //            IndirectList<face>
            //            (
            //                mesh.faces(),
            //                meshFaces
            //            ),
            //            mesh.points()
            //        )
            //    );
            //    const pointField thickness
            //    (
            //        wantedThickness(ppPtr(), 1e-1)     //cellSizeFraction
            //      * PatchTools::pointNormals(mesh, ppPtr())
            //    );
            //
            //    // Layer mesh modifier
            //    // - use intrusion, not extrusion
            //    addPatchCellLayer addLayer(mesh, true, false);
            //
            //    // Do mesh changes : introduce point, faces, cells
            //    autoPtr<mapPolyMesh> mapPtr = addBufferLayers
            //    (
            //        ppPtr(),
            //        -thickness,  //1e-3,   //cellSizeFraction,
            //        addLayer
            //    );
            //
            //    // Update numbering on baffles, pointToMaster. Note: uses
            //    // geometric tolerance - could be avoided since now we have
            //    // addPatchCellLayer still intact. However would like to avoid
            //    // use of addPatchCellLayer altogether.
            //    if (mapPtr)
            //    {
            //        // Invalidate extrusion (face numbering might have changed)
            //        ppPtr.clear();
            //
            //        labelList dummyPointToMaster;
            //        snappyLayerDriver::mapFaceZonePoints
            //        (
            //            meshRefiner_,
            //            mapPtr(),
            //            internalBaffles,
            //            dummyPointToMaster
            //        );
            //
            //        if (debug & meshRefinement::MESH)
            //        {
            //            const_cast<Time&>(mesh.time())++;
            //            Info<< "Writing INTERNAL ZONE buffer layer mesh"
            //                << " to time " << meshRefiner_.timeName() << endl;
            //            meshRefiner_.write
            //            (
            //                meshRefinement::debugType(debug),
            //                meshRefinement::writeType
            //                (
            //                    meshRefinement::writeLevel()
            //                  | meshRefinement::WRITEMESH
            //                ),
            //                meshRefiner_.timeName()
            //            );
            //        }
            //    }
            //}

            ////- re-split zone baffles
            //{
            //    labelList originatingFaceZone;
            //    meshRefiner_.createZoneBaffles
            //    (
            //        internalFaceZones,
            //        baffles,
            //        originatingFaceZone
            //    );
            //}


            Info<< "Adding buffer layers ..." << endl;

            // Remove references to pp
            meshMoverPtr.clear();
            // Invalidate extrusion (face numbering might have changed)
            ppPtr.clear();

            // Note: all the way up to addBufferLayers can probably be replaced
            // once addPatchCellLayer can add to both sides of faceZone ...

            // Duplicate points on faceZones that layers are added to
            labelList pointToMaster;
            {
                labelList numLayers(mesh.boundaryMesh().size(), 0);
                UIndirectList<label>(numLayers, bufPatchIDs) = 1;

                autoPtr<mapPolyMesh> mapPtr =
                    snappyLayerDriver::dupFaceZonePoints
                    (
                        meshRefiner_,
                        bufPatchIDs,  // patch indices
                        numLayers,      // num layers per patch
                        baffles,
                        pointToMaster
                    );
                if (mapPtr)
                {
                    // Update numbering on any baffles
                    meshRefinement::mapBaffles
                    (
                        mesh,
                        mapPtr().faceMap(),
                        internalBaffles
                    );


                    //// Shrink back mesh a bit to see if there are any
                    //// incorrect dupFaceZonePoints ...
                    //if (debug & meshRefinement::MESH)
                    //{
                    //    pointField newPoints(mesh.nPoints());
                    //    forAll(newPoints, pointi)
                    //    {
                    //        const auto& pCells = mesh.pointCells()[pointi];
                    //
                    //        point avg(Zero);
                    //        for (const label celli : pCells)
                    //        {
                    //            avg += mesh.cellCentres()[celli];
                    //        }
                    //        avg /= pCells.size();
                    //
                    //        newPoints[pointi] =
                    //            0.5*mesh.points()[pointi]
                    //          + 0.5*avg;
                    //    }
                    //    mesh.movePoints(newPoints);
                    //
                    //    const_cast<Time&>(mesh.time())++;
                    //    Info<< "Writing DEBUG-shrunk layer mesh to time "
                    //        << meshRefiner_.timeName() << endl;
                    //    meshRefiner_.write
                    //    (
                    //        meshRefinement::debugType(debug),
                    //        meshRefinement::writeType
                    //        (
                    //            meshRefinement::writeLevel()
                    //          | meshRefinement::WRITEMESH
                    //        ),
                    //        meshRefiner_.timeName()
                    //    );
                    //}
                }
            }

            //- Not needed: shrinking of mesh since now using intrusion ...
            // // Move mesh back with thickness. Two purposes:
            // // - avoid mapFaceZonePoints below merging points extraneously
            // //   (does not use addPatchCellLayer structure; uses geometric
            // //    tolerance)
            // // - see what is happening
            // {
            //     pointField newPoints(mesh.points());
            //     const auto& mp = ppPtr().meshPoints();
            //     forAll(mp, i)
            //     {
            //         newPoints[mp[i]] -= thickness[i];
            //     }
            //     mesh.movePoints(newPoints);
            //     ppPtr().movePoints(mesh.points());
            //
            //     if (debug & meshRefinement::MESH)
            //     {
            //         const_cast<Time&>(mesh.time())++;
            //         Info<< "Writing shrunk buffer layer mesh to time "
            //             << meshRefiner_.timeName() << endl;
            //         meshRefiner_.write
            //         (
            //             meshRefinement::debugType(debug),
            //             meshRefinement::writeType
            //             (
            //                 meshRefinement::writeLevel()
            //               | meshRefinement::WRITEMESH
            //             ),
            //             meshRefiner_.timeName()
            //         );
            //     }
            // }


            {
                // Layer mesh modifier
                // - use intrusion, not extrusion
                addPatchCellLayer addLayer(mesh, true, false);

                // Redo pp
                autoPtr<indirectPrimitivePatch> bufPatchPtr
                (
                    meshRefinement::makePatch(mesh, bufPatchIDs)
                );

                pointField thickness
                (
                    wantedThickness(bufPatchPtr(), 1e-1)  //cellSizeFraction
                  * PatchTools::pointNormals(mesh, bufPatchPtr())
                );

                // Make sure to adhere to constraints
                {

                    const pointMesh& pMesh = pointMesh::New(mesh);
                    const labelList& mp = bufPatchPtr().meshPoints();

                    tmp<pointVectorField> tdisp
                    (
                        meshRefinement::makeDisplacementField
                        (
                            pMesh,
                            adaptPatchIDs
                        )
                    );
                    // Set internal field
                    UIndirectList<point>(tdisp.ref(), mp) = thickness;
                    // Take over onto boundary field. Needed since constraint
                    // patch might be before the fixedValue patches so its
                    // value gets overwritten
                    for (auto& ppf : tdisp.ref().boundaryFieldRef())
                    {
                        ppf == ppf.patchInternalField();
                    }

                    // Adhere to multi-point constraints
                    const pointConstraints& pcs = pointConstraints::New(pMesh);
                    pcs.constrainDisplacement(tdisp.ref(), false);

                    thickness = UIndirectList<point>(tdisp(), mp);
                }



                // Print a bit
                {
                    // See snappyLayerDriver::calculateLayerThickness.
                    const auto& pbm = mesh.boundaryMesh();
                    label maxLen = 0;
                    for (const label patchi : bufPatchIDs)
                    {
                        maxLen = max(maxLen, label(pbm[patchi].name().size()));
                    }

                    const int oldPrecision = Info.stream().precision();

                    Info<< nl
                        << setf(ios_base::left) << setw(maxLen) << "patch"
                        << setw(0) << " faces    layers thickness[m]" << nl
                        << setf(ios_base::left) << setw(maxLen) << "-----"
                        << setw(0) << " -----    ------ ------------" << endl;

                    for (const label patchi : bufPatchIDs)
                    {
                        Info<< setf(ios_base::left) << setw(maxLen)
                            << pbm[patchi].name() << setprecision(3)
                            << " " << setw(8)
                            << returnReduce(pbm[patchi].size(), sumOp<scalar>())
                            << " " << setw(6) << 1
                            << " " << setw(8) << gAverage(mag(thickness))
                            << endl;
                    }
                    Info<< setprecision(oldPrecision) << endl;
                }


                // Do mesh changes : introduce point, faces, cells
                autoPtr<mapPolyMesh> mapPtr = addBufferLayers
                (
                    bufPatchPtr(),
                    -thickness,  //1e-3,   //cellSizeFraction,
                    addLayer
                );

                // Update numbering on baffles, pointToMaster. Note: uses
                // geometric tolerance - could be avoided since now we have
                // addPatchCellLayer still intact. However would like to avoid
                // use of addPatchCellLayer altogether.
                if (mapPtr)
                {
                    // Invalidate extrusion (face numbering might have changed)
                    ppPtr.clear();

                    snappyLayerDriver::mapFaceZonePoints
                    (
                        meshRefiner_,
                        mapPtr(),
                        internalBaffles,
                        pointToMaster
                    );
                }
            }



            {
                // Merge duplicated points (this creates illegal cells -
                // hopefully they will be smoothed out)
                autoPtr<mapPolyMesh> mapPtr =
                    meshRefiner_.mergePoints(pointToMaster);
                if (mapPtr)
                {
                    // Invalidate extrusion (point numbering might have changed)
                    ppPtr.clear();

                    // Update numbering on any baffles
                    meshRefinement::mapBaffles
                    (
                        mesh,
                        mapPtr().faceMap(),
                        internalBaffles
                    );
                    // Extract baffles across internal faceZones
                    internalBaffles = meshRefinement::subsetBaffles
                    (
                        mesh,
                        internalFaceZones,
                        internalBaffles
                    );

                    if (debug & meshRefinement::MESH)
                    {
                        const_cast<Time&>(mesh.time())++;
                        Info<< "Writing merged points buffer layer mesh"
                            << " to time " << meshRefiner_.timeName() << endl;
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


            Info<< "Inflating buffer layers ..." << endl;

            const pointMesh& pMesh = pointMesh::New(mesh);

            {
                autoPtr<displacementMotionSolver> motionPtr
                (
                    makeMotionSolver
                    (
                        pMesh,
                        snapDict,
                        bufPatchIDs
                        //patchConstraints
                    )
                );

                // Solve internal displacement
                tmp<pointField> tnewPoints(motionPtr->newPoints());

                // Move points
                mesh.movePoints(tnewPoints);

                // Reset moving flag to avoid problems with topo changes
                mesh.moving(false);

                if (debug & meshRefinement::MESH)
                {
                    const_cast<Time&>(mesh.time())++;
                    Info<< "Writing smoothed buffer layer mesh to time "
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

            // Update mesh mover
            ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
            meshMoverPtr.reset
            (
                new motionSmoother
                (
                    mesh,
                    ppPtr(),
                    adaptPatchIDs,
                    meshRefinement::makeDisplacementField
                    (
                        pointMesh::New(mesh),
                        adaptPatchIDs
                    ),
                    motionDict,
                    dryRun_
                )
            );
        }


        //- Only if in feature attraction mode:
        // Nearest feature
        vectorField patchAttraction;
        // Constraints at feature
        List<pointConstraint> patchConstraints;


        //- Any faces to split
        DynamicList<label> splitFaces;
        //- Indices in face to split across
        DynamicList<labelPair> splits;
        //- Patch for both sides of the face
        DynamicList<labelPair> splitPatches;


        for (label iter = 0; iter < nFeatIter; iter++)
        {
            Info<< nl
                << "Morph iteration " << iter << nl
                << "-----------------" << endl;

            // Splitting iteration?
            bool doSplit = false;
            if
            (
                doFeatures
             && snapParams.nFaceSplitInterval() > 0
             && (
                    (iter == nFeatIter-1)
                 || (iter > 0 && (iter%snapParams.nFaceSplitInterval()) == 0)
                )
            )
            {
                doSplit = true;
            }



            indirectPrimitivePatch& pp = ppPtr();
            motionSmoother& meshMover = meshMoverPtr();


            // Calculate displacement at every patch point if we need it:
            // - if automatic near-surface detection
            // - if face splitting active
            pointField nearestPoint;
            vectorField nearestNormal;

            if (snapParams.detectNearSurfacesSnap() || doSplit)
            {
                nearestPoint.setSize(pp.nPoints(), vector::max);
                nearestNormal.setSize(pp.nPoints(), Zero);
            }

            vectorField disp = calcNearestSurface
            (
                snapParams.strictRegionSnap(),  // attract points to region only
                meshRefiner_,
                globalToMasterPatch_,           // for if strictRegionSnap
                globalToSlavePatch_,            // for if strictRegionSnap
                pp,
                pp.localPoints(),
                snapDist,

                nearestPoint,
                nearestNormal
            );


            // Override displacement at thin gaps
            if (snapParams.detectNearSurfacesSnap())
            {
                detectNearSurfaces
                (
                    Foam::cos(degToRad(planarAngle)),// planar cos for gaps
                    pp,
                    pp.localPoints(),
                    nearestPoint,   // surfacepoint from nearest test
                    nearestNormal,  // surfacenormal from nearest test

                    disp
                );
            }

            // Override displacement with feature edge attempt
            if (doFeatures)
            {
                splitFaces.clear();
                splits.clear();
                splitPatches.clear();
                disp = calcNearestSurfaceFeature
                (
                    snapParams,
                    !doSplit,       // alignMeshEdges
                    false,          // no special handling for >=3 patch points
                    iter,
                    featureCos,
                    scalar(iter+1)/nFeatIter,

                    snapDist,
                    disp,
                    nearestNormal,
                    pp,
                    pp.localPoints(),

                    patchAttraction,
                    patchConstraints,

                    splitFaces,
                    splits,
                    splitPatches
                );
            }

            // Check for displacement being outwards.
            outwardsDisplacement(pp, disp);

            // Freeze points on exposed points/faces
            freezeExposedPoints
            (
                meshRefiner_,
                "frozenFaces",      // faceZone name
                "frozenPoints",     // pointZone name
                pp,
                disp
            );

            // Set initial distribution of displacement field (on patches)
            // from patchDisp and make displacement consistent with b.c.
            // on displacement pointVectorField.
            meshMover.setDisplacement(disp);


            if (debug&meshRefinement::ATTRACTION)
            {
                dumpMove
                (
                    mesh.time().path()
                  / "patchDisplacement_" + name(iter) + ".obj",
                    pp.localPoints(),
                    pp.localPoints() + disp
                );
            }

            // Get smoothly varying internal displacement field.
            smoothDisplacement(snapParams, meshMover);

            // Apply internal displacement to mesh.
            meshOk = scaleMesh
            (
                snapParams,
                nInitErrors,
                internalBaffles,
                meshMover
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

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing scaled mesh to time "
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
                Info<< "Writing displacement field ..." << endl;
                meshMover.displacement().write();
                tmp<pointScalarField> magDisp
                (
                    mag(meshMover.displacement())
                );
                magDisp().write();
            }

            // Use current mesh as base mesh
            meshMover.correct();



            // See if any faces need splitting
            label nTotalSplit = returnReduce(splitFaces.size(), sumOp<label>());
            if (nTotalSplit && doSplit)
            {
                // Filter out baffle faces from faceZones of type
                // internal/baffle

                labelList duplicateFace(getInternalOrBaffleDuplicateFace());

                {
                    labelList oldSplitFaces(std::move(splitFaces));
                    List<labelPair> oldSplits(std::move(splits));
                    List<labelPair> oldSplitPatches(std::move(splitPatches));
                    forAll(oldSplitFaces, i)
                    {
                        if (duplicateFace[oldSplitFaces[i]] == -1)
                        {
                            splitFaces.append(oldSplitFaces[i]);
                            splits.append(oldSplits[i]);
                            splitPatches.append(oldSplitPatches[i]);
                        }
                    }
                    nTotalSplit = returnReduce
                    (
                        splitFaces.size(),
                        sumOp<label>()
                    );
                }

                // Update mesh

                // Reset moving flag to avoid meshPhi problems with topo changes
                mesh.moving(false);

                meshRefiner_.splitFacesUndo
                (
                    splitFaces,
                    splits,
                    splitPatches,
                    motionDict,

                    duplicateFace,
                    internalBaffles
                );

                // Redo meshMover
                meshMoverPtr.clear();
                ppPtr.clear();

                // Update mesh mover
                ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
                meshMoverPtr.reset
                (
                    new motionSmoother
                    (
                        mesh,
                        ppPtr(),
                        adaptPatchIDs,
                        meshRefinement::makeDisplacementField
                        (
                            pointMesh::New(mesh),
                            adaptPatchIDs
                        ),
                        motionDict,
                        dryRun_
                    )
                );

                // Update snapping distance
                snapDist = calcSnapDistance(mesh, snapParams, ppPtr());


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
            }


            if (debug&meshRefinement::MESH)
            {
                forAll(internalBaffles, i)
                {
                    const labelPair& p = internalBaffles[i];
                    const point& fc0 = mesh.faceCentres()[p[0]];
                    const point& fc1 = mesh.faceCentres()[p[1]];

                    if (mag(fc0-fc1) > meshRefiner_.mergeDistance())
                    {
                        FatalErrorInFunction
                            << "Separated baffles : f0:" << p[0]
                            << " centre:" << fc0
                            << " f1:" << p[1] << " centre:" << fc1
                            << " distance:" << mag(fc0-fc1)
                            << exit(FatalError);
                    }
                }
            }
        }
    }


    // Reset moving flag to avoid any meshPhi mapping problems
    mesh.moving(false);


    // Merge any introduced baffles (from faceZones of faceType 'internal')
    {
        autoPtr<mapPolyMesh> mapPtr = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false   // baffle zones
        );

        if (mapPtr)
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

    if
    (
        mergeType == meshRefinement::FaceMergeType::GEOMETRIC
     || mergeType == meshRefinement::FaceMergeType::IGNOREPATCH
    )
    {
        labelList duplicateFace(getInternalOrBaffleDuplicateFace());

        // Repatching might have caused faces to be on same patch and hence
        // mergeable so try again to merge coplanar faces. Do not merge baffle
        // faces to ensure they both stay the same.
        label nChanged = meshRefiner_.mergePatchFacesUndo
        (
            featureCos,     // minCos
            featureCos,     // concaveCos
            meshRefiner_.meshedPatches(),
            motionDict,
            duplicateFace,  // faces not to merge
            mergeType
        );

        nChanged += meshRefiner_.mergeEdgesUndo(featureCos, motionDict);

        if (nChanged > 0 && debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing patchFace merged mesh to time "
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

    if (debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
    }
}


// ************************************************************************* //
