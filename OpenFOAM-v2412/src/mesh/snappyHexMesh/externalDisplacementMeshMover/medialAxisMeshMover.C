/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "medialAxisMeshMover.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"
#include "PointEdgeWave.H"
#include "meshRefinement.H"
#include "unitConversion.H"
#include "PatchTools.H"
#include "OBJstream.H"
#include "PointData.H"
#include "zeroFixedValuePointPatchFields.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(medialAxisMeshMover, 0);

    addToRunTimeSelectionTable
    (
        externalDisplacementMeshMover,
        medialAxisMeshMover,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::medialAxisMeshMover::isMaxEdge
(
    const List<PointData<vector>>& pointWallDist,
    const label edgeI,
    const scalar minCos,
    const bool disableWallEdges
) const
{
    const pointField& points = mesh().points();
    const edge& e = mesh().edges()[edgeI];

    if (disableWallEdges)
    {
        // 1. Do not mark edges with one side on moving wall.
        vector v0(points[e[0]] - pointWallDist[e[0]].origin());
        scalar magV0(mag(v0));
        if (magV0 < SMALL)
        {
            return false;
        }

        vector v1(points[e[1]] - pointWallDist[e[1]].origin());
        scalar magV1(mag(v1));
        if (magV1 < SMALL)
        {
            return false;
        }
    }

    //// 2. Do not mark edges with both sides on a moving wall.
    //vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    //scalar magV0(mag(v0));
    //vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    //scalar magV1(mag(v1));
    //if (magV0 < SMALL && magV1 < SMALL)
    //{
    //    return false;
    //}

    //// 3. Detect based on vector to nearest point differing for both endpoints
    //v0 /= magV0;
    //v1 /= magV1;
    //
    //// Test angle.
    //if ((v0 & v1) < minCos)
    //{
    //    return true;
    //}
    //else
    //{
    //    return false;
    //}

    //- Detect based on extrusion vector differing for both endpoints
    //  the idea is that e.g. a sawtooth wall can still be extruded
    //  successfully as long as it is done all to the same direction.
    if ((pointWallDist[e[0]].data() & pointWallDist[e[1]].data()) < minCos)
    {
        return true;
    }

    return false;
}


void Foam::medialAxisMeshMover::update(const dictionary& coeffDict)
{
    Info<< typeName
        << " : Calculating distance to Medial Axis ..." << endl;

    const pointField& points = mesh().points();

    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();


    // Read a few parameters
    // ~~~~~~~~~~~~~~~~~~~~~

    //- Smooth surface normals
    const label nSmoothSurfaceNormals
    (
        meshRefinement::get<label>
        (
            coeffDict,
            "nSmoothSurfaceNormals",
            dryRun_
        )
    );

    const scalar minMedialAxisAngle
    (
        meshRefinement::get<scalar>
        (
            coeffDict,
            "minMedialAxisAngle",
            dryRun_
        )
    );

    const scalar minMedialAxisAngleCos(Foam::cos(degToRad(minMedialAxisAngle)));

    //- Feature angle when to stop adding layers
    const scalar featureAngle
    (
        meshRefinement::get<scalar>(coeffDict, "featureAngle", dryRun_)
    );

    //- When to slip along wall
    const scalar slipFeatureAngle
    (
        coeffDict.getOrDefault<scalar>("slipFeatureAngle", (0.5*featureAngle))
    );

    //- Smooth internal normals
    const label nSmoothNormals
    (
        meshRefinement::get<label>(coeffDict, "nSmoothNormals", dryRun_)
    );

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.getOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );

    const bool disableWallEdges = coeffDict.getOrDefault<bool>
    (
        "disableWallEdges",
        false
    );



    // Predetermine mesh edges
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Precalulate (mesh) master point/edge (only relevant for shared pts/edges)
    const bitSet isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const bitSet isMeshMasterEdge(syncTools::getMasterEdges(mesh()));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const bitSet isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );
    const bitSet isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(PatchTools::pointNormals(mesh(), pp));

    // Smooth patch normal vectors
    fieldSmoother_.smoothPatchNormals
    (
        nSmoothSurfaceNormals,
        isPatchMasterPoint,
        isPatchMasterEdge,
        pp,
        pointNormals
    );


    // Calculate distance to pp points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Distance to wall
    List<PointData<vector>> pointWallDist(mesh().nPoints());

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;


    // 1. Calculate distance to points where displacement is specified.
    {
        // Seed data.
        List<PointData<vector>> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallInfo[patchPointI] = PointData<vector>
            (
                points[pointI],
                0.0,
                pointNormals[patchPointI]     // surface normals
            );
        }

        // Do all calculations
        List<PointData<vector>> edgeWallDist(mesh().nEdges());
        PointEdgeWave<PointData<vector>> wallDistCalc
        (
            mesh(),
            meshPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            0,   // max iterations
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);

        const label nUnvisit = returnReduce
        (
            wallDistCalc.nUnvisitedPoints(),
            sumOp<label>()
        );

        if (nUnvisit > 0)
        {
            if (nMedialAxisIter > 0)
            {
                Info<< typeName
                    << " : Limited walk to " << nMedialAxisIter
                    << " steps. Not visited " << nUnvisit
                    << " out of " << mesh().globalData().nTotalPoints()
                    << " points" << endl;
            }
            else
            {
                WarningInFunction
                    << "Walking did not visit all points." << nl
                    << "    Did not visit " << nUnvisit
                    << " out of " << mesh().globalData().nTotalPoints()
                    << " points. This is not necessarily a problem" << nl
                    << "    and might be due to faceZones splitting of part"
                    << " of the domain." << nl << endl;
            }
        }
    }


    // 2. Find points with max distance and transport information back to
    //    wall.
    {
        List<pointEdgePoint> pointMedialDist(mesh().nPoints());
        List<pointEdgePoint> edgeMedialDist(mesh().nEdges());

        // Seed point data.
        DynamicList<pointEdgePoint> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        // 1. Medial axis points

        const edgeList& edges = mesh().edges();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            if
            (
                !pointWallDist[e[0]].valid(dummyTrackData)
             || !pointWallDist[e[1]].valid(dummyTrackData)
            )
            {
                // Unvisited point. See above about nUnvisit warning
                forAll(e, ep)
                {
                    label pointI = e[ep];

                    if (!pointMedialDist[pointI].valid(dummyTrackData))
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointEdgePoint
                            (
                                points[pointI],
                                0.0
                            )
                        );
                        pointMedialDist[pointI] = maxInfo.last();
                    }
                }

            }
            else if
            (
                isMaxEdge
                (
                    pointWallDist,
                    edgeI,
                    minMedialAxisAngleCos,
                    disableWallEdges
                )
            )
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                //const point medialAxisPt = e.centre(points);
                vector eVec = e.vec(points);
                scalar eMag = mag(eVec);
                if (eMag > VSMALL)
                {
                    eVec /= eMag;

                    // Calculate distance along edge
                    const point& p0 = points[e[0]];
                    const point& origin0 = pointWallDist[e[0]].origin();
                    const point& p1 = points[e[1]];
                    const point& origin1 = pointWallDist[e[1]].origin();
                    scalar dist0 = (p0-origin0) & eVec;
                    scalar dist1 = (origin1-p1) & eVec;
                    scalar s = 0.5*(dist1+eMag+dist0);

                    point medialAxisPt(vector::max);
                    if (s <= dist0)
                    {
                        // Make sure point is not on wall. Note that this
                        // check used to be inside isMaxEdge.
                        if (magSqr((p0-origin0)) > Foam::sqr(SMALL))
                        {
                            medialAxisPt = p0;
                        }
                    }
                    else if (s >= dist0+eMag)
                    {
                        // Make sure point is not on wall. Note that this
                        // check used to be inside isMaxEdge.
                        if (magSqr((p1-origin1)) > Foam::sqr(SMALL))
                        {
                            medialAxisPt = p1;
                        }
                    }
                    else
                    {
                        medialAxisPt = p0+(s-dist0)*eVec;
                    }

                    if (medialAxisPt != vector::max)
                    {
                        forAll(e, ep)
                        {
                            label pointI = e[ep];

                            if (!pointMedialDist[pointI].valid(dummyTrackData))
                            {
                                maxPoints.append(pointI);
                                maxInfo.append
                                (
                                    pointEdgePoint
                                    (
                                        medialAxisPt,   //points[pointI],
                                        magSqr(points[pointI]-medialAxisPt)//0.0
                                    )
                                );
                                pointMedialDist[pointI] = maxInfo.last();
                            }
                        }
                    }
                }
            }
        }


        // 2. Seed non-adapt patches
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        labelHashSet adaptPatches(adaptPatchIDs_);


        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            const pointPatchVectorField& pvf =
                pointDisplacement().boundaryField()[patchI];

            if
            (
                !pp.coupled()
             && !isA<emptyPolyPatch>(pp)
             && !adaptPatches.found(patchI)
            )
            {
                const labelList& meshPoints = pp.meshPoints();

                // Check the type of the patchField. The types are
                //  - fixedValue (0 or more layers) but the >0 layers have
                //    already been handled in the adaptPatches loop
                //  - constraint (but not coupled) types, e.g. symmetryPlane,
                //    slip.
                if (pvf.fixesValue())
                {
                    // Disable all movement on fixedValue patchFields
                    Info<< typeName
                        << " : Inserting all points on patch " << pp.name()
                        << endl;

                    forAll(meshPoints, i)
                    {
                        label pointI = meshPoints[i];
                        if (!pointMedialDist[pointI].valid(dummyTrackData))
                        {
                            maxPoints.append(pointI);
                            maxInfo.append
                            (
                                pointEdgePoint
                                (
                                    points[pointI],
                                    0.0
                                )
                            );
                            pointMedialDist[pointI] = maxInfo.last();
                        }
                    }
                }
                else
                {
                    // Based on geometry: analyse angle w.r.t. nearest moving
                    // point. In the pointWallDist we transported the
                    // normal as the passive vector. Note that this points
                    // out of the originating wall so inside of the domain
                    // on this patch.
                    Info<< typeName
                        << " : Inserting points on patch " << pp.name()
                        << " if angle to nearest layer patch > "
                        << slipFeatureAngle << " degrees." << endl;

                    scalar slipFeatureAngleCos = Foam::cos
                    (
                        degToRad(slipFeatureAngle)
                    );
                    pointField pointNormals
                    (
                        PatchTools::pointNormals(mesh(), pp)
                    );

                    forAll(meshPoints, i)
                    {
                        label pointI = meshPoints[i];

                        if
                        (
                            pointWallDist[pointI].valid(dummyTrackData)
                        && !pointMedialDist[pointI].valid(dummyTrackData)
                        )
                        {
                            // Check if angle not too large.
                            scalar cosAngle =
                            (
                               -pointWallDist[pointI].data()
                              & pointNormals[i]
                            );
                            if (cosAngle > slipFeatureAngleCos)
                            {
                                // Extrusion direction practically perpendicular
                                // to the patch. Disable movement at the patch.

                                maxPoints.append(pointI);
                                maxInfo.append
                                (
                                    pointEdgePoint
                                    (
                                        points[pointI],
                                        0.0
                                    )
                                );
                                pointMedialDist[pointI] = maxInfo.last();
                            }
                            else
                            {
                                // Extrusion direction makes angle with patch
                                // so allow sliding - don't insert zero points
                            }
                        }
                    }
                }
            }
        }

        maxInfo.shrink();
        maxPoints.shrink();


        if (debug)
        {
            mkDir(mesh().time().timePath());
            OBJstream str(mesh().time().timePath()/"medialSurfacePoints.obj");

            pointSet seedPoints
            (
                mesh(),
                "medialSurfacePoints",
                maxPoints
            );

            Info<< typeName
                << " : Writing estimated medial surface:" << nl << incrIndent
                << indent << "locations : " << str.name() << nl
                << indent << "pointSet  : " << seedPoints.name() << nl
                << decrIndent << endl;

            for (const auto& info : maxInfo)
            {
                str.write(info.origin());
            }
            seedPoints.write();
        }


        // Do all calculations
        PointEdgeWave<pointEdgePoint> medialDistCalc
        (
            mesh(),
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            0,   // max iterations
            dummyTrackData
        );
        medialDistCalc.iterate(2*nMedialAxisIter);


        // Extract medial axis distance as pointScalarField
        forAll(pointMedialDist, pointI)
        {
            if (pointMedialDist[pointI].valid(dummyTrackData))
            {
                medialDist_[pointI] = Foam::sqrt
                (
                    pointMedialDist[pointI].distSqr()
                );
                medialVec_[pointI] = pointMedialDist[pointI].origin();
            }
            else
            {
                // Unvisited. Do as if on medial axis so unmoving
                medialDist_[pointI] = 0.0;
                medialVec_[pointI] = point(1, 0, 0);
            }
        }
    }

    // Extract transported surface normals as pointVectorField
    forAll(dispVec_, i)
    {
        if (!pointWallDist[i].valid(dummyTrackData))
        {
            dispVec_[i] = vector(1, 0, 0);
        }
        else
        {
            dispVec_[i] = pointWallDist[i].data();
        }
    }

    // Smooth normal vectors. Do not change normals on pp.meshPoints
    fieldSmoother_.smoothNormals
    (
        nSmoothNormals,
        isMeshMasterPoint,
        isMeshMasterEdge,
        meshPoints,
        dispVec_
    );

    // Calculate ratio point medial distance to point wall distance
    forAll(medialRatio_, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            medialRatio_[pointI] = 0.0;
        }
        else
        {
            scalar wDist2 = pointWallDist[pointI].distSqr();
            scalar mDist = medialDist_[pointI];

            if (wDist2 < sqr(SMALL) && mDist < SMALL)
            //- Note: maybe less strict:
            //(
            //    wDist2 < sqr(meshRefiner_.mergeDistance())
            // && mDist < meshRefiner_.mergeDistance()
            //)
            {
                medialRatio_[pointI] = 0.0;
            }
            else
            {
                medialRatio_[pointI] = mDist / (Foam::sqrt(wDist2) + mDist);
            }
        }
    }


    if (debug)
    {
        Info<< typeName
            << " : Writing medial axis fields:" << nl << incrIndent
            << indent << "ratio of medial distance to wall distance : "
            << medialRatio_.name() << nl
            << indent << "distance to nearest medial axis           : "
            << medialDist_.name() << nl
            << indent << "nearest medial axis location              : "
            << medialVec_.name() << nl
            << indent << "normal at nearest wall                    : "
            << dispVec_.name() << nl
            << decrIndent << endl;

        dispVec_.write();
        medialRatio_.write();
        medialDist_.write();
        medialVec_.write();
    }
}


bool Foam::medialAxisMeshMover::unmarkExtrusion
(
    const label patchPointI,
    pointField& patchDisp,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointI] == snappyLayerDriver::EXTRUDE)
    {
        extrudeStatus[patchPointI] = snappyLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = Zero;
        return true;
    }
    else if (extrudeStatus[patchPointI] == snappyLayerDriver::EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointI] = snappyLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = Zero;
        return true;
    }

    return false;
}


void Foam::medialAxisMeshMover::syncPatchDisplacement
(
    const scalarField& minThickness,
    pointField& patchDisp,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();

    //label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax           // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if (unmarkExtrusion(i, patchDisp, extrudeStatus))
                {
                    nChanged++;
                }
            }
        }

        //labelList syncPatchNLayers(patchNLayers);
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    minEqOp<label>(),
        //    labelMax            // null value
        //);
        //
        //// Reset if differs
        //// 1. take max
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    maxEqOp<label>(),
        //    labelMin            // null value
        //);
        //
        //// Reset if differs
        //// 2. take min
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}

        //nChangedTotal += nChanged;

        if (!returnReduceOr(nChanged))
        {
            break;
        }
    }

    //Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


// Stop layer growth where mesh wraps around edge with a
// large feature angle
void Foam::medialAxisMeshMover::
handleFeatureAngleLayerTerminations
(
    const scalar minCos,
    const bitSet& isPatchMasterPoint,
    const labelList& meshEdges,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    label& nPointCounter
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    // Mark faces that have all points extruded
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList extrudedFaces(pp.size(), true);

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE)
            {
                extrudedFaces[faceI] = false;
                break;
            }
        }
    }



    //label nOldPointCounter = nPointCounter;

    // Detect situation where two featureedge-neighbouring faces are partly or
    // not extruded and the edge itself is extruded. In this case unmark the
    // edge for extrusion.


    List<List<point>> edgeFaceNormals(pp.nEdges());
    List<List<bool>> edgeFaceExtrude(pp.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();
    const vectorField& faceNormals = pp.faceNormals();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        edgeFaceNormals[edgeI].setSize(eFaces.size());
        edgeFaceExtrude[edgeI].setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label faceI = eFaces[i];
            edgeFaceNormals[edgeI][i] = faceNormals[faceI];
            edgeFaceExtrude[edgeI][i] = extrudedFaces[faceI];
        }
    }

    syncTools::syncEdgeList
    (
        mesh(),
        meshEdges,
        edgeFaceNormals,
        ListOps::appendEqOp<point>(),
        List<point>()               // null value
    );

    syncTools::syncEdgeList
    (
        mesh(),
        meshEdges,
        edgeFaceExtrude,
        ListOps::appendEqOp<bool>(),
        List<bool>()                // null value
    );


    forAll(edgeFaceNormals, edgeI)
    {
        const List<point>& eFaceNormals = edgeFaceNormals[edgeI];
        const List<bool>& eFaceExtrude = edgeFaceExtrude[edgeI];

        if (eFaceNormals.size() == 2)
        {
            const edge& e = pp.edges()[edgeI];
            label v0 = e[0];
            label v1 = e[1];

            if
            (
                extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE
             || extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE
            )
            {
                if (!eFaceExtrude[0] || !eFaceExtrude[1])
                {
                    const vector& n0 = eFaceNormals[0];
                    const vector& n1 = eFaceNormals[1];

                    if ((n0 & n1) < minCos)
                    {
                        if (unmarkExtrusion(v0, patchDisp, extrudeStatus))
                        {
                            if (isPatchMasterPoint[v0])
                            {
                                nPointCounter++;
                            }
                        }
                        if (unmarkExtrusion(v1, patchDisp, extrudeStatus))
                        {
                            if (isPatchMasterPoint[v1])
                            {
                                nPointCounter++;
                            }
                        }
                    }
                }
            }
        }
    }

    //Info<< "Added "
    //    << returnReduce(nPointCounter-nOldPointCounter, sumOp<label>())
    //    << " point not to extrude due to minCos "
    //    << minCos << endl;
}


// Find isolated islands (points, edges and faces and layer terminations)
// in the layer mesh and stop any layer growth at these points.
void Foam::medialAxisMeshMover::findIsolatedRegions
(
    const scalar minCosLayerTermination,
    const bool detectExtrusionIsland,
    const bitSet& isPatchMasterPoint,
    const bitSet& isPatchMasterEdge,
    const labelList& meshEdges,
    const scalarField& minThickness,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();


    Info<< typeName << " : Removing isolated regions ..." << nl
        << indent << "- if partially extruded faces make angle < "
        << Foam::radToDeg(Foam::acos(minCosLayerTermination)) <<  nl;
    if (detectExtrusionIsland)
    {
        Info<< indent << "- if exclusively surrounded by non-extruded points"
            << nl;
    }
    else
    {
        Info<< indent << "- if exclusively surrounded by non-extruded faces"
            << nl;
    }

    // Keep count of number of points unextruded
    label nPointCounter = 0;

    while (true)
    {
        // Stop layer growth where mesh wraps around edge with a
        // large feature angle
        if (minCosLayerTermination > -1)
        {
            handleFeatureAngleLayerTerminations
            (
                minCosLayerTermination,
                isPatchMasterPoint,
                meshEdges,

                extrudeStatus,
                patchDisp,
                nPointCounter
            );

            syncPatchDisplacement(minThickness, patchDisp, extrudeStatus);
        }


        // Detect either:
        // - point where all surrounding points are not extruded
        //   (detectExtrusionIsland)
        // or
        // - point where all the faces surrounding it are not fully
        //   extruded

        boolList keptPoints(pp.nPoints(), false);

        if (detectExtrusionIsland)
        {
            // Do not extrude from point where all neighbouring
            // points are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            labelList islandPoint(pp.size(), -1);
            forAll(pp, faceI)
            {
                const face& f = pp.localFaces()[faceI];

                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] != snappyLayerDriver::NOEXTRUDE)
                    {
                        if (islandPoint[faceI] == -1)
                        {
                            // First point to extrude
                            islandPoint[faceI] = f[fp];
                        }
                        else if (islandPoint[faceI] != -2)
                        {
                            // Second or more point to extrude
                            islandPoint[faceI] = -2;
                        }
                    }
                }
            }

            // islandPoint:
            //  -1 : no point extruded on face
            //  -2 : >= 2 points extruded on face
            //  >=0: label of point extruded

            // Check all surrounding faces that I am the islandPoint
            forAll(pointFaces, patchPointI)
            {
                if (extrudeStatus[patchPointI] != snappyLayerDriver::NOEXTRUDE)
                {
                    const labelList& pFaces = pointFaces[patchPointI];

                    forAll(pFaces, i)
                    {
                        label faceI = pFaces[i];
                        if (islandPoint[faceI] != patchPointI)
                        {
                            keptPoints[patchPointI] = true;
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            // Do not extrude from point where all neighbouring
            // faces are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            boolList extrudedFaces(pp.size(), true);
            forAll(pp.localFaces(), faceI)
            {
                const face& f = pp.localFaces()[faceI];
                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE)
                    {
                        extrudedFaces[faceI] = false;
                        break;
                    }
                }
            }

            const labelListList& pointFaces = pp.pointFaces();

            forAll(keptPoints, patchPointI)
            {
                const labelList& pFaces = pointFaces[patchPointI];

                forAll(pFaces, i)
                {
                    label faceI = pFaces[i];
                    if (extrudedFaces[faceI])
                    {
                        keptPoints[patchPointI] = true;
                        break;
                    }
                }
            }
        }


        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            keptPoints,
            orEqOp<bool>(),
            false               // null value
        );

        label nChanged = 0;

        forAll(keptPoints, patchPointI)
        {
            if (!keptPoints[patchPointI])
            {
                if (unmarkExtrusion(patchPointI, patchDisp, extrudeStatus))
                {
                    nPointCounter++;
                    nChanged++;
                }
            }
        }


        if (!returnReduceOr(nChanged))
        {
            break;
        }
    }

    const edgeList& edges = pp.edges();


    // Count number of mesh edges using a point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList isolatedPoint(pp.nPoints(), Zero);

    forAll(edges, edgeI)
    {
        if (isPatchMasterEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            if (extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v0] += 1;
            }
            if (extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v1] += 1;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        isolatedPoint,
        plusEqOp<label>(),
        label(0)        // null value
    );

    // stop layer growth on isolated faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forAll(pp, faceI)
    {
        const face& f = pp.localFaces()[faceI];
        bool failed = false;
        forAll(f, fp)
        {
            if (isolatedPoint[f[fp]] > 2)
            {
                failed = true;
                break;
            }
        }
        bool allPointsExtruded = true;
        if (!failed)
        {
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE)
                {
                    allPointsExtruded = false;
                    break;
                }
            }

            if (allPointsExtruded)
            {
                forAll(f, fp)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            f[fp],
                            patchDisp,
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;
                    }
                }
            }
        }
    }

    Info<< typeName
        << " : Number of isolated points extrusion stopped : "
        << returnReduce(nPointCounter, sumOp<label>()) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::medialAxisMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement,
    const bool dryRun
)
:
    externalDisplacementMeshMover(dict, baffles, pointDisplacement, dryRun),
    adaptPatchIDs_(getFixedValueBCs(pointDisplacement)),
    adaptPatchPtr_(getPatch(mesh(), adaptPatchIDs_)),
    scale_
    (
        IOobject
        (
            "scale",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh().points()),
    meshMover_
    (
        const_cast<polyMesh&>(mesh()),
        const_cast<pointMesh&>(pMesh()),
        adaptPatchPtr_(),
        pointDisplacement,
        scale_,
        oldPoints_,
        adaptPatchIDs_,
        dict,
        dryRun
    ),
    fieldSmoother_(mesh()),
    dispVec_
    (
        IOobject
        (
            "dispVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pMesh(),
        dimensionedVector(dimLength, Zero)
    ),
    medialRatio_
    (
        IOobject
        (
            "medialRatio",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pMesh(),
        dimensionedScalar(dimless, Zero)
    ),
    medialDist_
    (
        IOobject
        (
            "pointMedialDist",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    medialVec_
    (
        IOobject
        (
            "medialVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pMesh(),
        dimensionedVector(dimLength, Zero)
    )
{
    update(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::~medialAxisMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::medialAxisMeshMover::calculateDisplacement
(
    const dictionary& coeffDict,
    const scalarField& minThickness,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
)
{
    Info<< typeName << " : Smoothing using Medial Axis ..." << endl;

    const indirectPrimitivePatch& pp = *adaptPatchPtr_;
    const labelList& meshPoints = pp.meshPoints();


    // Read settings
    // ~~~~~~~~~~~~~

    //- (lambda-mu) smoothing of internal displacement
    const label nSmoothDisplacement =
        coeffDict.getOrDefault("nSmoothDisplacement", 0);

    //- Layer thickness too big
    const scalar maxThicknessToMedialRatio =
        coeffDict.get<scalar>("maxThicknessToMedialRatio");

    //- Feature angle when to stop adding layers
    const scalar featureAngle = coeffDict.get<scalar>("featureAngle");

    //- Stop layer growth where mesh wraps around sharp edge
    scalar layerTerminationAngle = coeffDict.getOrDefault<scalar>
    (
        "layerTerminationAngle",
        0.5*featureAngle
    );
    scalar minCosLayerTermination = Foam::cos(degToRad(layerTerminationAngle));

    //- Smoothing wanted patch thickness
    const label nSmoothPatchThickness =
        coeffDict.get<label>("nSmoothThickness");

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.getOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );

    //- Use strict extrusionIsland detection
    const bool detectExtrusionIsland = coeffDict.getOrDefault
    (
        "detectExtrusionIsland",
        false
    );


    // Precalulate master points/edge (only relevant for shared points/edges)
    const bitSet isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const bitSet isMeshMasterEdge(syncTools::getMasterEdges(mesh()));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const bitSet isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );
    const bitSet isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );


    scalarField thickness(mag(patchDisp));

    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == snappyLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }

    label numThicknessRatioExclude = 0;

    // reduce thickness where thickness/medial axis distance large
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<OBJstream> str;
    if (debug)
    {
        str.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "thicknessRatioExcludePoints_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing points with too large an extrusion distance to "
            << str().name() << endl;
    }

    autoPtr<OBJstream> medialVecStr;
    if (debug)
    {
        medialVecStr.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "thicknessRatioExcludeMedialVec_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing medial axis vectors on points with too large"
            << " an extrusion distance to " << medialVecStr().name() << endl;
    }

    forAll(meshPoints, patchPointI)
    {
        if (extrudeStatus[patchPointI] != snappyLayerDriver::NOEXTRUDE)
        {
            label pointI = meshPoints[patchPointI];

            //- Option 1: look only at extrusion thickness v.s. distance
            //  to nearest (medial axis or static) point.
            scalar mDist = medialDist_[pointI];
            scalar thicknessRatio = thickness[patchPointI]/(mDist+VSMALL);

            //- Option 2: Look at component in the direction
            //  of nearest (medial axis or static) point.
            const vector n = normalised(patchDisp[patchPointI]);
            const vector mVec =
                normalised
                (
                    medialVec_[pointI] - mesh().points()[pointI]
                );

            thicknessRatio *= (n & mVec);

            if (thicknessRatio > maxThicknessToMedialRatio)
            {
                // Truncate thickness.
                if (debug&2)
                {
                    Pout<< "truncating displacement at "
                        << mesh().points()[pointI]
                        << " from " << thickness[patchPointI]
                        << " to "
                        <<  0.5
                           *(
                                minThickness[patchPointI]
                               +thickness[patchPointI]
                            )
                        << " medial direction:" << mVec
                        << " extrusion direction:" << n
                        << " with thicknessRatio:" << thicknessRatio
                        << endl;
                }

                thickness[patchPointI] =
                    0.5*(minThickness[patchPointI]+thickness[patchPointI]);

                patchDisp[patchPointI] = thickness[patchPointI]*n;

                if (isPatchMasterPoint[patchPointI])
                {
                    numThicknessRatioExclude++;
                }

                if (str)
                {
                    const point& pt = mesh().points()[pointI];
                    str().writeLine(pt, pt+patchDisp[patchPointI]);
                }
                if (medialVecStr)
                {
                    const point& pt = mesh().points()[pointI];
                    medialVecStr().writeLine(pt, medialVec_[pointI]);
                }
            }
        }
    }

    Info<< typeName << " : Reducing layer thickness at "
        << returnReduce(numThicknessRatioExclude, sumOp<label>())
        << " nodes where thickness to medial axis distance is large " << endl;


    // find points where layer growth isolated to a lone point, edge or face

    findIsolatedRegions
    (
        minCosLayerTermination,
        detectExtrusionIsland,

        isPatchMasterPoint,
        isPatchMasterEdge,
        meshEdges,
        minThickness,

        extrudeStatus,
        patchDisp
    );

    // Update thickness for changed extrusion
    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == snappyLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }


    // Smooth layer thickness on moving patch. Since some locations will have
    // disabled the extrusion this might cause big jumps in wanted displacement
    // for neighbouring patch points. So smooth the wanted displacement
    // before actually trying to move the mesh.
    fieldSmoother_.minSmoothField
    (
        nSmoothPatchThickness,
        isPatchMasterPoint,
        isPatchMasterEdge,
        pp,
        minThickness,
        thickness
    );


    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    List<PointData<scalar>> pointWallDist(mesh().nPoints());

    const pointField& points = mesh().points();
    // 1. Calculate distance to points where displacement is specified.
    // This wave is used to transport layer thickness
    {
        // Distance to wall and medial axis on edges.
        List<PointData<scalar>> edgeWallDist(mesh().nEdges());
        labelList wallPoints(meshPoints.size());

        // Seed data.
        List<PointData<scalar>> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallPoints[patchPointI] = pointI;
            wallInfo[patchPointI] = PointData<scalar>
            (
                points[pointI],
                0.0,
                thickness[patchPointI]        // transport layer thickness
            );
        }

        // Do all calculations
        PointEdgeWave<PointData<scalar>> wallDistCalc
        (
            mesh(),
            wallPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            0,   // max iterations
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);
    }


    // Calculate scaled displacement vector
    pointField& displacement = pointDisplacement_;

    forAll(displacement, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            displacement[pointI] = Zero;
        }
        else
        {
            // 1) displacement on nearest wall point, scaled by medialRatio
            //    (wall distance / medial distance)
            // 2) pointWallDist[pointI].data() is layer thickness transported
            //    from closest wall point.
            // 3) shrink in opposite direction of addedPoints
            displacement[pointI] =
                -medialRatio_[pointI]
                *pointWallDist[pointI].data()
                *dispVec_[pointI];
        }
    }


    // Smear displacement away from fixed values (medialRatio=0 or 1)
    if (nSmoothDisplacement > 0)
    {
        bitSet isToBeSmoothed(displacement.size(), false);

        forAll(displacement, i)
        {
            if (medialRatio_[i] > SMALL && medialRatio_[i] < 1-SMALL)
            {
                isToBeSmoothed.set(i);
            }
        }

        fieldSmoother_.smoothLambdaMuDisplacement
        (
            nSmoothDisplacement,
            isMeshMasterPoint,
            isMeshMasterEdge,
            isToBeSmoothed,
            displacement
        );
    }

    if (str)
    {
        Info<< typeName
            << " : Written " << returnReduce(str().nVertices(), sumOp<label>())
            << " points with too large an extrusion distance to "
            << str().name() << endl;
    }
    if (medialVecStr)
    {
        Info<< typeName
            << " : Written "
            << returnReduce(medialVecStr().nVertices(), sumOp<label>())
            << " medial axis vectors on points with too large"
            << " an extrusion distance to " << medialVecStr().name() << endl;
    }
}


bool Foam::medialAxisMeshMover::shrinkMesh
(
    const dictionary& meshQualityDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    //- Number of attempts shrinking the mesh
    const label nSnap = meshQualityDict.get<label>("nRelaxIter");


    // Make sure displacement boundary conditions is uptodate with
    // internal field
    meshMover_.setDisplacementPatchFields();

    Info<< typeName << " : Moving mesh ..." << endl;
    scalar oldErrorReduction = -1;

    bool meshOk = false;

    for (label iter = 0; iter < 2*nSnap ; iter++)
    {
        Info<< typeName
            << " : Iteration " << iter << endl;
        if (iter == nSnap)
        {
            Info<< typeName
                << " : Displacement scaling for error reduction set to 0."
                << endl;
            oldErrorReduction = meshMover_.setErrorReduction(0.0);
        }

        if
        (
            meshMover_.scaleMesh
            (
                checkFaces,
                baffles_,
                meshMover_.paramDict(),
                meshQualityDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< typeName << " : Successfully moved mesh" << endl;
            meshOk = true;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover_.setErrorReduction(oldErrorReduction);
    }

    Info<< typeName << " : Finished moving mesh ..." << endl;

    return meshOk;
}


bool Foam::medialAxisMeshMover::move
(
    const dictionary& moveDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    // Read a few settings
    // ~~~~~~~~~~~~~~~~~~~

    //- Name of field specifying min thickness
    const word minThicknessName = moveDict.get<word>("minThicknessName");


    // Extract out patch-wise displacement
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    scalarField zeroMinThickness;
    if (minThicknessName == "none")
    {
        zeroMinThickness = scalarField(pp.nPoints(), Zero);
    }
    const scalarField& minThickness =
    (
        (minThicknessName == "none")
      ? zeroMinThickness
      : mesh().lookupObject<scalarField>(minThicknessName)
    );


    pointField patchDisp
    (
        pointDisplacement_.internalField(),
        pp.meshPoints()
    );

    List<snappyLayerDriver::extrudeMode> extrudeStatus
    (
        pp.nPoints(),
        snappyLayerDriver::EXTRUDE
    );
    forAll(extrudeStatus, pointI)
    {
        if (mag(patchDisp[pointI]) <= minThickness[pointI]+SMALL)
        {
            extrudeStatus[pointI] = snappyLayerDriver::NOEXTRUDE;
        }
    }


    // Solve displacement
    calculateDisplacement(moveDict, minThickness, extrudeStatus, patchDisp);

    //- Move mesh according to calculated displacement
    return shrinkMesh
    (
        moveDict,           // meshQualityDict,
        nAllowableErrors,   // nAllowableErrors
        checkFaces
    );
}


void Foam::medialAxisMeshMover::movePoints(const pointField& p)
{
    externalDisplacementMeshMover::movePoints(p);

    // Update motionSmoother for new geometry (moves adaptPatchPtr_)
    meshMover_.movePoints();

    // Assume current mesh location is correct (reset oldPoints, scale)
    meshMover_.correct();
}


// ************************************************************************* //
