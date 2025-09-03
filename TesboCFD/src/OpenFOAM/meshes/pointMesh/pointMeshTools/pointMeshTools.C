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

\*---------------------------------------------------------------------------*/

#include "pointMeshTools.H"
#include "syncTools.H"
#include "PatchTools.H"
#include "dummyTransform.H"
#include "unitConversion.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::pointMeshTools::featurePointsEdges
(
    const polyMesh& mesh,

    const uindirectPrimitivePatch& allBoundary,

    // Per boundary face to zone
    const labelUList& faceToZone,
    // Number of zones
    const label nZones,

    const scalar edgeFeatureAngle,
    //const scalar pointFeatureAngle, //not yet done

    // Feature edge(points) internal to a zone
    labelListList& zoneToMeshPoints,
    List<pointConstraintList>& zoneToConstraints,

    // Feature edge(points) in between zones
    labelList& twoZoneMeshPoints,
    pointConstraintList& twoZoneConstraints,

    // Feature points on > 2 zones
    labelList& multiZoneMeshPoints,
    pointConstraintList& multiZoneConstraints
)
{
    // Analyses edges on mesh faces and splits them up according to topology
    // and geometry:
    // - edges in between faces on same zone but making a feature angle
    // - edges in between faces on two different zones
    // - points on faces on > 2 zones


    const globalMeshData& globalData = mesh.globalData();
    const indirectPrimitivePatch& cpp = globalData.coupledPatch();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const auto& mp = allBoundary.meshPoints();

    const vector nullVector(vector::uniform(0));
    const auto assignNonNull = [&](vector& x, const vector& y)
    {
        if (x == nullVector && y != nullVector)
        {
            x = y;
        }
    };

    // Calculate parallel-consistent point normals (as unweighted average
    // of faceNormals). Note: only valid on patch points, not on mesh points
    // that are coupled to these.
    const pointField pointNormals
    (
        PatchTools::pointNormals
        (
            mesh,
            allBoundary
        )
    );


    // Find correspondence between allBoundary and coupled edges
    labelList allEdges;
    labelList coupledEdges;
    bitSet sameEdgeOrientation;
    PatchTools::matchEdges
    (
        allBoundary,
        cpp,
        allEdges,
        coupledEdges,
        sameEdgeOrientation
    );



    // To construct the patches we need to know per edge:
    //  - patch on either side (if topological feature edge)
    //  - faceNormal on either side (if feature angle)
    // We need to know per point:
    //  - patches on all connected faces
    //  - faceNormals on all connected faces? And compare to average?
    //    or edge normals on all connected edges


    typedef Tuple2<label, vector> PN;
    const PN nullPN(-1, nullVector);



    // Point-based analysis
    // ~~~~~~~~~~~~~~~~~~~~

    // Collect per (mesh)point the zones (1, 2 or >2). Note: per mesh to
    // make it easier to sync. See edge-based code below where we explicitly
    // have to transfer from patch-edge to mesh-point etc. Note sure which one
    // fits better....
    // State:
    //  (-1, -1)    : initial state
    //  (-2, -2)    : more than 2 zones
    //  (>=0, >=0)  : zones from connected faces
    labelPairList pointToZones(mesh.nPoints(), labelPair(-1, -1));
    {
        // Combine zones.
        const auto combineZones = [&](labelPair& x, const labelPair& y)
        {
            if (x == labelPair(-2, -2))
            {
                // Already marked as having >2 zones
            }
            else if (y == labelPair(-2, -2))
            {
                x = y;
            }
            else
            {
                // Find first free slot
                if (x[0] == -1)
                {
                    if (y[0] != -1)
                    {
                        x[0] = y[0];
                    }
                    else
                    {
                        x[0] = y[1];
                    }
                }
                else if (x[1] == -1)
                {
                    if (y[0] != -1 && y[0] != x[0])
                    {
                        x[1] = y[0];
                    }
                    else if (y[1] != -1 && y[1] != x[1])
                    {
                        x[1] = y[1];
                    }
                }
                else
                {
                    // Both x slots filled. See if y adds a 3rd element
                    if (y[0] != -1 && y[0] != x[0] && y[0] != x[1])
                    {
                        x = labelPair(-2, -2);
                    }
                    else if (y[1] != -1 && y[1] != x[0] && y[1] != x[1])
                    {
                        x = labelPair(-2, -2);
                    }
                }
            }
        };


        forAll(allBoundary, facei)
        {
            const auto& f = allBoundary[facei];
            const label zonei = faceToZone[facei];
            for (const label pointi : f)
            {
                auto& pZones = pointToZones[pointi];

                if (pZones != labelPair(-2, -2) && !pZones.contains(zonei))
                {
                    if (pZones.first() == -1)
                    {
                        pZones.first() = zonei;
                    }
                    else if (pZones.second() == -1)
                    {
                        pZones.second() = zonei;
                    }
                    else
                    {
                        // Mark as >2 zones
                        pZones = labelPair(-2, -2);
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointToZones,
            combineZones,
            labelPair(-1, -1),
            dummyTransform()
        );
    }




    // Edge-based analysis
    // ~~~~~~~~~~~~~~~~~~~~

    // 1. Local analysis

    List<Pair<PN>> allEdgeToFaces
    (
        allBoundary.nEdges(),
        Pair<PN>(nullPN, nullPN)
    );
    {
        const auto& edgeFaces = allBoundary.edgeFaces();
        const auto& faceNormals = allBoundary.faceNormals();

        forAll(edgeFaces, edgei)
        {
            const auto& eFaces = edgeFaces[edgei];
            const vector& n0 = faceNormals[eFaces[0]];
            const label zone0 = faceToZone[eFaces[0]];
            if (eFaces.size() == 1)
            {
                allEdgeToFaces[edgei] = Pair<PN>(PN(zone0, n0), nullPN);
            }
            else
            {
                const vector& n1 = faceNormals[eFaces[1]];
                const label zone1 = faceToZone[eFaces[1]];
                allEdgeToFaces[edgei] = Pair<PN>
                (
                    PN(zone0, n0),
                    PN(zone1, n1)
                );
            }
        }
    }


    // 2. Sync across coupled patches

    {
        // Combine pair of normals
        const auto vectorPairMax = [&](Pair<PN>& x, const Pair<PN>& y)
        {
            if (x[0] == nullPN)
            {
                if (y[0] != nullPN)
                {
                    x[0] = y[0];
                }
                else
                {
                    x[0] = y[1];
                }
            }
            else if (x[1] == nullPN)
            {
                if (y[0] != nullPN && y[0] != x[0])
                {
                    x[1] = y[0];
                }
                else
                {
                    x[1] = y[1];
                }
            }
        };

        List<Pair<PN>> cppEdgeData
        (
            map.constructSize(),
            Pair<PN>(nullPN, nullPN)
        );
        UIndirectList<Pair<PN>>(cppEdgeData, coupledEdges) =
            UIndirectList<Pair<PN>>(allEdgeToFaces, allEdges);

        globalData.syncData
        (
            cppEdgeData,
            globalData.globalEdgeSlaves(),
            globalData.globalEdgeTransformedSlaves(),
            map,
            globalData.globalTransforms(),
            vectorPairMax,
            dummyTransform()
        );

        UIndirectList<Pair<PN>>(allEdgeToFaces, allEdges) =
            UIndirectList<Pair<PN>>(cppEdgeData, coupledEdges);
    }


    // Now we have all the per-patch edge information
    // - do inter-patch edges
    // - do feature-angle edges
    // Store on mesh points

    const auto assignNonNullPN = [&](PN& x, const PN& y)
    {
        if (x.second() == nullVector && y.second() != nullVector)
        {
            x = y;
        }
    };


    const scalar featEdgeCos = Foam::cos(degToRad(edgeFeatureAngle));


    //OFstream os("allBoundary.obj");
    //Pout<< "Dumping feature edges to " << os.name() << endl;
    //OFstream interOs("interZoneBoundary.obj");
    //Pout<< "Dumping inter-zone edges to " << os.name() << endl;

    // Storing the normal for points that are on inter-patch edges
    vectorField patchEdgeNormal(mesh.nPoints(), nullVector);
    // Storing the constraint for points that are on patch-internal feat edges
    List<PN> featEdgeNormal(mesh.nPoints(), nullPN);

    // Use point-based sync
    {
        forAll(allEdgeToFaces, edgei)
        {
            const edge& e = allBoundary.edges()[edgei];
            const label mp0 = mp[e[0]];
            const label mp1 = mp[e[1]];

            const Pair<PN>& facesInfo = allEdgeToFaces[edgei];

            if (facesInfo[1] == nullPN)
            {
                // Real boundary edge. On single patch only so no need
                // to put in separate point patch ... (? tbd)
            }
            else
            {
                const label zone0 = facesInfo[0].first();
                const label zone1 = facesInfo[1].first();

                const point& p0 = allBoundary.points()[mp0];
                const point& p1 = allBoundary.points()[mp1];
                vector eVec(p1-p0);
                eVec.normalise();

                if (zone0 != zone1)
                {
                    // Inter-patch. TBD: check for feature angle as well?

                    //patchEdgeNormal[mp0] = pointNormals[e[0]];
                    //patchEdgeNormal[mp1] = pointNormals[e[1]];

                    patchEdgeNormal[mp0] = eVec;
                    patchEdgeNormal[mp1] = eVec;
                }
                else
                {
                    // Same patch - check for feature angle

                    const vector& n0 = facesInfo[0].second();
                    const vector& n1 = facesInfo[1].second();

                    if ((n0 & n1) < featEdgeCos)
                    {
                        //Pout<< "** feature edge:" << edgei
                        //    << " points:" << allBoundary.points()[mp0]
                        //    << allBoundary.points()[mp1]
                        //    << nl
                        //    << " faceNormal0:" << n0
                        //    << " faceNormal1:" << n1 << nl
                        //    << " zone0:" << zone0
                        //    << " zone1:" << zone1 << nl
                        //    << " pointNormals0:" << pointNormals[e[0]]
                        //    << " pointNormals1:" << pointNormals[e[1]]
                        //    << nl
                        //    << endl;

                        if (patchEdgeNormal[mp0] == nullVector)
                        {
                            //featEdgeNormal[mp0] = PN
                            //(
                            //    zone0,   // zone
                            //    pointNormals[e[0]]
                            //);
                            featEdgeNormal[mp0].first() = zone0;
                            // Assign edge direction. TBD: average & parallel
                            featEdgeNormal[mp0].second() = eVec;
                        }
                        if (patchEdgeNormal[mp1] == nullVector)
                        {
                            //featEdgeNormal[mp1] = PN
                            //(
                            //    zone1,   // zone
                            //    pointNormals[e[1]]
                            //);
                            featEdgeNormal[mp1].first() = zone1;
                            // Assign edge direction. TBD: average & parallel
                            featEdgeNormal[mp1].second() = eVec;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            patchEdgeNormal,
            assignNonNull,
            nullVector
        );
        syncTools::syncPointList
        (
            mesh,
            featEdgeNormal,
            assignNonNullPN,
            nullPN,
            dummyTransform()
        );
    }

    // Make sure that inter-patch points are not also in feature-edge
    // points. Note: not absolutely necessary since all inter-patch points
    // will also be in the 'normal' facePointPatches.

    DynamicList<label> dynMultiZoneMeshPoints(allBoundary.nPoints());
    DynamicList<pointConstraint> dynMultiZoneConstraints(allBoundary.nPoints());
    forAll(pointToZones, pointi)
    {
        if (pointToZones[pointi] == labelPair(-2, -2))
        {
            dynMultiZoneMeshPoints.append(pointi);
            // Note: pointConstraint just a slip constraint for now
            dynMultiZoneConstraints.append
            (
                pointConstraint
                (
                    3,                  // feature point
                    Zero                //meshPointNormals[pointi]
                )
            );
            // Unmark as feature angle point
            patchEdgeNormal[pointi] = nullVector;
            featEdgeNormal[pointi] = nullPN;
        }
    }
    multiZoneMeshPoints = std::move(dynMultiZoneMeshPoints);
    multiZoneConstraints = std::move(dynMultiZoneConstraints);


    DynamicList<label> dynTwoZoneMeshPoints(allBoundary.nPoints());
    DynamicList<pointConstraint> dynTwoZoneConstraints(allBoundary.nPoints());
    forAll(patchEdgeNormal, pointi)
    {
        if (patchEdgeNormal[pointi] != nullVector)
        {
            dynTwoZoneMeshPoints.append(pointi);
            // Note: pointConstraint just a slip constraint for now
            dynTwoZoneConstraints.append
            (
                pointConstraint
                (
                    2,                  // feature edge
                    patchEdgeNormal[pointi]
                )
            );

            // Unmark as feature angle point
            featEdgeNormal[pointi] = nullPN;
        }
    }
    twoZoneMeshPoints = std::move(dynTwoZoneMeshPoints);
    twoZoneConstraints = std::move(dynTwoZoneConstraints);


    // Sort featEdgeNormal according to zone
    zoneToMeshPoints.resize_nocopy(nZones);
    zoneToConstraints.resize_nocopy(nZones);
    {
        labelList sizes(nZones, 0);
        forAll(featEdgeNormal, pointi)
        {
            const auto& pInfo = featEdgeNormal[pointi];
            if (pInfo != nullPN)
            {
                const label zonei = pInfo.first();
                sizes[zonei]++;
            }
        }
        forAll(zoneToMeshPoints, zonei)
        {
            zoneToMeshPoints[zonei].setSize(sizes[zonei]);
            zoneToConstraints[zonei].setSize(sizes[zonei]);
        }
        sizes = 0;
        forAll(featEdgeNormal, pointi)
        {
            const auto& pInfo = featEdgeNormal[pointi];
            if (pInfo != nullPN)
            {
                const label zonei = pInfo.first();
                const label index = sizes[zonei]++;

                zoneToMeshPoints[zonei][index] = pointi;
                zoneToConstraints[zonei][index] =
                    pointConstraint(2, pInfo.second()); // constrained to slide
            }
        }
    }
}


// ************************************************************************* //
