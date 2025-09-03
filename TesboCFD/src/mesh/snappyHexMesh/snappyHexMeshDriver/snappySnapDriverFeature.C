/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2022,2024 OpenCFD Ltd.
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

#include "snappySnapDriver.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "OBJstream.H"
#include "motionSmoother.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "unitConversion.H"
#include "plane.H"
#include "featureEdgeMesh.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "snapParameters.H"
#include "PatchTools.H"
#include "pyramid.H"
#include "localPointRegion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class T>
    class listPlusEqOp
    {
    public:

        void operator()(List<T>& x, const List<T>& y) const
        {
            label sz = x.size();
            x.setSize(sz+y.size());
            forAll(y, i)
            {
                x[sz++] = y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//bool Foam::snappySnapDriver::isFeaturePoint
//(
//    const scalar featureCos,
//    const indirectPrimitivePatch& pp,
//    const bitSet& isFeatureEdge,
//    const label pointi
//) const
//{
//    const pointField& points = ppLocalPoints;
//    const edgeList& edges = pp.edges();
//    const labelList& pEdges = pp.pointEdges()[pointi];
//
//    label nFeatEdges = 0;
//
//    forAll(pEdges, i)
//    {
//        if (isFeatureEdge[pEdges[i]])
//        {
//            nFeatEdges++;
//
//            for (label j = i+1; j < pEdges.size(); j++)
//            {
//                if (isFeatureEdge[pEdges[j]])
//                {
//                    const edge& ei = edges[pEdges[i]];
//                    const edge& ej = edges[pEdges[j]];
//
//                    const point& p = points[pointi];
//                    const point& pi = points[ei.otherVertex(pointi)];
//                    const point& pj = points[ej.otherVertex(pointi)];
//
//                    vector vi = p-pi;
//                    scalar viMag = mag(vi);
//
//                    vector vj = pj-p;
//                    scalar vjMag = mag(vj);
//
//                    if
//                    (
//                        viMag > SMALL
//                     && vjMag > SMALL
//                     && ((vi/viMag & vj/vjMag) < featureCos)
//                    )
//                    {
//                        return true;
//                    }
//                }
//            }
//        }
//    }
//
//    if (nFeatEdges == 1)
//    {
//        // End of feature-edge string
//        return true;
//    }
//
//    return false;
//}


void Foam::snappySnapDriver::smoothAndConstrain
(
    const bitSet& isPatchMasterEdge,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        // Calculate average displacement of neighbours
        // - unconstrained (i.e. surface) points use average of all
        //   neighbouring points
        // - from testing it has been observed that it is not beneficial
        //   to have edge constrained points use average of all edge or point
        //   constrained neighbours since they're already attracted to
        //   the nearest point on the feature.
        //   Having them attract to point-constrained neighbours does not
        //   make sense either since there is usually just one of them so
        //   it severely distorts it.
        // - same for feature points. They are already attracted to the
        //   nearest feature point.

        vectorField dispSum(pp.nPoints(), Zero);
        labelList dispCount(pp.nPoints(), Zero);

        const labelListList& pointEdges = pp.pointEdges();
        const edgeList& edges = pp.edges();

        forAll(pointEdges, pointi)
        {
            const labelList& pEdges = pointEdges[pointi];

            label nConstraints = constraints[pointi].first();

            if (nConstraints <= 1)
            {
                forAll(pEdges, i)
                {
                    label edgei = pEdges[i];

                    if (isPatchMasterEdge[edgei])
                    {
                        label nbrPointi = edges[edgei].otherVertex(pointi);
                        if (constraints[nbrPointi].first() >= nConstraints)
                        {
                            dispSum[pointi] += disp[nbrPointi];
                            dispCount[pointi]++;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispSum,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispCount,
            plusEqOp<label>(),
            label(0),
            mapDistribute::transform()
        );

        // Constraints
        forAll(constraints, pointi)
        {
            if (dispCount[pointi] > 0)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointi] =
                    0.5
                   *(disp[pointi] + dispSum[pointi]/dispCount[pointi]);
            }
        }
    }
}


void Foam::snappySnapDriver::calcNearestFace
(
    const label iter,
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,

    const scalarField& faceSnapDist,
    vectorField& faceDisp,
    vectorField& faceSurfaceNormal,
    labelList& faceSurfaceGlobalRegion
    //vectorField& faceRotation
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    const pointField ppFaceCentres
    (
        primitivePatch
        (
            SubList<face>(pp.localFaces()),
            ppLocalPoints
        ).faceCentres()
    );


    // Displacement and orientation per pp face.
    faceDisp.setSize(pp.size());
    faceDisp = Zero;
    faceSurfaceNormal.setSize(pp.size());
    faceSurfaceNormal = Zero;
    faceSurfaceGlobalRegion.setSize(pp.size());
    faceSurfaceGlobalRegion = -1;

    // Divide surfaces into zoned and unzoned
    const labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones());
    const labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones());

    // Per pp face the current surface snapped to
    labelList snapSurf(pp.size(), -1);


    // Do zoned surfaces
    // ~~~~~~~~~~~~~~~~~
    // Zoned faces only attract to corresponding surface

    // Extract faces per zone
    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfi = zonedSurfaces[i];

        const wordList& faceZoneNames = surfZones[zoneSurfi].faceZoneNames();

        // Get indices of faces on pp that are also in zone
        DynamicList<label> ppFaces;
        DynamicList<label> meshFaces;
        DynamicField<point> fc;
        forAll(faceZoneNames, fzi)
        {
            const word& faceZoneName = faceZoneNames[fzi];
            label zonei = mesh.faceZones().findZoneID(faceZoneName);
            if (zonei == -1)
            {
                FatalErrorInFunction
                    << "Problem. Cannot find zone " << faceZoneName
                    << exit(FatalError);
            }
            const faceZone& fZone = mesh.faceZones()[zonei];
            const bitSet isZonedFace(mesh.nFaces(), fZone);

            ppFaces.reserve(ppFaces.capacity()+fZone.size());
            meshFaces.reserve(meshFaces.capacity()+fZone.size());
            fc.reserve(meshFaces.capacity()+fZone.size());

            forAll(pp.addressing(), i)
            {
                if (isZonedFace[pp.addressing()[i]])
                {
                    snapSurf[i] = zoneSurfi;
                    ppFaces.append(i);
                    meshFaces.append(pp.addressing()[i]);
                    fc.append(ppFaceCentres[i]);
                }
            }

            //Pout<< "For faceZone " << fZone.name()
            //    << " found " << ppFaces.size() << " out of " << pp.size()
            //    << endl;
        }

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            labelList(1, zoneSurfi),
            fc,
            sqr(scalarField(faceSnapDist, ppFaces)),// sqr of attract dist
            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, hiti)
        {
            if (hitInfo[hiti].hit())
            {
                label facei = ppFaces[hiti];
                faceDisp[facei] = hitInfo[hiti].point() - fc[hiti];
                faceSurfaceNormal[facei] = hitNormal[hiti];
                faceSurfaceGlobalRegion[facei] = surfaces.globalRegion
                (
                    hitSurface[hiti],
                    hitRegion[hiti]
                );
            }
            else
            {
                static label nWarn = 0;

                if (nWarn < 100)
                {
                    WarningInFunction
                        << "Did not find surface near face centre " << fc[hiti]
                    << endl;
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
    }


    // Do unzoned surfaces
    // ~~~~~~~~~~~~~~~~~~~
    // Unzoned faces attract to any unzoned surface

    DynamicList<label> ppFaces(pp.size());
    DynamicList<label> meshFaces(pp.size());
    DynamicField<point> fc(pp.size());
    forAll(pp.addressing(), i)
    {
        if (snapSurf[i] == -1)
        {
            ppFaces.append(i);
            meshFaces.append(pp.addressing()[i]);
            fc.append(ppFaceCentres[i]);
        }
    }
    //Pout<< "Found " << ppFaces.size() << " unzoned faces out of "
    //   << pp.size() << endl;

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    labelList hitRegion;
    vectorField hitNormal;
    surfaces.findNearestRegion
    (
        unzonedSurfaces,
        fc,
        sqr(scalarField(faceSnapDist, ppFaces)),// sqr of attract dist
        hitSurface,
        hitInfo,
        hitRegion,
        hitNormal
    );

    forAll(hitInfo, hiti)
    {
        if (hitInfo[hiti].hit())
        {
            label facei = ppFaces[hiti];
            faceDisp[facei] = hitInfo[hiti].point() - fc[hiti];
            faceSurfaceNormal[facei] = hitNormal[hiti];
            faceSurfaceGlobalRegion[facei] = surfaces.globalRegion
            (
                hitSurface[hiti],
                hitRegion[hiti]
            );
        }
        else
        {
            static label nWarn = 0;

            if (nWarn < 100)
            {
                WarningInFunction
                    << "Did not find surface near face centre " << fc[hiti]
                    << endl;

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


    //// Determine rotation
    //// ~~~~~~~~~~~~~~~~~~
    //
    //// Determine rotation axis
    //faceRotation.setSize(pp.size());
    //faceRotation = Zero;
    //
    //forAll(faceRotation, facei)
    //{
    //    // Note: extend to >180 degrees checking
    //    faceRotation[facei] =
    //        pp.faceNormals()[facei]
    //      ^ faceSurfaceNormal[facei];
    //}
    //
    //if (debug&meshRefinement::ATTRACTION)
    //{
    //    dumpMove
    //    (
    //        mesh.time().path()
    //      / "faceDisp_" + name(iter) + ".obj",
    //        pp.faceCentres(),
    //        pp.faceCentres() + faceDisp
    //    );
    //    dumpMove
    //    (
    //        mesh.time().path()
    //      / "faceRotation_" + name(iter) + ".obj",
    //        pp.faceCentres(),
    //        pp.faceCentres() + faceRotation
    //    );
    //}
}


// Collect (possibly remote) per point data of all surrounding faces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// - faceSurfaceNormal
// - faceDisp
// - faceCentres&faceNormal
void Foam::snappySnapDriver::calcNearestFacePointProperties
(
    const label iter,
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,

    const vectorField& faceDisp,
    const vectorField& faceSurfaceNormal,
    const labelList& faceSurfaceGlobalRegion,

    List<List<point>>& pointFaceSurfNormals,
    List<List<point>>& pointFaceDisp,
    List<List<point>>& pointFaceCentres,
    List<labelList>&  pointFacePatchID
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const bitSet isMasterFace(syncTools::getMasterFaces(mesh));


    // For now just get all surrounding face data. Expensive - should just
    // store and sync data on coupled points only
    // (see e.g PatchToolsNormals.C)

    pointFaceSurfNormals.setSize(pp.nPoints());
    pointFaceDisp.setSize(pp.nPoints());
    pointFaceCentres.setSize(pp.nPoints());
    pointFacePatchID.setSize(pp.nPoints());

    // Fill local data
    forAll(pp.pointFaces(), pointi)
    {
        const labelList& pFaces = pp.pointFaces()[pointi];

        // Count valid face normals
        label nFaces = 0;
        forAll(pFaces, i)
        {
            label facei = pFaces[i];
            label globalRegioni = faceSurfaceGlobalRegion[facei];

            if (isMasterFace[pp.addressing()[facei]] && globalRegioni != -1)
            {
                nFaces++;
            }
        }


        List<point>& pNormals = pointFaceSurfNormals[pointi];
        pNormals.setSize(nFaces);
        List<point>& pDisp = pointFaceDisp[pointi];
        pDisp.setSize(nFaces);
        List<point>& pFc = pointFaceCentres[pointi];
        pFc.setSize(nFaces);
        labelList& pFid = pointFacePatchID[pointi];
        pFid.setSize(nFaces);

        nFaces = 0;
        forAll(pFaces, i)
        {
            label facei = pFaces[i];
            label globalRegioni = faceSurfaceGlobalRegion[facei];

            if (isMasterFace[pp.addressing()[facei]] && globalRegioni != -1)
            {
                pNormals[nFaces] = faceSurfaceNormal[facei];
                pDisp[nFaces] = faceDisp[facei];
                pFc[nFaces] = pp.localFaces()[facei].centre(ppLocalPoints);
                pFid[nFaces] = globalToMasterPatch_[globalRegioni];
                nFaces++;
            }
        }
    }


    // Collect additionally 'normal' boundary faces for boundaryPoints of pp
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // points on the boundary of pp should pick up non-pp normals
    // as well for the feature-reconstruction to behave correctly.
    // (the movement is already constrained outside correctly so it
    //  is only that the unconstrained attraction vector is calculated
    //  correctly)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList patchID(pbm.patchID());

        // Unmark all non-coupled boundary faces
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (pp.coupled() || isA<emptyPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label meshFacei = pp.start()+i;
                    patchID[meshFacei-mesh.nInternalFaces()] = -1;
                }
            }
        }

        // Remove any meshed faces
        forAll(pp.addressing(), i)
        {
            label meshFacei = pp.addressing()[i];
            patchID[meshFacei-mesh.nInternalFaces()] = -1;
        }



        // See if edge of pp uses any non-meshed boundary faces. If so add the
        // boundary face as additional constraint. Note that we account for
        // both 'real' boundary edges and boundary edge of baffles

        const labelList bafflePair
        (
            localPointRegion::findDuplicateFaces(mesh, pp.addressing())
        );


        // Mark all points on 'boundary' edges
        bitSet isBoundaryPoint(pp.nPoints());

        const labelListList& edgeFaces = pp.edgeFaces();
        const edgeList& edges = pp.edges();

        forAll(edgeFaces, edgei)
        {
            const edge& e = edges[edgei];
            const labelList& eFaces = edgeFaces[edgei];

            if (eFaces.size() == 1)
            {
                // 'real' boundary edge
                isBoundaryPoint.set(e[0]);
                isBoundaryPoint.set(e[1]);
            }
            else if (eFaces.size() == 2 && bafflePair[eFaces[0]] == eFaces[1])
            {
                // 'baffle' boundary edge
                isBoundaryPoint.set(e[0]);
                isBoundaryPoint.set(e[1]);
            }
        }


        // Construct labelList equivalent of meshPointMap
        labelList meshToPatchPoint(mesh.nPoints(), -1);
        forAll(pp.meshPoints(), pointi)
        {
            meshToPatchPoint[pp.meshPoints()[pointi]] = pointi;
        }

        forAll(patchID, bFacei)
        {
            label patchi = patchID[bFacei];

            if (patchi != -1)
            {
                label facei = mesh.nInternalFaces()+bFacei;
                const face& f = mesh.faces()[facei];

                forAll(f, fp)
                {
                    label pointi = meshToPatchPoint[f[fp]];

                    if (pointi != -1 && isBoundaryPoint.test(pointi))
                    {
                        List<point>& pNormals = pointFaceSurfNormals[pointi];
                        List<point>& pDisp = pointFaceDisp[pointi];
                        List<point>& pFc = pointFaceCentres[pointi];
                        labelList& pFid = pointFacePatchID[pointi];


                        // Note: these are points which are on the patch but
                        //       the face itself is not. Could recalculate the
                        //       facearea with the modified point but hopefully
                        //       not necessary.
                        //const point& pt = mesh.points()[f[fp]];
                        const point& pt = ppLocalPoints[pointi];
                        vector fn = mesh.faceAreas()[facei];

                        pNormals.append(fn/mag(fn));
                        pDisp.append(mesh.faceCentres()[facei]-pt);
                        pFc.append(mesh.faceCentres()[facei]);
                        pFid.append(patchi);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceSurfNormals,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceDisp,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );

    {
        // Make into displacement before synchronising to avoid any problems
        // with parallel cyclics
        //pointField localPoints(pp.points(), pp.meshPoints());
        forAll(pointFaceCentres, pointi)
        {
            //const point& pt = pp.points()[pp.meshPoints()[pointi]];
            const point& pt = ppLocalPoints[pointi];

            List<point>& pFc = pointFaceCentres[pointi];
            for (point& p : pFc)
            {
                p -= pt;
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pointFaceCentres,
            listPlusEqOp<point>(),
            List<point>(),
            mapDistribute::transform()
        );
        forAll(pointFaceCentres, pointi)
        {
            //const point& pt = pp.points()[pp.meshPoints()[pointi]];
            const point& pt = ppLocalPoints[pointi];

            List<point>& pFc = pointFaceCentres[pointi];
            for (point& p : pFc)
            {
                p += pt;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFacePatchID,
        listPlusEqOp<label>(),
        List<label>()
    );


    // Sort the data according to the face centres. This is only so we get
    // consistent behaviour serial and parallel.
    labelList visitOrder;
    forAll(pointFaceDisp, pointi)
    {
        List<point>& pNormals = pointFaceSurfNormals[pointi];
        List<point>& pDisp = pointFaceDisp[pointi];
        List<point>& pFc = pointFaceCentres[pointi];
        labelList& pFid = pointFacePatchID[pointi];

        sortedOrder(mag(pFc)(), visitOrder);

        pNormals = List<point>(pNormals, visitOrder);
        pDisp = List<point>(pDisp, visitOrder);
        pFc = List<point>(pFc, visitOrder);
        pFid = labelUIndList(pFid, visitOrder)();
    }
}


// Gets passed in offset to nearest point on feature edge. Calculates
// if the point has a different number of faces on either side of the feature
// and if so attracts the point to that non-dominant plane.
void Foam::snappySnapDriver::correctAttraction
(
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceCounts,
    const point& edgePt,
    const vector& edgeNormal,       // normalised normal
    const point& pt,

    vector& edgeOffset              // offset from pt to point on edge
) const
{
    // Tangential component along edge
    scalar tang = ((pt-edgePt)&edgeNormal);

    labelList order(sortedOrder(surfaceCounts));

    if (order[0] < order[1])
    {
        // There is a non-dominant plane. Use the point on the plane to
        // attract to.
        vector attractD = surfacePoints[order[0]]-edgePt;
        // Tangential component along edge
        scalar tang2 = (attractD&edgeNormal);
        // Normal component
        attractD -= tang2*edgeNormal;
        // Calculate fraction of normal distances
        scalar magAttractD = mag(attractD);
        scalar fraction = magAttractD/(magAttractD+mag(edgeOffset));

        point linePt =
            edgePt
          + ((1.0-fraction)*tang2 + fraction*tang)*edgeNormal;
        edgeOffset = linePt-pt;
    }
}


Foam::pointIndexHit Foam::snappySnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const List<point>& faceCentres
) const
{
    // Determine if multiple patchIDs
    if (patchIDs.size())
    {
        label patch0 = patchIDs[0];

        for (label i = 1; i < patchIDs.size(); i++)
        {
            if (patchIDs[i] != patch0)
            {
                return pointIndexHit(true, pt, labelMax);
            }
        }
    }
    return pointIndexHit(false, Zero, labelMax);
}


Foam::label Foam::snappySnapDriver::findNormal
(
    const scalar featureCos,
    const vector& n,
    const DynamicList<vector>& surfaceNormals
) const
{
    label index = -1;

    forAll(surfaceNormals, j)
    {
        scalar cosAngle = (n&surfaceNormals[j]);

        if
        (
            (cosAngle >= featureCos)
         || (cosAngle < (-1+0.001)) // triangle baffles
        )
        {
            index = j;
            break;
        }
    }
    return index;
}


// Detect multiple patches. Returns pointIndexHit:
// - false, index=-1 : single patch
// - true , index=0  : multiple patches but on different normals planes
//                     (so geometric feature edge is also a region edge)
// - true , index=1  : multiple patches on same normals plane i.e. flat region
//                     edge
Foam::pointIndexHit Foam::snappySnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const DynamicList<vector>& surfaceNormals,
    const labelList& faceToNormalBin
) const
{
    if (patchIDs.empty())
    {
        return pointIndexHit(false, pt, -1);
    }

    // Detect single patch situation (to avoid allocation)
    label patch0 = patchIDs[0];

    for (label i = 1; i < patchIDs.size(); i++)
    {
        if (patchIDs[i] != patch0)
        {
            patch0 = -1;
            break;
        }
    }

    if (patch0 >= 0)
    {
        // Single patch
        return pointIndexHit(false, pt, -1);
    }
    else
    {
        if (surfaceNormals.size() == 1)
        {
            // Same normals plane, flat region edge.
            return pointIndexHit(true, pt, 1);
        }
        else
        {
            // Detect per normals bin
            labelList normalToPatch(surfaceNormals.size(), -1);
            forAll(faceToNormalBin, i)
            {
                if (faceToNormalBin[i] != -1)
                {
                    label& patch = normalToPatch[faceToNormalBin[i]];
                    if (patch == -1)
                    {
                        // First occurrence
                        patch = patchIDs[i];
                    }
                    else if (patch == -2)
                    {
                        // Already marked as being on multiple patches
                    }
                    else if (patch != patchIDs[i])
                    {
                        // Mark as being on multiple patches
                        patch = -2;
                    }
                }
            }

            forAll(normalToPatch, normali)
            {
                if (normalToPatch[normali] == -2)
                {
                    // Multiple patches on same normals plane, flat region
                    // edge
                    return pointIndexHit(true, pt, 1);
                }
            }

            // All patches on either side of geometric feature anyway
            return pointIndexHit(true, pt, 0);
        }
    }
}


void Foam::snappySnapDriver::writeStats
(
    const indirectPrimitivePatch& pp,
    const bitSet& isPatchMasterPoint,
    const List<pointConstraint>& patchConstraints
) const
{
    label nMasterPoints = 0;
    label nPlanar = 0;
    label nEdge = 0;
    label nPoint = 0;

    forAll(patchConstraints, pointi)
    {
        if (isPatchMasterPoint[pointi])
        {
            nMasterPoints++;

            if (patchConstraints[pointi].first() == 1)
            {
                nPlanar++;
            }
            else if (patchConstraints[pointi].first() == 2)
            {
                nEdge++;
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                nPoint++;
            }
        }
    }

    reduce(nMasterPoints, sumOp<label>());
    reduce(nPlanar, sumOp<label>());
    reduce(nEdge, sumOp<label>());
    reduce(nPoint, sumOp<label>());
    Info<< "total master points :" << nMasterPoints
        << " of which attracted to :" << nl
        << "    feature point   : " << nPoint << nl
        << "    feature edge    : " << nEdge << nl
        << "    nearest surface : " << nPlanar << nl
        << "    rest            : " << nMasterPoints-nPoint-nEdge-nPlanar
        << nl
        << endl;
}


void Foam::snappySnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const label pointi,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    labelList& faceToNormalBin,

    vector& patchAttraction,
    pointConstraint& patchConstraint
) const
{
    patchAttraction = Zero;
    patchConstraint = pointConstraint();

    const List<point>& pfSurfNormals = pointFaceSurfNormals[pointi];
    const List<point>& pfDisp = pointFaceDisp[pointi];
    const List<point>& pfCentres = pointFaceCentres[pointi];

    // Bin according to surface normal
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Bins of differing normals:
    //  - one normal   : flat(tish) surface
    //  - two normals  : geometric feature edge
    //  - three normals: geometric feature point
    //  - four normals : too complex a feature
    surfacePoints.clear();
    surfaceNormals.clear();

    //- From face to above normals bin
    faceToNormalBin.setSize(pfDisp.size());
    faceToNormalBin = -1;

    forAll(pfSurfNormals, i)
    {
        const point& fc = pfCentres[i];
        const vector& fSNormal = pfSurfNormals[i];
        const vector& fDisp = pfDisp[i];

        // What to do with very far attraction? For now just ignore the face
        if (magSqr(fDisp) < sqr(snapDist[pointi]) && mag(fSNormal) > VSMALL)
        {
            const point pt = fc + fDisp;

            // Do we already have surface normal?
            faceToNormalBin[i] = findNormal
            (
                featureCos,
                fSNormal,
                surfaceNormals
            );

            if (faceToNormalBin[i] != -1)
            {
                // Same normal
            }
            else
            {
                // Now check if the planes go through the same edge or point

                if (surfacePoints.size() <= 1)
                {
                    surfacePoints.append(pt);
                    faceToNormalBin[i] = surfaceNormals.size();
                    surfaceNormals.append(fSNormal);
                }
                else if (surfacePoints.size() == 2)
                {
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());

                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        // Definitely makes a feature point
                        surfacePoints.append(pt);
                        faceToNormalBin[i] = surfaceNormals.size();
                        surfaceNormals.append(fSNormal);
                    }
                }
                else if (surfacePoints.size() == 3)
                {
                    // Have already feature point. See if this new plane is
                    // the same point or not.
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane pl2(surfacePoints[2], surfaceNormals[2]);
                    point p012(pl0.planePlaneIntersect(pl1, pl2));

                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());
                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        plane pl3(pt, fSNormal);
                        point p013(pl0.planePlaneIntersect(pl1, pl3));

                        if (mag(p012-p013) > snapDist[pointi])
                        {
                            // Different feature point
                            surfacePoints.append(pt);
                            faceToNormalBin[i] = surfaceNormals.size();
                            surfaceNormals.append(fSNormal);
                        }
                    }
                }
            }
        }
    }


    //const point& pt = pp.localPoints()[pointi];
    const point& pt = ppLocalPoints[pointi];

    // Check the number of directions
    if (surfaceNormals.size() == 1)
    {
        // Normal distance to plane
        vector d =
            ((surfacePoints[0]-pt) & surfaceNormals[0])
           *surfaceNormals[0];

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
    }
    else if (surfaceNormals.size() == 2)
    {
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane::ray r(pl0.planeIntersect(pl1));
        vector n = r.dir() / mag(r.dir());

        // Get nearest point on infinite ray
        vector d = r.refPoint()-pt;
        d -= (d&n)*n;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
    }
    else if (surfaceNormals.size() == 3)
    {
        // Calculate point from the faces.
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane pl2(surfacePoints[2], surfaceNormals[2]);
        point cornerPt(pl0.planePlaneIntersect(pl1, pl2));
        vector d = cornerPt - pt;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
        patchConstraint.applyConstraint(surfaceNormals[2]);
    }
}


// Special version that calculates attraction in one go
void Foam::snappySnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> feStr;
    autoPtr<OBJstream> fpStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        feStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-edge direction to "
            << feStr().name() << endl;

        fpStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeaturePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-point direction to "
            << fpStr().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(ppLocalPoints, pointi)
    {
        vector attraction = Zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,
            nearestDisp,

            pointi,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );

        if
        (
            (constraint.first() > patchConstraints[pointi].first())
         || (
                (constraint.first() == patchConstraints[pointi].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointi]))
            )
        )
        {
            patchAttraction[pointi] = attraction;
            patchConstraints[pointi] = constraint;

            //const point& pt = pp.localPoints()[pointi];
            const point& pt = ppLocalPoints[pointi];

            if (feStr && patchConstraints[pointi].first() == 2)
            {
                feStr().writeLine(pt, pt+patchAttraction[pointi]);
            }
            else if (fpStr && patchConstraints[pointi].first() == 3)
            {
                fpStr().writeLine(pt, pt+patchAttraction[pointi]);
            }
        }
    }
}


void Foam::snappySnapDriver::stringFeatureEdges
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)

    // What this does is fill in any faces where not all points
    // on the face are being attracted:
    /*
           +
          / \
         /   \
      ---+    +---
         \   /
          \ /
           +
    */
    // so the top and bottom will never get attracted since the nearest
    // back from the feature edge will always be one of the left or right
    // points since the face is diamond like. So here we walk the feature edges
    // and add any non-attracted points.


    while (true)
    {
        label nChanged = 0;

        const labelListList& pointEdges = pp.pointEdges();
        forAll(pointEdges, pointi)
        {
            if (patchConstraints[pointi].first() == 2)
            {
                //const point& pt = pp.localPoints()[pointi];
                const point& pt = ppLocalPoints[pointi];
                const labelList& pEdges = pointEdges[pointi];
                const vector& featVec = patchConstraints[pointi].second();

                // Detect whether there are edges in both directions.
                // (direction along the feature edge that is)
                bool hasPos = false;
                bool hasNeg = false;

                forAll(pEdges, pEdgei)
                {
                    const edge& e = pp.edges()[pEdges[pEdgei]];
                    label nbrPointi = e.otherVertex(pointi);

                    if (patchConstraints[nbrPointi].first() > 1)
                    {
                        //const point& nbrPt = pp.localPoints()[nbrPointi];
                        const point& nbrPt = ppLocalPoints[nbrPointi];
                        const point featPt =
                            nbrPt + patchAttraction[nbrPointi];
                        const scalar cosAngle = (featVec & (featPt-pt));

                        if (cosAngle > 0)
                        {
                            hasPos = true;
                        }
                        else
                        {
                            hasNeg = true;
                        }
                    }
                }

                if (!hasPos || !hasNeg)
                {
                    //Pout<< "**Detected feature string end at  "
                    //    << ppLocalPoints[pointi] << endl;

                    // No string. Assign best choice on either side
                    label bestPosPointi = -1;
                    scalar minPosDistSqr = GREAT;
                    label bestNegPointi = -1;
                    scalar minNegDistSqr = GREAT;

                    forAll(pEdges, pEdgei)
                    {
                        const edge& e = pp.edges()[pEdges[pEdgei]];
                        label nbrPointi = e.otherVertex(pointi);

                        if
                        (
                            patchConstraints[nbrPointi].first() <= 1
                         && rawPatchConstraints[nbrPointi].first() > 1
                        )
                        {
                            const vector& nbrFeatVec =
                                rawPatchConstraints[pointi].second();

                            if (mag(featVec&nbrFeatVec) > featureCos)
                            {
                                // nbrPointi attracted to sameish feature
                                // Note: also check on position.

                                scalar d2 = magSqr
                                (
                                    rawPatchAttraction[nbrPointi]
                                );

                                const point featPt =
                                    //pp.localPoints()[nbrPointi]
                                    ppLocalPoints[nbrPointi]
                                  + rawPatchAttraction[nbrPointi];
                                const scalar cosAngle =
                                    (featVec & (featPt-pt));

                                if (cosAngle > 0)
                                {
                                    if (!hasPos && d2 < minPosDistSqr)
                                    {
                                        minPosDistSqr = d2;
                                        bestPosPointi = nbrPointi;
                                    }
                                }
                                else
                                {
                                    if (!hasNeg && d2 < minNegDistSqr)
                                    {
                                        minNegDistSqr = d2;
                                        bestNegPointi = nbrPointi;
                                    }
                                }
                            }
                        }
                    }

                    if (bestPosPointi != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    ppLocalPoints[bestPosPointi];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestPosPointi]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestPosPointi] =
                            0.5*rawPatchAttraction[bestPosPointi];
                        patchConstraints[bestPosPointi] =
                            rawPatchConstraints[bestPosPointi];

                        nChanged++;
                    }
                    if (bestNegPointi != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    ppLocalPoints[bestNegPointi];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestNegPointi]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestNegPointi] =
                            0.5*rawPatchAttraction[bestNegPointi];
                        patchConstraints[bestNegPointi] =
                            rawPatchConstraints[bestNegPointi];

                        nChanged++;
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        Info<< "Stringing feature edges : changed " << nChanged << " points"
            << endl;
        if (nChanged == 0)
        {
            break;
        }
    }
}


void Foam::snappySnapDriver::releasePointsNextToMultiPatch
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,

    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> multiPatchStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        multiPatchStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "multiPatch_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping removed constraints due to same-face"
            << " multi-patch points to "
            << multiPatchStr().name() << endl;
    }


    // 1. Mark points on multiple patches
    bitSet isMultiPatchPoint(pp.size());

    forAll(pointFacePatchID, pointi)
    {
        pointIndexHit multiPatchPt = findMultiPatchPoint
        (
            ppLocalPoints[pointi],
            pointFacePatchID[pointi],
            pointFaceCentres[pointi]
        );
        isMultiPatchPoint.set(pointi, multiPatchPt.hit());
    }

    // 2. Make sure multi-patch points are also attracted
    forAll(isMultiPatchPoint, pointi)
    {
        if (isMultiPatchPoint.test(pointi))
        {
            if
            (
                patchConstraints[pointi].first() <= 1
             && rawPatchConstraints[pointi].first() > 1
            )
            {
                patchAttraction[pointi] = rawPatchAttraction[pointi];
                patchConstraints[pointi] = rawPatchConstraints[pointi];

                //if (multiPatchStr)
                //{
                //    Pout<< "Adding constraint on multiPatchPoint:"
                //        << ppLocalPoints[pointi]
                //        << " constraint:" << patchConstraints[pointi]
                //        << " attraction:" << patchAttraction[pointi]
                //        << endl;
                //}
            }
        }
    }

    // Up to here it is all parallel ok.


    // 3. Knock out any attraction on faces with multi-patch points
    label nChanged = 0;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        label nMultiPatchPoints = 0;
        forAll(f, fp)
        {
            label pointi = f[fp];
            if
            (
                isMultiPatchPoint.test(pointi)
             && patchConstraints[pointi].first() > 1
            )
            {
                ++nMultiPatchPoints;
            }
        }

        if (nMultiPatchPoints > 0)
        {
            forAll(f, fp)
            {
                label pointi = f[fp];
                if
                (
                   !isMultiPatchPoint.test(pointi)
                 && patchConstraints[pointi].first() > 1
                )
                {
                    //Pout<< "Knocking out constraint"
                    //    << " on non-multiPatchPoint:"
                    //    << ppLocalPoints[pointi] << endl;
                    patchAttraction[pointi] = Zero;
                    patchConstraints[pointi] = pointConstraint();
                    nChanged++;

                    if (multiPatchStr)
                    {
                        multiPatchStr().write(ppLocalPoints[pointi]);
                    }
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());
    Info<< "Removing constraints near multi-patch points : changed "
        << nChanged << " points" << endl;
}


Foam::labelPair Foam::snappySnapDriver::findDiagonalAttraction
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchAttraction,
    const List<pointConstraint>& patchConstraints,
    const label facei
) const
{
    const face& f = pp.localFaces()[facei];
    // For now just detect any attraction. Improve this to look at
    // actual attraction position and orientation

    labelPair attractIndices(-1, -1);

    if (f.size() >= 4)
    {
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
                if
                (
                    patchConstraints[f[startFp]].first() >= 2
                 && patchConstraints[f[endFp]].first() >= 2
                )
                {
                    attractIndices = labelPair(startFp, endFp);
                    break;
                }
            }
        }
    }
    return attractIndices;
}


bool Foam::snappySnapDriver::isSplitAlignedWithFeature
(
    const scalar featureCos,
    const point& p0,
    const pointConstraint& pc0,
    const point& p1,
    const pointConstraint& pc1
) const
{
    vector d(p1-p0);
    scalar magD = mag(d);
    if (magD < VSMALL)
    {
        // Two diagonal points already colocated?
        return false;
    }
    else
    {
        d /= magD;

        // Is diagonal d aligned with at least one of the feature
        // edges?

        if (pc0.first() == 2 && mag(d & pc0.second()) > featureCos)
        {
            return true;
        }
        else if (pc1.first() == 2 && mag(d & pc1.second()) > featureCos)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


// Is situation very concave
bool Foam::snappySnapDriver::isConcave
(
    const point& c0,
    const vector& area0,
    const point& c1,
    const vector& area1,
    const scalar concaveCos
) const
{
    vector n0 = area0;
    scalar magN0 = mag(n0);
    if (magN0 < VSMALL)
    {
        // Zero area face. What to return? For now disable splitting.
        return true;
    }
    n0 /= magN0;

    // Distance from c1 to plane of face0
    scalar d = (c1-c0)&n0;

    if (d <= 0)
    {
        // Convex (face1 centre on 'inside' of face0)
        return false;
    }
    else
    {
        // Is a bit or very concave?
        vector n1 = area1;
        scalar magN1 = mag(n1);
        if (magN1 < VSMALL)
        {
            // Zero area face. See above
            return true;
        }
        n1 /= magN1;

        if ((n0&n1) < concaveCos)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


Foam::labelPair Foam::snappySnapDriver::findDiagonalAttraction
(
    const scalar featureCos,
    const scalar concaveCos,
    const scalar minAreaRatio,
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const vectorField& patchAttr,
    const List<pointConstraint>& patchConstraints,
    const vectorField& nearestAttr,
    const vectorField& nearestNormal,
    const label facei,

    DynamicField<point>& points0,
    DynamicField<point>& points1
) const
{
    const face& localF = pp.localFaces()[facei];

    labelPair attractIndices(-1, -1);

    if (localF.size() >= 4)
    {
        //// Estimate cell centre taking patchAttraction into account
        //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //// (is this necessary?)
        //const polyMesh& mesh = meshRefiner_.mesh();
        //label meshFacei = pp.addressing()[facei];
        //const face& meshF = mesh.faces()[meshFacei];
        //label celli = mesh.faceOwner()[meshFacei];
        //const labelList& cPoints = mesh.cellPoints(celli);
        //
        //point cc(mesh.points()[meshF[0]]);
        //for (label i = 1; i < meshF.size(); i++)
        //{
        //    cc += mesh.points()[meshF[i]]+patchAttr[localF[i]];
        //}
        //forAll(cPoints, i)
        //{
        //    label pointi = cPoints[i];
        //    if (!meshF.found(pointi))
        //    {
        //        cc += mesh.points()[pointi];
        //    }
        //}
        //cc /= cPoints.size();
        ////const point& cc = mesh.cellCentres()[celli];
        //
        //const scalar vol = pyrVol(pp, patchAttr, localF, cc);
        //const scalar area = localF.mag(localPts);



        // Try all diagonal cuts
        // ~~~~~~~~~~~~~~~~~~~~~

        face f0(3);
        face f1(3);

        for (label startFp = 0; startFp < localF.size()-2; startFp++)
        {
            label minFp = localF.rcIndex(startFp);

            for
            (
                label endFp = localF.fcIndex(localF.fcIndex(startFp));
                endFp < localF.size() && endFp != minFp;
                endFp++
            )
            {
                label startPti = localF[startFp];
                label endPti = localF[endFp];

                const pointConstraint& startPc = patchConstraints[startPti];
                const pointConstraint& endPc = patchConstraints[endPti];

                if (startPc.first() >= 2 && endPc.first() >= 2)
                {
                    if (startPc.first() == 2 || endPc.first() == 2)
                    {
                        // Check if
                        // - sameish feature edge normal
                        // - diagonal aligned with feature edge normal
                        point start =
                            ppLocalPoints[startPti]+patchAttr[startPti];
                        point end =
                            ppLocalPoints[endPti]+patchAttr[endPti];

                        if
                        (
                           !isSplitAlignedWithFeature
                            (
                                featureCos,
                                start,
                                startPc,
                                end,
                                endPc
                            )
                        )
                        {
                            // Attract to different features. No need to
                            // introduce split
                            continue;
                        }
                    }



                    // Form two faces
                    // ~~~~~~~~~~~~~~
                    // Predict position of faces. End points of the faces
                    // attract to the feature
                    // and all the other points just attract to the nearest

                    // face0

                    f0.setSize(endFp-startFp+1);
                    label i0 = 0;
                    for (label fp = startFp; fp <= endFp; fp++)
                    {
                        f0[i0++] = localF[fp];
                    }

                    // Get compact face and points
                    const face compact0(identity(f0.size()));
                    points0.clear();
                    points0.append(ppLocalPoints[f0[0]] + patchAttr[f0[0]]);
                    for (label fp=1; fp < f0.size()-1; fp++)
                    {
                        label pi = f0[fp];
                        points0.append(ppLocalPoints[pi] + nearestAttr[pi]);
                    }
                    points0.append
                    (
                        ppLocalPoints[f0.last()] + patchAttr[f0.last()]
                    );


                    // face1

                    f1.setSize(localF.size()+2-f0.size());
                    label i1 = 0;
                    for
                    (
                        label fp = endFp;
                        fp != startFp;
                        fp = localF.fcIndex(fp)
                    )
                    {
                        f1[i1++] = localF[fp];
                    }
                    f1[i1++] = localF[startFp];


                    // Get compact face and points
                    const face compact1(identity(f1.size()));
                    points1.clear();
                    points1.append(ppLocalPoints[f1[0]] + patchAttr[f1[0]]);
                    for (label fp=1; fp < f1.size()-1; fp++)
                    {
                        label pi = f1[fp];
                        points1.append(ppLocalPoints[pi] + nearestAttr[pi]);
                    }
                    points1.append
                    (
                        ppLocalPoints[f1.last()] + patchAttr[f1.last()]
                    );



                    // Avoid splitting concave faces
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if
                    (
                        isConcave
                        (
                            compact0.centre(points0),
                            compact0.areaNormal(points0),
                            compact1.centre(points1),
                            compact1.areaNormal(points1),
                            concaveCos
                        )
                    )
                    {
                        //// Dump to obj
                        //Pout<< "# f0:" << f0 << endl;
                        //forAll(p, i)
                        //{
                        //    meshTools::writeOBJ(Pout, points0[i]);
                        //}
                        //Pout<< "# f1:" << f1 << endl;
                        //forAll(p, i)
                        //{
                        //    meshTools::writeOBJ(Pout, points1[i]);
                        //}
                    }
                    else
                    {
                        // Existing areas
                        const scalar area0 = f0.mag(ppLocalPoints);
                        const scalar area1 = f1.mag(ppLocalPoints);

                        if
                        (
                            area0/area1 >= minAreaRatio
                         && area1/area0 >= minAreaRatio
                        )
                        {
                            attractIndices = labelPair(startFp, endFp);
                        }
                    }
                }
            }
        }
    }
    return attractIndices;
}


void Foam::snappySnapDriver::splitDiagonals
(
    const scalar featureCos,
    const scalar concaveCos,
    const scalar minAreaRatio,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const vectorField& nearestAttraction,
    const vectorField& nearestNormal,
    const List<labelList>& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits,
    DynamicList<labelPair>& splitPatches
) const
{
    const labelList& bFaces = pp.addressing();

    splitFaces.clear();
    splitFaces.setCapacity(bFaces.size());
    splits.clear();
    splits.setCapacity(bFaces.size());
    splitPatches.clear();
    splitPatches.setCapacity(bFaces.size());

    // Work arrays for storing points of face
    DynamicField<point> facePoints0;
    DynamicField<point> facePoints1;
    face f0, f1;

    forAll(bFaces, facei)
    {
        const labelPair split
        (
            findDiagonalAttraction
            (
                featureCos,
                concaveCos,
                minAreaRatio,

                pp,
                ppLocalPoints,
                patchAttraction,
                patchConstraints,

                nearestAttraction,
                nearestNormal,
                facei,

                facePoints0,
                facePoints1
            )
        );

        if (split != labelPair(-1, -1))
        {
            splitFaces.append(bFaces[facei]);
            splits.append(split);
            splitPatches.append(labelPair(-1, -1));

            const face& f = pp.localFaces()[facei];

            // Knock out other attractions on face
            forAll(f, fp)
            {
                if
                (
                    fp != split[0]
                 && fp != split[1]
                 && patchConstraints[f[fp]].first() >= 2
                )
                {
                    patchConstraints[f[fp]] = pointConstraint
                    (
                        Tuple2<label, vector>
                        (
                            1,
                            nearestNormal[f[fp]]
                        )
                    );
                    patchAttraction[f[fp]] = nearestAttraction[f[fp]];
                }
            }


            // Detect any patches to give to the two faces
            labelPair& twoPatches = splitPatches.last();

            // Split into two faces:
            //  - face0 is from split[0] up to split[1]
            //  - face1 is from split[1] up to split[0]
            meshRefinement::splitFace(f, split, f0, f1);

            // Bit heuristic:
            // - find any inbetween vertex (i.e. not on diagonal) that
            //   has single patch only
            // - if we find such a point on both sides we use these patches
            // - if only found a point on one side check the other side
            //   for any other patch

            // Find single patch
            for (label fp = 1; fp < f0.size()-1; fp++)
            {
                const auto& patches = pointFacePatchID[f0[fp]];
                if (patches.uniform())
                {
                    twoPatches[0] = patches[0];
                }
            }
            for (label fp = 1; fp < f1.size()-1; fp++)
            {
                const auto& patches = pointFacePatchID[f1[fp]];
                if (patches.uniform())
                {
                    twoPatches[1] = patches[0];
                }
            }

            if (twoPatches[0] == -1)
            {
                // Find any patch on any point that differs from twoPatches[1]
                // (can be -1 or a unform patch)

                for (label fp = 1; fp < f0.size()-1; fp++)
                {
                    const auto& patches = pointFacePatchID[f0[fp]];

                    // Pick any patch which is not the unique patch of the
                    // other face. TBD: count occurences and choose most used?
                    for (const label patchi : patches)
                    {
                        if (patchi != twoPatches[1])
                        {
                            twoPatches[0] = patchi;
                            break;
                        }
                    }
                }
            }

            if (twoPatches[1] == -1)
            {
                // Find any patch on any point that differs from twoPatches[0]
                // (can be -1 or a unform patch)

                for (label fp = 1; fp < f1.size()-1; fp++)
                {
                    const auto& patches = pointFacePatchID[f1[fp]];

                    // Pick any patch which is not the unique patch of the
                    // other face. TBD: count occurences and choose most used?
                    for (const label patchi : patches)
                    {
                        if (patchi != twoPatches[0])
                        {
                            twoPatches[1] = patchi;
                            break;
                        }
                    }
                }
            }

            // Fall-back : choose any connected patch
            if (twoPatches[0] == -1)
            {
                twoPatches[0] = pointFacePatchID[f0[1]][0];
            }
            if (twoPatches[1] == -1)
            {
                twoPatches[1] = pointFacePatchID[f1[1]][0];
            }

            //if (twoPatches != labelPair(-1, -1))
            //{
            //    Pout<< "Face:" << bFaces[facei]
            //        << " at:" << pp.faceCentres()[facei]
            //        << " split into f0:" << flatOutput(f0)
            //        << " patch:" << twoPatches[0]
            //        << " and f1:" << flatOutput(f1)
            //        << " patch:" << twoPatches[1]
            //        << endl;
            //}
        }
    }
}


void Foam::snappySnapDriver::avoidDiagonalAttraction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        labelPair diag = findDiagonalAttraction
        (
            pp,
            patchAttraction,
            patchConstraints,
            facei
        );

        if (diag[0] != -1 && diag[1] != -1)
        {
            // Found two diagonal points that being attracted.
            // For now just attract my one to the average of those.
            const label i0 = f[diag[0]];
            const point pt0 =
                ppLocalPoints[i0]+patchAttraction[i0];
            const label i1 = f[diag[1]];
            const point pt1 =
                ppLocalPoints[i1]+patchAttraction[i1];
            const point mid = 0.5*(pt0+pt1);

            const scalar cosAngle = mag
            (
                patchConstraints[i0].second()
              & patchConstraints[i1].second()
            );

            //Pout<< "Found diagonal attraction at indices:"
            //    << diag[0]
            //    << " and " << diag[1]
            //    << " with cosAngle:" << cosAngle
            //    << " mid:" << mid << endl;

            if (cosAngle > featureCos)
            {
                // The two feature edges are roughly in the same direction.
                // Add the nearest of the other points in the face as
                // attractor
                label minFp = -1;
                scalar minDistSqr = GREAT;
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (patchConstraints[pointi].first() <= 1)
                    {
                        scalar distSqr = mid.distSqr(ppLocalPoints[pointi]);
                        if (distSqr < minDistSqr)
                        {
                            minFp = fp;
                        }
                    }
                }
                if (minFp != -1)
                {
                    label minPointi = f[minFp];
                    patchAttraction[minPointi] =
                        mid-ppLocalPoints[minPointi];
                    patchConstraints[minPointi] = patchConstraints[f[diag[0]]];
                }
            }
            else
            {
                //Pout<< "Diagonal attractors at" << nl
                //    << "    pt0:" << pt0
                //    << "    constraint:"
                //    << patchConstraints[i0].second() << nl
                //    << "    pt1:" << pt1
                //    << "    constraint:"
                //    << patchConstraints[i1].second() << nl
                //    << "    make too large an angle:"
                //    <<  mag
                //        (
                //            patchConstraints[i0].second()
                //          & patchConstraints[i1].second()
                //        )
                //    << endl;
            }
        }
    }
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::snappySnapDriver::findNearFeatureEdge
(
    const bool isRegionEdge,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const label pointi,
    const point& estimatedPt,

    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearEdgeFeat;
    List<pointIndexHit> nearEdgeInfo;
    vectorField nearNormal;

    if (isRegionEdge)
    {
        features.findNearestRegionEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointi])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }
    else
    {
        features.findNearestEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointi])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }

    const pointIndexHit& nearInfo = nearEdgeInfo[0];
    label feati = nearEdgeFeat[0];

    if (nearInfo.hit())
    {
        // So we have a point on the feature edge. Use this
        // instead of our estimate from planes.
        edgeAttractors[feati][nearInfo.index()].append(nearInfo.point());
        pointConstraint c(Tuple2<label, vector>(2, nearNormal[0]));
        edgeConstraints[feati][nearInfo.index()].append(c);

        // Store for later use
        patchAttraction[pointi] = nearInfo.point()-ppLocalPoints[pointi];
        patchConstraints[pointi] = c;
    }
    return Tuple2<label, pointIndexHit>(feati, nearInfo);
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::snappySnapDriver::findNearFeaturePoint
(
    const bool isRegionPoint,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const label pointi,
    const point& estimatedPt,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Search for for featurePoints only! This ignores any non-feature points.

    labelList nearFeat;
    List<pointIndexHit> nearInfo;
    features.findNearestPoint
    (
        pointField(1, estimatedPt),
        scalarField(1, sqr(snapDist[pointi])),
        nearFeat,
        nearInfo
    );

    label feati = nearFeat[0];

    if (feati != -1)
    {
        const point& pt = ppLocalPoints[pointi];

        label featPointi = nearInfo[0].index();
        const point& featPt = nearInfo[0].hitPoint();
        scalar distSqr = featPt.distSqr(pt);

        // Check if already attracted
        label oldPointi = pointAttractor[feati][featPointi];

        if (oldPointi != -1)
        {
            // Check distance
            if (distSqr >= featPt.distSqr(ppLocalPoints[oldPointi]))
            {
                // oldPointi nearest. Keep.
                feati = -1;
                featPointi = -1;
            }
            else
            {
                // Current pointi nearer.
                pointAttractor[feati][featPointi] = pointi;
                pointConstraints[feati][featPointi].first() = 3;
                pointConstraints[feati][featPointi].second() = Zero;

                // Store for later use
                patchAttraction[pointi] = featPt-pt;
                patchConstraints[pointi] = pointConstraints[feati][featPointi];

                // Reset oldPointi to nearest on feature edge
                patchAttraction[oldPointi] = Zero;
                patchConstraints[oldPointi] = pointConstraint();

                findNearFeatureEdge
                (
                    isRegionPoint,      // search region edges only

                    pp,
                    ppLocalPoints,
                    snapDist,
                    oldPointi,
                    ppLocalPoints[oldPointi],

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );
            }
        }
        else
        {
            // Current pointi nearer.
            pointAttractor[feati][featPointi] = pointi;
            pointConstraints[feati][featPointi].first() = 3;
            pointConstraints[feati][featPointi].second() = Zero;

            // Store for later use
            patchAttraction[pointi] = featPt-pt;
            patchConstraints[pointi] = pointConstraints[feati][featPointi];
        }
    }

    return Tuple2<label, pointIndexHit>(feati, nearInfo[0]);
}


// Determines for every pp point - that is on multiple faces that form
// a feature - the nearest feature edge/point.
void Foam::snappySnapDriver::determineFeatures
(
    const label iter,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,
    const bool strictRegionFeatureSnap, // special >=3 patch points handling

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> featureEdgeStr;
    autoPtr<OBJstream> missedEdgeStr;
    autoPtr<OBJstream> featurePointStr;
    autoPtr<OBJstream> missedMP0Str;
    autoPtr<OBJstream> missedMP1Str;

    if (debug&meshRefinement::ATTRACTION)
    {
        featureEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edge sampling to "
            << featureEdgeStr().name() << endl;

        missedEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edges that are too far away to "
            << missedEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-point sampling to "
            << featurePointStr().name() << endl;

        missedMP0Str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdgeFromMPEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping region-edges that are too far away to "
            << missedMP0Str().name() << endl;

        missedMP1Str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdgeFromMPPoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping region-points that are too far away to "
            << missedMP1Str().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(ppLocalPoints, pointi)
    {
        const point& pt = ppLocalPoints[pointi];


        // Determine the geometric planes the point is (approximately) on.
        // This is returned as a
        // - attraction vector
        // - and a constraint
        //   (1: attract to surface, constraint is normal of plane
        //    2: attract to feature line, constraint is feature line direction
        //    3: attract to feature point, constraint is zero)

        vector attraction = Zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,
            nearestDisp,

            pointi,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,   // currently not used

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );


        if (strictRegionFeatureSnap)
        {
            // Bit tricky: if a point is on more than 2 points and only
            // attracted-to-surface upgrade it to attract-to-edge so it follows
            // more complex logic below. This is easier than complicating
            // the nearestFeatureEdge/Point below.

            if (constraint.first() == 1)
            {
                const auto& patches = pointFacePatchID[pointi];
                label patch1 = -1;
                for (label i = 1; i < patches.size(); i++)
                {
                    if (patches[i] != patches[0])
                    {
                        if (patch1 == -1)
                        {
                            patch1 = patches[i];
                        }
                        else if (patches[i] != patch1)
                        {
                            // >= on 2 patches. Upgrade to edge-attract
                            constraint.first() = 2;
                            //Pout<< "** 3 patch point:" << pointi
                            //    << " at:" << pt
                            //    << " patches:" << flatOutput(patches)
                            //    << " upgraded to feature-edge"
                            //    << " attracted to " << attraction << endl;
                            break;
                        }
                    }
                }
            }
        }


        // Now combine the reconstruction with the current state of the
        // point. The logic is quite complicated:
        // - the new constraint (from reconstruction) will only win if
        //   - the constraint is higher (feature-point wins from feature-edge
        //     etc.)
        //   - or the constraint is the same but the attraction distance is less
        //
        // - then this will be combined with explicit searching on the
        //   features and optionally the analysis of the patches using the
        //   point. This analysis can do three thing:
        //      - the point is not on multiple patches
        //      - the point is on multiple patches but these are also
        //        different planes (so the region feature is also a geometric
        //        feature)
        //      - the point is on multiple patches some of which are on
        //        the same plane. This is the problem one - do we assume it is
        //        an additional constraint (feat edge upgraded to region point,
        //        see below)?
        //
        //      Reconstruction  MultiRegionFeatureSnap          Attraction
        //      -------         ----------------------          -----------
        //      surface         false                           surface
        //      surface         true                            region edge
        //      feat edge       false                           feat edge
        //      feat edge       true and no planar regions      feat edge
        //      feat edge       true and yes planar regions     region point
        //      feat point      false                           feat point
        //      feat point      true                            region point


        if
        (
            (constraint.first() > patchConstraints[pointi].first())
         || (
                (constraint.first() == patchConstraints[pointi].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointi]))
            )
        )
        {
            patchAttraction[pointi] = attraction;
            patchConstraints[pointi] = constraint;

            // Check the number of directions
            if (patchConstraints[pointi].first() == 1)
            {
                // Flat surface. Check for different patchIDs
                if (multiRegionFeatureSnap)
                {
                    const point estimatedPt(pt + nearestDisp[pointi]);
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );

                    if (multiPatchPt.hit())
                    {
                        // Behave like when having two surface normals so
                        // attract to nearest feature edge (with a guess for
                        // the multipatch point as starting point)
                        Tuple2<label, pointIndexHit> nearInfo =
                        findNearFeatureEdge
                        (
                            true,                       // isRegionEdge
                            pp,
                            ppLocalPoints,
                            snapDist,
                            pointi,
                            multiPatchPt.point(),       // estimatedPt

                            edgeAttractors,
                            edgeConstraints,

                            patchAttraction,
                            patchConstraints
                        );

                        const pointIndexHit& info = nearInfo.second();
                        if (info.hit())
                        {
                            // Dump
                            if (featureEdgeStr)
                            {
                                featureEdgeStr().writeLine
                                (
                                    pt,
                                    info.point()
                                );
                            }
                        }
                        else
                        {
                            if (missedEdgeStr)
                            {
                                missedEdgeStr().writeLine
                                (
                                    pt,
                                    multiPatchPt.point()
                                );
                            }
                        }
                    }
                }
            }
            else if (patchConstraints[pointi].first() == 2)
            {
                // Mark point on the nearest feature edge. Note that we
                // only search within the surrounding since the plane
                // reconstruction might find a feature where there isn't one.
                const point estimatedPt(pt + patchAttraction[pointi]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                // Geometric feature edge. Check for different patchIDs
                bool hasSnapped = false;
                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        if (multiPatchPt.index() == 0)
                        {
                            // Region edge is also a geometric feature edge
                            nearInfo = findNearFeatureEdge
                            (
                                true,               // isRegionEdge
                                pp,
                                ppLocalPoints,
                                snapDist,
                                pointi,
                                estimatedPt,

                                edgeAttractors,
                                edgeConstraints,

                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;

                            // Debug: dump missed feature point
                            if (missedMP0Str && !nearInfo.second().hit())
                            {
                                missedMP0Str().writeLine(pt, estimatedPt);
                            }
                        }
                        else
                        {
                            // One of planes of feature contains multiple
                            // regions. We assume (contentious!) that the
                            // separation between
                            // the regions is not aligned with the geometric
                            // feature so is an additional constraint on the
                            // point -> is region-feature-point.
                            nearInfo = findNearFeaturePoint
                            (
                                true,           // isRegionPoint
                                pp,
                                ppLocalPoints,
                                snapDist,
                                pointi,
                                estimatedPt,

                                // Feature-point to pp point
                                pointAttractor,
                                pointConstraints,
                                // Feature-edge to pp point
                                edgeAttractors,
                                edgeConstraints,
                                // pp point to nearest feature
                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;

                            // More contentious: if we don't find
                            // a near feature point we will never find the
                            // attraction to a feature edge either since
                            // the edgeAttractors/edgeConstraints do not get
                            // filled and we're using reverse attraction
                            // Note that we're in multiRegionFeatureSnap which
                            // where findMultiPatchPoint can decide the
                            // wrong thing. So: if failed finding a near
                            // feature point try for a feature edge
                            if (!nearInfo.second().hit())
                            {
                                nearInfo = findNearFeatureEdge
                                (
                                    true,           // isRegionEdge
                                    pp,
                                    ppLocalPoints,
                                    snapDist,
                                    pointi,
                                    estimatedPt,

                                    // Feature-edge to pp point
                                    edgeAttractors,
                                    edgeConstraints,
                                    // pp point to nearest feature
                                    patchAttraction,
                                    patchConstraints
                                );
                            }

                            // Debug: dump missed feature point
                            if (missedMP1Str && !nearInfo.second().hit())
                            {
                                missedMP1Str().writeLine(pt, estimatedPt);
                            }
                        }
                    }
                }

                if (!hasSnapped)
                {
                    // Determine nearest point on feature edge. Store
                    // constraint
                    // (calculated from feature edge, alternative would be to
                    //  use constraint calculated from both surfaceNormals)
                    nearInfo = findNearFeatureEdge
                    (
                        false,      // isRegionPoint
                        pp,
                        ppLocalPoints,
                        snapDist,
                        pointi,
                        estimatedPt,

                        edgeAttractors,
                        edgeConstraints,

                        patchAttraction,
                        patchConstraints
                    );
                    hasSnapped = true;
                }

                // Dump to obj
                const pointIndexHit& info = nearInfo.second();
                if (info.hit())
                {
                    if
                    (
                        featurePointStr
                     && patchConstraints[pointi].first() == 3
                    )
                    {
                        featurePointStr().writeLine(pt, info.point());
                    }
                    else if
                    (
                        featureEdgeStr
                     && patchConstraints[pointi].first() == 2
                    )
                    {
                        featureEdgeStr().writeLine(pt, info.point());
                    }
                }
                else
                {
                    if (missedEdgeStr)
                    {
                        missedEdgeStr().writeLine(pt, estimatedPt);
                    }
                }
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                // Mark point on the nearest feature point.
                const point estimatedPt(pt + patchAttraction[pointi]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        // Multiple regions
                        nearInfo = findNearFeaturePoint
                        (
                            true,           // isRegionPoint
                            pp,
                            ppLocalPoints,
                            snapDist,
                            pointi,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                    else
                    {
                        nearInfo = findNearFeaturePoint
                        (
                            false,              // isRegionPoint
                            pp,
                            ppLocalPoints,
                            snapDist,
                            pointi,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // No multi-patch snapping
                    nearInfo = findNearFeaturePoint
                    (
                        false,              // isRegionPoint
                        pp,
                        ppLocalPoints,
                        snapDist,
                        pointi,
                        estimatedPt,

                        // Feature-point to pp point
                        pointAttractor,
                        pointConstraints,
                        // Feature-edge to pp point
                        edgeAttractors,
                        edgeConstraints,
                        // pp point to nearest feature
                        patchAttraction,
                        patchConstraints
                    );
                }

                const pointIndexHit& info = nearInfo.second();
                if (featurePointStr && info.hit())
                {
                    featurePointStr().writeLine(pt, info.point());
                }
            }
        }
    }
}


// Baffle handling
// ~~~~~~~~~~~~~~~
// Override pointAttractor, edgeAttractor, patchAttration etc. to
// implement 'baffle' handling.
// Baffle: the mesh pp point originates from a loose standing
// baffle.
// Sampling the surface with the surrounding face-centres only picks up
// a single triangle normal so above determineFeatures will not have
// detected anything. So explicitly pick up feature edges on the pp
// (after duplicating points & smoothing so will already have been
// expanded) and match these to the features.
void Foam::snappySnapDriver::determineBaffleFeatures
(
    const label iter,
    const bool baffleFeaturePoints,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementFeatures& features = meshRefiner_.features();

    // Calculate edge-faces
    List<List<point>> edgeFaceNormals(pp.nEdges());

    vectorField ppFaceNormals(pp.size());
    forAll(ppFaceNormals, facei)
    {
        ppFaceNormals[facei] = pp.localFaces()[facei].unitNormal(ppLocalPoints);
    }


    // Fill local data
    forAll(pp.edgeFaces(), edgei)
    {
        const labelList& eFaces = pp.edgeFaces()[edgei];
        List<point>& eFc = edgeFaceNormals[edgei];
        eFc.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label facei = eFaces[i];
            //eFc[i] = pp.faceNormals()[facei];
            eFc[i] = ppFaceNormals[facei];
        }
    }

    {
        // Precalculate mesh edges for pp.edges.
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        // Collect all coupled edges. Does not filter duplicates/order
        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            edgeFaceNormals,
            listPlusEqOp<point>(),
            List<point>(),
            mapDistribute::transform(),
            identityOp()  // No flipping
        );
    }

    // Detect baffle edges. Assume initial mesh will have 0,90 or 180
    // (baffle) degree angles so smoothing should make 0,90
    // to be less than 90. Choose reasonable value
    const scalar baffleFeatureCos = Foam::cos(degToRad(110.0));


    autoPtr<OBJstream> baffleEdgeStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        baffleEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "baffleEdge_" + name(iter) + ".obj"
            )
        );
        Info<< nl << "Dumping baffle-edges to "
            << baffleEdgeStr().name() << endl;
    }


    // Is edge on baffle
    bitSet isBaffleEdge(pp.nEdges());
    label nBaffleEdges = 0;
    // Is point on
    //  0 : baffle-edge (0)
    //  1 : baffle-feature-point (1)
    // -1 : rest
    labelList pointStatus(pp.nPoints(), -1);

    forAll(edgeFaceNormals, edgei)
    {
        const List<point>& efn = edgeFaceNormals[edgei];

        if (efn.size() == 2 && (efn[0]&efn[1]) < baffleFeatureCos)
        {
            isBaffleEdge.set(edgei);
            ++nBaffleEdges;
            const edge& e = pp.edges()[edgei];
            pointStatus[e[0]] = 0;
            pointStatus[e[1]] = 0;

            if (baffleEdgeStr)
            {
                baffleEdgeStr().write(e, ppLocalPoints);
            }
        }
    }

    reduce(nBaffleEdges, sumOp<label>());

    Info<< "Detected " << nBaffleEdges
        << " baffle edges out of "
        << returnReduce(pp.nEdges(), sumOp<label>())
        << " edges." << endl;


    //- Baffle edges will be too ragged to sensibly determine feature points
    //forAll(pp.pointEdges(), pointi)
    //{
    //    if
    //    (
    //        isFeaturePoint
    //        (
    //            featureCos,
    //            pp,
    //            isBaffleEdge,
    //            pointi
    //        )
    //    )
    //    {
    //        //Pout<< "Detected feature point:" << ppLocalPoints[pointi]
    //        //    << endl;
    //        //-TEMPORARILY DISABLED:
    //        //pointStatus[pointi] = 1;
    //    }
    //}


    label nBafflePoints = 0;
    forAll(pointStatus, pointi)
    {
        if (pointStatus[pointi] != -1)
        {
            nBafflePoints++;
        }
    }
    reduce(nBafflePoints, sumOp<label>());


    label nPointAttract = 0;
    label nEdgeAttract = 0;

    forAll(pointStatus, pointi)
    {
        const point& pt = ppLocalPoints[pointi];

        if (pointStatus[pointi] == 0)   // baffle edge
        {
            // 1: attract to near feature edge first

            Tuple2<label, pointIndexHit> nearInfo = findNearFeatureEdge
            (
                false,          // isRegionPoint?
                pp,
                ppLocalPoints,
                snapDist,
                pointi,
                pt,

                edgeAttractors,
                edgeConstraints,
                patchAttraction,
                patchConstraints
            );


            //- MEJ:
            // 2: optionally override with nearest feature point.
            //    On baffles we don't have enough normals to construct a feature
            //    point so assume all feature edges are close to feature points
            if (nearInfo.second().hit())
            {
                nEdgeAttract++;

                if (baffleFeaturePoints)
                {
                    nearInfo = findNearFeaturePoint
                    (
                        false,          // isRegionPoint,

                        pp,
                        ppLocalPoints,
                        snapDist,
                        pointi,
                        pt,             // estimatedPt,

                        // Feature-point to pp point
                        pointAttractor,
                        pointConstraints,
                        // Feature-edge to pp point
                        edgeAttractors,
                        edgeConstraints,
                        // pp point to nearest feature
                        patchAttraction,
                        patchConstraints
                    );

                    if (nearInfo.first() != -1)
                    {
                        nEdgeAttract--;
                        nPointAttract++;
                    }
                }
            }
        }
        else if (pointStatus[pointi] == 1)   // baffle point
        {
            labelList nearFeat;
            List<pointIndexHit> nearInfo;
            features.findNearestPoint
            (
                pointField(1, pt),
                scalarField(1, sqr(snapDist[pointi])),
                nearFeat,
                nearInfo
            );

            label feati = nearFeat[0];

            if (feati != -1)
            {
                nPointAttract++;

                label featPointi = nearInfo[0].index();
                const point& featPt = nearInfo[0].hitPoint();
                scalar distSqr = featPt.distSqr(pt);

                // Check if already attracted
                label oldPointi = pointAttractor[feati][featPointi];

                if
                (
                    oldPointi == -1
                 || (
                        distSqr
                      < featPt.distSqr(ppLocalPoints[oldPointi])
                    )
                )
                {
                    pointAttractor[feati][featPointi] = pointi;
                    pointConstraints[feati][featPointi].first() = 3;
                    pointConstraints[feati][featPointi].second() = Zero;

                    // Store for later use
                    patchAttraction[pointi] = featPt-pt;
                    patchConstraints[pointi] =
                        pointConstraints[feati][featPointi];

                    if (oldPointi != -1)
                    {
                        // The current point is closer so wins. Reset
                        // the old point to attract to nearest edge
                        // instead.
                        findNearFeatureEdge
                        (
                            false,              // isRegionPoint
                            pp,
                            ppLocalPoints,
                            snapDist,
                            oldPointi,
                            ppLocalPoints[oldPointi],

                            edgeAttractors,
                            edgeConstraints,
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // Make it fall through to check below
                    feati = -1;
                }
            }

            // Not found a feature point or another point is already
            // closer to that feature
            if (feati == -1)
            {
                //Pout<< "*** Falling back to finding nearest feature"
                //    << " edge"
                //    << " for baffle-feature-point " << pt
                //    << endl;

                Tuple2<label, pointIndexHit> nearInfo = findNearFeatureEdge
                (
                    false,                  // isRegionPoint
                    pp,
                    ppLocalPoints,
                    snapDist,
                    pointi,
                    pt,                     // starting point

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );

                if (nearInfo.first() != -1)
                {
                    nEdgeAttract++;
                }
            }
        }
    }

    reduce(nPointAttract, sumOp<label>());
    reduce(nEdgeAttract, sumOp<label>());

    Info<< "Baffle points     : " << nBafflePoints
        << " of which attracted to :" << nl
        << "    feature point : " << nPointAttract << nl
        << "    feature edge  : " << nEdgeAttract << nl
        << "    rest          : " << nBafflePoints-nPointAttract-nEdgeAttract
        << nl
        << endl;
}


void Foam::snappySnapDriver::reverseAttractMeshPoints
(
    const label iter,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,

    // Feature-point to pp point
    const List<labelList>& pointAttractor,
    const List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    const List<List<DynamicList<point>>>& edgeAttractors,
    const List<List<DynamicList<pointConstraint>>>& edgeConstraints,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reverse lookup : go through all edgeAttractors and find the
    // nearest point on pp

    // Get search domain and extend it a bit
    treeBoundBox bb(ppLocalPoints);
    {
        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb.inflate(rndGen, 1e-4, ROOTVSMALL);
    }

    // Collect candidate points for attraction
    DynamicList<label> attractPoints(pp.nPoints());
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        boolList isFeatureEdgeOrPoint(pp.nPoints(), false);
        label nFeats = 0;
        forAll(rawPatchConstraints, pointi)
        {
            if (rawPatchConstraints[pointi].first() >= 2)
            {
                isFeatureEdgeOrPoint[pointi] = true;
                nFeats++;
            }
        }

        Info<< "Initially selected " << returnReduce(nFeats, sumOp<label>())
            << " mesh points out of "
            << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;

        // Make sure is synchronised (note: check if constraint is already
        // synced in which case this is not needed here)
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            isFeatureEdgeOrPoint,
            orEqOp<bool>(),         // combine op
            false
        );

        for (label nGrow = 0; nGrow < 1; nGrow++)
        {
            boolList newIsFeatureEdgeOrPoint(isFeatureEdgeOrPoint);

            forAll(pp.localFaces(), facei)
            {
                const face& f = pp.localFaces()[facei];

                forAll(f, fp)
                {
                    if (isFeatureEdgeOrPoint[f[fp]])
                    {
                        // Mark all points on face
                        forAll(f, fp)
                        {
                            newIsFeatureEdgeOrPoint[f[fp]] = true;
                        }
                        break;
                    }
                }
            }

            isFeatureEdgeOrPoint = newIsFeatureEdgeOrPoint;

            syncTools::syncPointList
            (
                mesh,
                pp.meshPoints(),
                isFeatureEdgeOrPoint,
                orEqOp<bool>(),         // combine op
                false
            );
        }


        // Collect attractPoints
        forAll(isFeatureEdgeOrPoint, pointi)
        {
            if (isFeatureEdgeOrPoint[pointi])
            {
                attractPoints.append(pointi);
            }
        }

        Info<< "Selected "
            << returnReduce(attractPoints.size(), sumOp<label>())
            << " mesh points out of "
            << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;
    }


    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(ppLocalPoints, attractPoints),
        bb,                             // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    // Per mesh point the point on nearest feature edge.
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = Zero;
    patchConstraints.setSize(pp.nPoints());
    patchConstraints = pointConstraint();

    forAll(edgeAttractors, feati)
    {
        const List<DynamicList<point>>& edgeAttr = edgeAttractors[feati];
        const List<DynamicList<pointConstraint>>& edgeConstr =
            edgeConstraints[feati];

        forAll(edgeAttr, featEdgei)
        {
            const DynamicList<point>& attr = edgeAttr[featEdgei];
            forAll(attr, i)
            {
                // Find nearest pp point
                const point& featPt = attr[i];
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointi =
                        ppTree.shapes().objectIndex(nearInfo.index());

                    const point attraction
                    (
                        featPt
                      - ppTree.shapes()[nearInfo.index()]
                    );

                    // Check if this point is already being attracted. If so
                    // override it only if nearer.
                    if
                    (
                        patchConstraints[pointi].first() <= 1
                     || magSqr(attraction) < magSqr(patchAttraction[pointi])
                    )
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = edgeConstr[featEdgei][i];
                    }
                }
                else
                {
                    static label nWarn = 0;

                    if (nWarn < 100)
                    {
                        WarningInFunction
                            << "Did not find pp point near " << featPt
                            << endl;
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
        }
    }


    // Different procs might have different patchAttraction,patchConstraints
    // however these only contain geometric information, no topology
    // so as long as we synchronise after overriding with feature points
    // there is no problem, just possibly a small error.


    // Find nearest mesh point to feature point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (overrides attraction to feature edge)
    forAll(pointAttractor, feati)
    {
        const labelList& pointAttr = pointAttractor[feati];
        const List<pointConstraint>& pointConstr = pointConstraints[feati];

        forAll(pointAttr, featPointi)
        {
            if (pointAttr[featPointi] != -1)
            {
                const point& featPt = features[feati].points()[featPointi];

                // Find nearest pp point
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointi =
                        ppTree.shapes().objectIndex(nearInfo.index());

                    const point attraction
                    (
                        featPt
                      - ppTree.shapes()[nearInfo.index()]
                    );

                    // - already attracted to feature edge : point always wins
                    // - already attracted to feature point: nearest wins

                    if (patchConstraints[pointi].first() <= 1)
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = pointConstr[featPointi];
                    }
                    else if (patchConstraints[pointi].first() == 2)
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = pointConstr[featPointi];
                    }
                    else if (patchConstraints[pointi].first() == 3)
                    {
                        // Only if nearer
                        if
                        (
                            magSqr(attraction)
                          < magSqr(patchAttraction[pointi])
                        )
                        {
                            patchAttraction[pointi] = attraction;
                            patchConstraints[pointi] =
                                pointConstr[featPointi];
                        }
                    }
                }
            }
        }
    }
}


void Foam::snappySnapDriver::featureAttractionUsingFeatureEdges
(
    const label iter,
    const bool multiRegionFeatureSnap,
    const bool strictRegionFeatureSnap, // special >=3 patch points handling

    const bool detectBaffles,
    const bool baffleFeaturePoints,

    const bool releasePoints,
    const bool stringFeatures,
    const bool avoidDiagonal,

    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const vectorField& nearestNormal,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();
    const fvMesh& mesh = meshRefiner_.mesh();

    const bitSet isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            pp.meshPoints()
        )
    );


    // Collect ordered attractions on feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per feature, per feature-edge a list of attraction points and their
    // originating vertex.
    List<List<DynamicList<point>>> edgeAttractors(features.size());
    List<List<DynamicList<pointConstraint>>> edgeConstraints
    (
        features.size()
    );
    forAll(features, feati)
    {
        label nFeatEdges = features[feati].edges().size();
        edgeAttractors[feati].setSize(nFeatEdges);
        edgeConstraints[feati].setSize(nFeatEdges);
    }

    // Per feature, per feature-point the pp point that is attracted to it.
    // This list is only used to subset the feature-points that are actually
    // used.
    List<labelList> pointAttractor(features.size());
    List<List<pointConstraint>> pointConstraints(features.size());
    forAll(features, feati)
    {
        label nFeatPoints = features[feati].points().size();
        pointAttractor[feati].setSize(nFeatPoints, -1);
        pointConstraints[feati].setSize(nFeatPoints);
    }

    // Reverse: from pp point to nearest feature
    vectorField rawPatchAttraction(pp.nPoints(), Zero);
    List<pointConstraint> rawPatchConstraints(pp.nPoints());

    determineFeatures
    (
        iter,
        featureCos,
        multiRegionFeatureSnap,
        strictRegionFeatureSnap, // extra logic for points on >=3 patches 

        pp,
        ppLocalPoints,
        snapDist,               // per point max distance and nearest surface
        nearestDisp,

        pointFaceSurfNormals,   // per face nearest surface
        pointFaceDisp,
        pointFaceCentres,
        pointFacePatchID,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        rawPatchAttraction,
        rawPatchConstraints
    );

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "Raw geometric feature analysis : ";
        writeStats(pp, isPatchMasterPoint, rawPatchConstraints);
    }

    // Baffle handling
    // ~~~~~~~~~~~~~~~
    // Override pointAttractor, edgeAttractor, rawPatchAttration etc. to
    // implement 'baffle' handling.
    // Baffle: the mesh pp point originates from a loose standing
    // baffle.
    // Sampling the surface with the surrounding face-centres only picks up
    // a single triangle normal so above determineFeatures will not have
    // detected anything. So explicitly pick up feature edges on the pp
    // (after duplicating points & smoothing so will already have been
    // expanded) and match these to the features.
    if (detectBaffles)
    {
        determineBaffleFeatures
        (
            iter,
            baffleFeaturePoints,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,

            // Feature-point to pp point
            pointAttractor,
            pointConstraints,
            // Feature-edge to pp point
            edgeAttractors,
            edgeConstraints,
            // pp point to nearest feature
            rawPatchAttraction,
            rawPatchConstraints
        );
    }

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "After baffle feature analysis : ";
        writeStats(pp, isPatchMasterPoint, rawPatchConstraints);
    }


    // Reverse lookup: Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~----------------~~~~~~~~~~~~~~~~~~~~~~~
    // go through all edgeAttractors and find the nearest point on pp

    reverseAttractMeshPoints
    (
        iter,

        pp,
        ppLocalPoints,
        snapDist,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,

        // Estimated feature point
        rawPatchAttraction,
        rawPatchConstraints,

        // pp point to nearest feature
        patchAttraction,
        patchConstraints
    );

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "Reverse attract feature analysis : ";
        writeStats(pp, isPatchMasterPoint, patchConstraints);
    }

    // Dump
    if (debug&meshRefinement::ATTRACTION)
    {
        OBJstream featureEdgeStr
        (
            meshRefiner_.mesh().time().path()
          / "edgeAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-edge attraction to "
            << featureEdgeStr.name() << endl;

        OBJstream featurePointStr
        (
            meshRefiner_.mesh().time().path()
          / "pointAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-point attraction to "
            << featurePointStr.name() << endl;

        forAll(patchConstraints, pointi)
        {
            const point& pt = ppLocalPoints[pointi];
            const vector& attr = patchAttraction[pointi];

            if (patchConstraints[pointi].first() == 2)
            {
                featureEdgeStr.writeLine(pt, pt+attr);
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                featurePointStr.writeLine(pt, pt+attr);
            }
        }
    }


    //MEJ: any faces that have multi-patch points only keep the
    //     multi-patch
    //     points. The other points on the face will be dragged along
    //     (hopefully)
    if (releasePoints)
    {
        releasePointsNextToMultiPatch
        (
            iter,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,

            pointFaceCentres,
            pointFacePatchID,

            rawPatchAttraction,
            rawPatchConstraints,

            patchAttraction,
            patchConstraints
        );
    }


    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)
    if (stringFeatures)
    {
        stringFeatureEdges
        (
            iter,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,

            rawPatchAttraction,
            rawPatchConstraints,

            patchAttraction,
            patchConstraints
        );
    }


    // Avoid diagonal attraction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Attract one of the non-diagonal points.
    if (avoidDiagonal)
    {
        avoidDiagonalAttraction
        (
            iter,
            featureCos,
            pp,
            ppLocalPoints,
            patchAttraction,
            patchConstraints
        );
    }


    if (debug&meshRefinement::ATTRACTION)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            ppLocalPoints,
            ppLocalPoints + patchAttraction
        );
    }
}


// Correct for squeezing of face
void Foam::snappySnapDriver::preventFaceSqueeze
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const scalarField& snapDist,
    const vectorField& nearestAttraction,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> strPtr;
    if (debug&meshRefinement::ATTRACTION)
    {
        strPtr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "faceSqueeze_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping faceSqueeze corrections to "
            << strPtr().name() << endl;
    }

    pointField points;
    face singleF;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        if (f.size() != points.size())
        {
            points.setSize(f.size());
            singleF.setSize(f.size());
            for (label i = 0; i < f.size(); i++)
            {
                singleF[i] = i;
            }
        }
        label nConstraints = 0;
        forAll(f, fp)
        {
            label pointi = f[fp];
            const point& pt = ppLocalPoints[pointi];

            if (patchConstraints[pointi].first() > 1)
            {
                points[fp] = pt + patchAttraction[pointi];
                nConstraints++;
            }
            else
            {
                points[fp] = pt;
            }
        }

        if (nConstraints == f.size())
        {
            if (f.size() == 3)
            {
                // Triangle: knock out attraction altogether

                // For now keep the points on the longest edge
                label maxFp = -1;
                scalar maxS = -1;
                forAll(f, fp)
                {
                    const point& pt = ppLocalPoints[f[fp]];
                    const point& nextPt = ppLocalPoints[f.nextLabel(fp)];

                    scalar s = pt.distSqr(nextPt);
                    if (s > maxS)
                    {
                        maxS = s;
                        maxFp = fp;
                    }
                }
                if (maxFp != -1)
                {
                    label pointi = f.prevLabel(maxFp);

                    // Reset attraction on pointi to nearest

                    const point& pt = ppLocalPoints[pointi];

                    //Pout<< "** on triangle " << pp.faceCentres()[facei]
                    //    << " knocking out attraction to " << pointi
                    //    << " at:" << pt
                    //    << endl;

                    patchAttraction[pointi] = nearestAttraction[pointi];

                    if (strPtr)
                    {
                        strPtr().writeLine(pt, pt+patchAttraction[pointi]);
                    }
                }
            }
            else
            {
                scalar oldArea = f.mag(ppLocalPoints);
                scalar newArea = singleF.mag(points);
                if (newArea < 0.1*oldArea)
                {
                    // For now remove the point with largest distance
                    label maxFp = -1;
                    scalar maxS = -1;
                    forAll(f, fp)
                    {
                        scalar s = magSqr(patchAttraction[f[fp]]);
                        if (s > maxS)
                        {
                            maxS = s;
                            maxFp = fp;
                        }
                    }
                    if (maxFp != -1)
                    {
                        label pointi = f[maxFp];
                        // Lower attraction on pointi
                        patchAttraction[pointi] *= 0.5;
                    }
                }
            }
        }
    }
}


Foam::vectorField Foam::snappySnapDriver::calcNearestSurfaceFeature
(
    const snapParameters& snapParams,
    const bool alignMeshEdges,
    const bool strictRegionFeatureSnap, // use patches
    const label iter,
    const scalar featureCos,
    const scalar featureAttract,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const vectorField& nearestNormal,
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits,
    DynamicList<labelPair>& splitPatches
) const
{
    if (dryRun_)
    {
        return nearestDisp;
    }


    const Switch implicitFeatureAttraction = snapParams.implicitFeatureSnap();
    const Switch explicitFeatureAttraction = snapParams.explicitFeatureSnap();
    const Switch multiRegionFeatureSnap = snapParams.multiRegionFeatureSnap();

    Info<< "Overriding displacement on features :" << nl
        << "   implicit features    : " << implicitFeatureAttraction << nl
        << "   explicit features    : " << explicitFeatureAttraction << nl
        << "   multi-patch features : " << multiRegionFeatureSnap << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();


    //const bitSet isMasterPoint(syncTools::getMasterPoints(mesh));
    const bitSet isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            pp.meshPoints()
        )
    );

    // Per point, per surrounding face:
    // - faceSurfaceNormal
    // - faceDisp
    // - faceCentres
    List<List<point>> pointFaceSurfNormals;
    List<List<point>> pointFaceDisp;
    List<List<point>> pointFaceCentres;
    List<labelList> pointFacePatchID;

    {
        // Calculate attraction distance per face (from the attraction distance
        // per point)
        scalarField faceSnapDist(pp.size(), -GREAT);
        forAll(pp.localFaces(), facei)
        {
            const face& f = pp.localFaces()[facei];
            forAll(f, fp)
            {
                faceSnapDist[facei] = max(faceSnapDist[facei], snapDist[f[fp]]);
            }
        }


        // Displacement and orientation per pp face
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // vector from point on surface back to face centre
        vectorField faceDisp(pp.size(), Zero);
        // normal of surface at point on surface
        vectorField faceSurfaceNormal(pp.size(), Zero);
        labelList faceSurfaceGlobalRegion(pp.size(), -1);
        //vectorField faceRotation(pp.size(), Zero);

        calcNearestFace
        (
            iter,
            pp,
            ppLocalPoints,
            faceSnapDist,
            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion
            //faceRotation
        );


        // Collect (possibly remote) per point data of all surrounding faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // - faceSurfaceNormal
        // - faceDisp
        // - faceCentres
        calcNearestFacePointProperties
        (
            iter,
            pp,
            ppLocalPoints,
            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID
        );
    }


    // Start off with nearest point on surface
    vectorField patchDisp = nearestDisp;


    // Main calculation
    // ~~~~~~~~~~~~~~~~
    // This is the main intelligence which calculates per point the vector to
    // attract it to the nearest surface. There are lots of possibilities
    // here.

    // Nearest feature
    patchAttraction.setSize(ppLocalPoints.size());
    patchAttraction = Zero;
    // Constraints at feature
    patchConstraints.setSize(ppLocalPoints.size());
    patchConstraints = pointConstraint();

    if (implicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. Reconstruct a feature edge/point if possible and snap to
        // it.
        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            ppLocalPoints,
            snapDist,
            nearestDisp,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (explicitFeatureAttraction)
    {
        // Only do fancy stuff if alignMeshEdges
        bool releasePoints = false;
        bool stringFeatures = false;
        bool avoidDiagonal = false;
        if (alignMeshEdges)
        {
            releasePoints = snapParams.releasePoints();
            stringFeatures = snapParams.stringFeatures();
            avoidDiagonal = snapParams.avoidDiagonal();
        }


        // Sample faces around each point and see if nearest surface normal
        // differs. For those find the nearest real feature edge/point and
        // store the correspondence. Then loop over feature edge/point
        // and attract those nearest mesh point. (the first phase just is
        // a subsetting of candidate points, the second makes sure that only
        // one mesh point gets attracted per feature)
        featureAttractionUsingFeatureEdges
        (
            iter,
            multiRegionFeatureSnap,
            strictRegionFeatureSnap,    // special logic for >=3 patches points?

            snapParams.detectBaffles(),
            snapParams.baffleFeaturePoints(),   // all points on baffle edges
                                                // are attracted to feature pts

            releasePoints,
            stringFeatures,
            avoidDiagonal,

            featureCos,

            pp,
            ppLocalPoints,
            snapDist,
            nearestDisp,
            nearestNormal,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (!alignMeshEdges)
    {
        const scalar concaveCos = Foam::cos
        (
            degToRad(snapParams.concaveAngle())
        );
        const scalar minAreaRatio = snapParams.minAreaRatio();

        Info<< "Introducing face splits to avoid rotating"
            << " mesh edges. Splitting faces when" << nl
            << indent << "- angle not concave by more than "
            << snapParams.concaveAngle() << " degrees" << nl
            << indent << "- resulting triangles of similar area "
            << " (ratio within " << minAreaRatio << ")" << nl
            << endl;

        splitDiagonals
        (
            featureCos,
            concaveCos,
            minAreaRatio,
            pp,
            ppLocalPoints,

            nearestDisp,
            nearestNormal,
            pointFacePatchID,   // per point the connected patches. Is used
                                // to re-patch newly split faces
            patchAttraction,
            patchConstraints,

            splitFaces,
            splits,
            splitPatches
        );

        if (debug)
        {
            Info<< "Diagonal attraction feature correction : ";
            writeStats(pp, isPatchMasterPoint, patchConstraints);
        }
    }


    preventFaceSqueeze
    (
        iter,
        featureCos,

        pp,
        ppLocalPoints,
        snapDist,
        nearestDisp,

        patchAttraction,
        patchConstraints
    );

    {
        vector avgPatchDisp = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchDisp
        );
        vector avgPatchAttr = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchAttraction
        );

        Info<< "Attraction:" << endl
            << "    linear   : max:" << gMaxMagSqr(patchDisp)
            << " avg:" << avgPatchDisp << endl
            << "    feature  : max:" << gMaxMagSqr(patchAttraction)
            << " avg:" << avgPatchAttr << endl;
    }

    // So now we have:
    // - patchDisp          : point movement to go to nearest point on surface
    //                       (either direct or through interpolation of
    //                        face nearest)
    // - patchAttraction    : direct attraction to features
    // - patchConstraints   : type of features

    // Use any combination of patchDisp and direct feature attraction.


    // Mix with direct feature attraction
    forAll(patchConstraints, pointi)
    {
        if (patchConstraints[pointi].first() > 1)
        {
            patchDisp[pointi] =
                (1.0-featureAttract)*patchDisp[pointi]
              + featureAttract*patchAttraction[pointi];
        }
    }



    // Count
    {
        Info<< "Feature analysis : ";
        writeStats(pp, isPatchMasterPoint, patchConstraints);
    }


    // Now we have the displacement per patch point to move onto the surface
    // Split into tangential and normal direction.
    // - start off with all non-constrained points following the constrained
    //   ones since point normals not relevant.
    // - finish with only tangential component smoothed.
    // Note: tangential is most
    // likely to come purely from face-centre snapping, not face rotation.
    // Note: could use the constraints here (constraintTransformation())
    //       but this is not necessarily accurate and we're smoothing to
    //       get out of problems.

    if (featureAttract < 1-0.001)
    {
        //const bitSet isMasterEdge(syncTools::getMasterEdges(mesh));
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        const bitSet isPatchMasterEdge
        (
            meshRefinement::getMasterEdges
            (
                mesh,
                meshEdges
            )
        );

        const vectorField pointNormals
        (
            PatchTools::pointNormals
            (
                mesh,
                pp
            )
        );

        // 1. Smoothed all displacement
        vectorField smoothedPatchDisp = patchDisp;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            smoothedPatchDisp
        );


        // 2. Smoothed tangential component
        vectorField tangPatchDisp = patchDisp;
        tangPatchDisp -= (pointNormals & patchDisp) * pointNormals;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            tangPatchDisp
        );

        // Re-add normal component
        tangPatchDisp += (pointNormals & patchDisp) * pointNormals;

        if (debug&meshRefinement::ATTRACTION)
        {
            dumpMove
            (
                mesh.time().path()
              / "tangPatchDispConstrained_" + name(iter) + ".obj",
                ppLocalPoints,
                ppLocalPoints + tangPatchDisp
            );
        }

        patchDisp =
             (1.0-featureAttract)*smoothedPatchDisp
           + featureAttract*tangPatchDisp;
    }


    const scalar relax = featureAttract;
    patchDisp *= relax;


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


// ************************************************************************* //
