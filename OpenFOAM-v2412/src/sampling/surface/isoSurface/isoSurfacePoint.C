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

\*---------------------------------------------------------------------------*/

#include "isoSurfacePoint.H"
#include "mergePoints.H"
#include "slicedVolFields.H"
#include "volFields.H"
#include "triSurfaceTools.H"
#include "triSurface.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isoSurfaceBaseMethods.H"
defineIsoSurfaceInterpolateMethods(Foam::isoSurfacePoint);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfacePoint, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Helper class for slicing triangles
struct storeOp
{
    DynamicList<triPoints>& list;

    storeOp(DynamicList<triPoints>& tris)
    :
        list(tris)
    {}

    void operator()(const triPoints& tri)
    {
        list.append(tri);
    }

    void operator()(triPoints&& tri)
    {
        list.append(std::move(tri));
    }
};


// Is patch co-located (i.e. non-separated) coupled patch?
static inline bool collocatedPatch(const polyPatch& pp)
{
    const auto* cpp = isA<coupledPolyPatch>(pp);
    return cpp && cpp->parallel() && !cpp->separated();
}


// Collocated patch, with extra checks
#if 0
static bool isCollocatedPatch(const coupledPolyPatch& pp)
{
    if (isA<processorPolyPatch>(pp) || isA<cyclicPolyPatch>(pp))
    {
        return collocatedPatch(pp);
    }

    FatalErrorInFunction
        << "Unhandled coupledPolyPatch type " << pp.type() << nl
        << abort(FatalError);

    return false;
}
#endif

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::bitSet Foam::isoSurfacePoint::collocatedFaces
(
    const coupledPolyPatch& pp
)
{
    // Initialise to false
    bitSet collocated(pp.size());

    if (isA<processorPolyPatch>(pp) || isA<cyclicPolyPatch>(pp))
    {
        if (collocatedPatch(pp))
        {
            // All on
            collocated = true;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unhandled coupledPolyPatch type " << pp.type()
            << abort(FatalError);
    }
    return collocated;
}


void Foam::isoSurfacePoint::syncUnseparatedPoints
(
    pointField& pointValues,
    const point& nullValue
) const
{
    // Until syncPointList handles separated coupled patches with multiple
    // transforms do our own synchronisation of non-separated patches only
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        const labelList& procPatches = mesh_.globalData().processorPatches();

        // Send
        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && collocatedPatch(pp))
            {
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                pointField patchInfo(meshPts.size());

                forAll(nbrPts, pointi)
                {
                    const label nbrPointi = nbrPts[pointi];
                    patchInfo[nbrPointi] = pointValues[meshPts[pointi]];
                }

                OPstream toNbr
                (
                    UPstream::commsTypes::buffered,
                    procPatch.neighbProcNo()
                );
                toNbr << patchInfo;
            }
        }

        // Receive and combine.
        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && collocatedPatch(pp))
            {
                // We do not know the number of points on the other side
                // so cannot use UIPstream::read

                pointField nbrPatchInfo;
                IPstream::recv(nbrPatchInfo, procPatch.neighbProcNo());

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointi)
                {
                    const label meshPointi = meshPts[pointi];
                    minEqOp<point>()
                    (
                        pointValues[meshPointi],
                        nbrPatchInfo[pointi]
                    );
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner() && collocatedPatch(*cpp))
        {
            // Owner does all.

            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();
            const labelList& nbrMeshPoints = nbrPatch.meshPoints();

            pointField half0Values(coupledPoints.size());
            pointField half1Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                half0Values[i] = pointValues[meshPts[e[0]]];
                half1Values[i] = pointValues[nbrMeshPoints[e[1]]];
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                const label p0 = meshPts[e[0]];
                const label p1 = nbrMeshPoints[e[1]];

                minEqOp<point>()(pointValues[p0], half1Values[i]);
                minEqOp<point>()(pointValues[p1], half0Values[i]);
            }
        }
    }

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh_.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            const label meshPointi = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointi];
        }

        // Globally consistent
        Pstream::listCombineReduce(sharedPts, minEqOp<point>());

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            const label meshPointi = pd.sharedPointLabels()[i];
            pointValues[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


Foam::scalar Foam::isoSurfacePoint::isoFraction
(
    const scalar s0,
    const scalar s1
) const
{
    const scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        return (iso_-s0)/d;
    }

    return -1.0;
}


bool Foam::isoSurfacePoint::isEdgeOfFaceCut
(
    const scalarField& pVals,
    const face& f,
    const bool ownLower,
    const bool neiLower
) const
{
    forAll(f, fp)
    {
        const bool fpLower = (pVals[f[fp]] < iso_);

        if
        (
            fpLower != ownLower
         || fpLower != neiLower
         || fpLower != (pVals[f[f.fcIndex(fp)]] < iso_)
        )
        {
            return true;
        }
    }

    return false;
}


void Foam::isoSurfacePoint::getNeighbour
(
    const labelList& boundaryRegion,
    const volVectorField& meshC,
    const volScalarField& cVals,
    const label celli,
    const label facei,
    scalar& nbrValue,
    point& nbrPoint
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    if (mesh_.isInternalFace(facei))
    {
        const label nbr = (own[facei] == celli ? nei[facei] : own[facei]);
        nbrValue = cVals[nbr];
        nbrPoint = meshC[nbr];
    }
    else
    {
        const label bFacei = facei-mesh_.nInternalFaces();
        const label patchi = boundaryRegion[bFacei];
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const label patchFacei = facei-pp.start();

        nbrValue = cVals.boundaryField()[patchi][patchFacei];
        nbrPoint = meshC.boundaryField()[patchi][patchFacei];
    }
}


void Foam::isoSurfacePoint::calcCutTypes
(
    const labelList& boundaryRegion,
    const volVectorField& meshC,
    const volScalarField& cVals,
    const scalarField& pVals
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    faceCutType_.resize(mesh_.nFaces());
    faceCutType_ = cutType::NOTCUT;

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        const scalar ownValue = cVals[own[facei]];
        const bool ownLower = (ownValue < iso_);

        scalar nbrValue;
        point nbrPoint;
        getNeighbour
        (
            boundaryRegion,
            meshC,
            cVals,
            own[facei],
            facei,
            nbrValue,
            nbrPoint
        );

        const bool neiLower = (nbrValue < iso_);

        if (ownLower != neiLower)
        {
            faceCutType_[facei] = cutType::CUT;
        }
        else
        {
            // Is mesh edge cut?
            // - by looping over all the edges of the face.
            const face& f = mesh_.faces()[facei];

            if (isEdgeOfFaceCut(pVals, f, ownLower, neiLower))
            {
                faceCutType_[facei] = cutType::CUT;
            }
        }
    }

    for (const polyPatch& pp : patches)
    {
        for (const label facei : pp.range())
        {
            const scalar ownValue = cVals[own[facei]];
            const bool ownLower = (ownValue < iso_);

            scalar nbrValue;
            point nbrPoint;
            getNeighbour
            (
                boundaryRegion,
                meshC,
                cVals,
                own[facei],
                facei,
                nbrValue,
                nbrPoint
            );

            const bool neiLower = (nbrValue < iso_);

            if (ownLower != neiLower)
            {
                faceCutType_[facei] = cutType::CUT;
            }
            else
            {
                // Is mesh edge cut?
                // - by looping over all the edges of the face.
                const face& f = mesh_.faces()[facei];

                if (isEdgeOfFaceCut(pVals, f, ownLower, neiLower))
                {
                    faceCutType_[facei] = cutType::CUT;
                }
            }
        }
    }

    nCutCells_ = 0;
    cellCutType_.resize(mesh_.nCells());
    cellCutType_ = cutType::NOTCUT;


    // Propagate internal face cuts into the cells.

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (faceCutType_[facei] == cutType::NOTCUT)
        {
            continue;
        }

        if (cellCutType_[own[facei]] == cutType::NOTCUT)
        {
            cellCutType_[own[facei]] = cutType::CUT;
            ++nCutCells_;
        }
        if (cellCutType_[nei[facei]] == cutType::NOTCUT)
        {
            cellCutType_[nei[facei]] = cutType::CUT;
            ++nCutCells_;
        }
    }


    // Propagate boundary face cuts into the cells.

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (faceCutType_[facei] == cutType::NOTCUT)
        {
            continue;
        }

        if (cellCutType_[own[facei]] == cutType::NOTCUT)
        {
            cellCutType_[own[facei]] = cutType::CUT;
            ++nCutCells_;
        }
    }

    if (debug)
    {
        Pout<< "isoSurfacePoint : candidate cut cells "
            << nCutCells_ << " / " << mesh_.nCells() << endl;
    }
}


Foam::point Foam::isoSurfacePoint::calcCentre(const triSurface& s)
{
    // Calculate centre of surface.

    vector sum = Zero;

    forAll(s, i)
    {
        sum += s[i].centre(s.points());
    }
    return sum/s.size();
}


void Foam::isoSurfacePoint::calcSnappedCc
(
    const labelList& boundaryRegion,
    const volVectorField& meshC,
    const volScalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedCc
) const
{
    const pointField& pts = mesh_.points();
    const pointField& cc = mesh_.cellCentres();

    snappedCc.setSize(mesh_.nCells());
    snappedCc = -1;

    // Work arrays
    DynamicList<point, 64> localTriPoints(64);

    forAll(mesh_.cells(), celli)
    {
        if (cellCutType_[celli] == cutType::CUT)
        {
            const scalar cVal = cVals[celli];

            localTriPoints.clear();
            label nOther = 0;
            point otherPointSum = Zero;

            // Create points for all intersections close to cell centre
            // (i.e. from pyramid edges)

            for (const label facei : mesh_.cells()[celli])
            {
                const face& f = mesh_.faces()[facei];

                scalar nbrValue;
                point nbrPoint;
                getNeighbour
                (
                    boundaryRegion,
                    meshC,
                    cVals,
                    celli,
                    facei,
                    nbrValue,
                    nbrPoint
                );

                // Calculate intersection points of edges to cell centre
                FixedList<scalar, 3> s;
                FixedList<point, 3> pt;

                // From cc to neighbour cc.
                s[2] = isoFraction(cVal, nbrValue);
                pt[2] = (1.0-s[2])*cc[celli] + s[2]*nbrPoint;

                forAll(f, fp)
                {
                    // From cc to fp
                    label p0 = f[fp];
                    s[0] = isoFraction(cVal, pVals[p0]);
                    pt[0] = (1.0-s[0])*cc[celli] + s[0]*pts[p0];

                    // From cc to fp+1
                    label p1 = f[f.fcIndex(fp)];
                    s[1] = isoFraction(cVal, pVals[p1]);
                    pt[1] = (1.0-s[1])*cc[celli] + s[1]*pts[p1];

                    if
                    (
                        (s[0] >= 0.0 && s[0] <= 0.5)
                     && (s[1] >= 0.0 && s[1] <= 0.5)
                     && (s[2] >= 0.0 && s[2] <= 0.5)
                    )
                    {
                        localTriPoints.append(pt[0]);
                        localTriPoints.append(pt[1]);
                        localTriPoints.append(pt[2]);
                    }
                    else
                    {
                        // Get average of all other points
                        forAll(s, i)
                        {
                            if (s[i] >= 0.0 && s[i] <= 0.5)
                            {
                                otherPointSum += pt[i];
                                ++nOther;
                            }
                        }
                    }
                }
            }

            if (localTriPoints.size() == 0)
            {
                // No complete triangles. Get average of all intersection
                // points.
                if (nOther > 0)
                {
                    snappedCc[celli] = snappedPoints.size();
                    snappedPoints.append(otherPointSum/nOther);

                    //Pout<< "    point:" << pointi
                    //    << " replacing coord:" << mesh_.points()[pointi]
                    //    << " by average:" << collapsedPoint[pointi] << endl;
                }
            }
            else if (localTriPoints.size() == 3)
            {
                // Single triangle. No need for any analysis. Average points.
                snappedCc[celli] = snappedPoints.size();
                snappedPoints.append(sum(localTriPoints)/3);
                localTriPoints.clear();

                //Pout<< "    point:" << pointi
                //    << " replacing coord:" << mesh_.points()[pointi]
                //    << " by average:" << collapsedPoint[pointi] << endl;
            }
            else
            {
                // Convert points into triSurface.

                // Merge points and compact out non-valid triangles
                labelList triMap;               // merged to unmerged triangle
                labelList triPointReverseMap;   // unmerged to merged point
                triSurface surf
                (
                    stitchTriPoints
                    (
                        false,              // do not check for duplicate tris
                        localTriPoints,
                        triPointReverseMap,
                        triMap
                    )
                );

                labelList faceZone;
                label nZones = surf.markZones
                (
                    boolList(surf.nEdges(), false),
                    faceZone
                );

                if (nZones == 1)
                {
                    snappedCc[celli] = snappedPoints.size();
                    snappedPoints.append(calcCentre(surf));
                    //Pout<< "    point:" << pointi << " nZones:" << nZones
                    //    << " replacing coord:" << mesh_.points()[pointi]
                    //    << " by average:" << collapsedPoint[pointi] << endl;
                }
            }
        }
    }
}


void Foam::isoSurfacePoint::calcSnappedPoint
(
    const bitSet& isBoundaryPoint,
    const labelList& boundaryRegion,
    const volVectorField& meshC,
    const volScalarField& cVals,
    const scalarField& pVals,

    DynamicList<point>& snappedPoints,
    labelList& snappedPoint
) const
{
    const pointField& pts = mesh_.points();
    const pointField& cc = mesh_.cellCentres();

    pointField collapsedPoint(mesh_.nPoints(), point::max);

    // Work arrays
    DynamicList<point, 64> localTriPoints(100);

    forAll(mesh_.pointFaces(), pointi)
    {
        if (isBoundaryPoint.test(pointi))
        {
            continue;
        }

        const labelList& pFaces = mesh_.pointFaces()[pointi];

        bool anyCut = false;

        for (const label facei : pFaces)
        {
            if (faceCutType_[facei] == cutType::CUT)
            {
                anyCut = true;
                break;
            }
        }

        if (!anyCut)
        {
            continue;
        }


        localTriPoints.clear();
        label nOther = 0;
        point otherPointSum = Zero;

        for (const label facei : pFaces)
        {
            // Create points for all intersections close to point
            // (i.e. from pyramid edges)
            const face& f = mesh_.faces()[facei];
            const label own = mesh_.faceOwner()[facei];

            // Get neighbour value
            scalar nbrValue;
            point nbrPoint;
            getNeighbour
            (
                boundaryRegion,
                meshC,
                cVals,
                own,
                facei,
                nbrValue,
                nbrPoint
            );

            // Calculate intersection points of edges emanating from point
            FixedList<scalar, 4> s;
            FixedList<point, 4> pt;

            label fp = f.find(pointi);
            s[0] = isoFraction(pVals[pointi], cVals[own]);
            pt[0] = (1.0-s[0])*pts[pointi] + s[0]*cc[own];

            s[1] = isoFraction(pVals[pointi], nbrValue);
            pt[1] = (1.0-s[1])*pts[pointi] + s[1]*nbrPoint;

            label nextPointi = f[f.fcIndex(fp)];
            s[2] = isoFraction(pVals[pointi], pVals[nextPointi]);
            pt[2] = (1.0-s[2])*pts[pointi] + s[2]*pts[nextPointi];

            label prevPointi = f[f.rcIndex(fp)];
            s[3] = isoFraction(pVals[pointi], pVals[prevPointi]);
            pt[3] = (1.0-s[3])*pts[pointi] + s[3]*pts[prevPointi];

            if
            (
                (s[0] >= 0.0 && s[0] <= 0.5)
             && (s[1] >= 0.0 && s[1] <= 0.5)
             && (s[2] >= 0.0 && s[2] <= 0.5)
            )
            {
                localTriPoints.append(pt[0]);
                localTriPoints.append(pt[1]);
                localTriPoints.append(pt[2]);
            }
            if
            (
                (s[0] >= 0.0 && s[0] <= 0.5)
             && (s[1] >= 0.0 && s[1] <= 0.5)
             && (s[3] >= 0.0 && s[3] <= 0.5)
            )
            {
                localTriPoints.append(pt[3]);
                localTriPoints.append(pt[0]);
                localTriPoints.append(pt[1]);
            }

            // Get average of points as fallback
            forAll(s, i)
            {
                if (s[i] >= 0.0 && s[i] <= 0.5)
                {
                    otherPointSum += pt[i];
                    ++nOther;
                }
            }
        }

        if (localTriPoints.size() == 0)
        {
            // No complete triangles. Get average of all intersection
            // points.
            if (nOther > 0)
            {
                collapsedPoint[pointi] = otherPointSum/nOther;
            }
        }
        else if (localTriPoints.size() == 3)
        {
            // Single triangle. No need for any analysis. Average points.
            pointField points;
            points.transfer(localTriPoints);
            collapsedPoint[pointi] = sum(points)/points.size();
        }
        else
        {
            // Convert points into triSurface.

            // Merge points and compact out non-valid triangles
            labelList triMap;               // merged to unmerged triangle
            labelList triPointReverseMap;   // unmerged to merged point
            triSurface surf
            (
                stitchTriPoints
                (
                    false,                  // do not check for duplicate tris
                    localTriPoints,
                    triPointReverseMap,
                    triMap
                )
            );

            labelList faceZone;
            label nZones = surf.markZones
            (
                boolList(surf.nEdges(), false),
                faceZone
            );

            if (nZones == 1)
            {
                collapsedPoint[pointi] = calcCentre(surf);
            }
        }
    }


    // Synchronise snap location
    syncUnseparatedPoints(collapsedPoint, point::max);


    snappedPoint.setSize(mesh_.nPoints());
    snappedPoint = -1;

    forAll(collapsedPoint, pointi)
    {
        if (collapsedPoint[pointi] != point::max)
        {
            snappedPoint[pointi] = snappedPoints.size();
            snappedPoints.append(collapsedPoint[pointi]);
        }
    }
}


Foam::triSurface Foam::isoSurfacePoint::stitchTriPoints
(
    const bool checkDuplicates,
    const List<point>& triPoints,
    labelList& triPointReverseMap,  // unmerged to merged point
    labelList& triMap               // merged to unmerged triangle
) const
{
    label nTris = triPoints.size()/3;

    if ((triPoints.size() % 3) != 0)
    {
        FatalErrorInFunction
            << "Problem: number of points " << triPoints.size()
            << " not a multiple of 3." << abort(FatalError);
    }

    pointField newPoints;
    mergePoints
    (
        triPoints,
        mergeDistance_,
        false,
        triPointReverseMap,
        newPoints
    );

    // Check that enough merged.
    if (debug)
    {
        Pout<< "isoSurfacePoint : merged from " << triPointReverseMap.size()
            << " down to " << newPoints.size() << " unique points." << endl;

        pointField newNewPoints;
        labelList oldToNew;
        bool hasMerged = mergePoints
        (
            newPoints,
            mergeDistance_,
            true,
            oldToNew,
            newNewPoints
        );

        if (hasMerged)
        {
            FatalErrorInFunction
                << "Merged points contain duplicates"
                << " when merging with distance " << mergeDistance_ << endl
                << "merged:" << oldToNew.size() << " re-merged:"
                << newNewPoints.size()
                << abort(FatalError);
        }
    }


    List<labelledTri> tris;
    {
        DynamicList<labelledTri> dynTris(nTris);
        label rawPointi = 0;
        DynamicList<label> newToOldTri(nTris);

        for (label oldTriI = 0; oldTriI < nTris; oldTriI++)
        {
            labelledTri tri
            (
                triPointReverseMap[rawPointi],
                triPointReverseMap[rawPointi+1],
                triPointReverseMap[rawPointi+2],
                0
            );
            rawPointi += 3;

            if (tri.good())
            {
                newToOldTri.append(oldTriI);
                dynTris.append(tri);
            }
        }

        triMap.transfer(newToOldTri);
        tris.transfer(dynTris);
    }



    // Determine 'flat hole' situation (see RMT paper).
    // Two unconnected triangles get connected because (some of) the edges
    // separating them get collapsed. Below only checks for duplicate triangles,
    // not non-manifold edge connectivity.
    if (checkDuplicates)
    {
        labelListList pointFaces;
        invertManyToMany(newPoints.size(), tris, pointFaces);

        // Filter out duplicates.
        DynamicList<label> newToOldTri(tris.size());

        forAll(tris, triI)
        {
            const labelledTri& tri = tris[triI];
            const labelList& pFaces = pointFaces[tri[0]];

            // Find the maximum of any duplicates. Maximum since the tris
            // below triI
            // get overwritten so we cannot use them in a comparison.
            label dupTriI = -1;
            forAll(pFaces, i)
            {
                label nbrTriI = pFaces[i];

                if (nbrTriI > triI && (tris[nbrTriI] == tri))
                {
                    //Pout<< "Duplicate : " << triI << " verts:" << tri
                    //    << " to " << nbrTriI << " verts:" << tris[nbrTriI]
                    //    << endl;
                    dupTriI = nbrTriI;
                    break;
                }
            }

            if (dupTriI == -1)
            {
                // There is no (higher numbered) duplicate triangle
                label newTriI = newToOldTri.size();
                newToOldTri.append(triMap[triI]);
                tris[newTriI] = tris[triI];
            }
        }

        triMap.transfer(newToOldTri);
        tris.setSize(triMap.size());

        if (debug)
        {
            Pout<< "isoSurfacePoint : merged from " << nTris
                << " down to " << tris.size() << " unique triangles." << endl;
        }


        if (debug)
        {
            triSurface surf(tris, geometricSurfacePatchList(), newPoints);

            forAll(surf, facei)
            {
                const labelledTri& f = surf[facei];
                const labelList& fFaces = surf.faceFaces()[facei];

                forAll(fFaces, i)
                {
                    label nbrFacei = fFaces[i];

                    if (nbrFacei <= facei)
                    {
                        // lower numbered faces already checked
                        continue;
                    }

                    const labelledTri& nbrF = surf[nbrFacei];

                    if (f == nbrF)
                    {
                        FatalErrorInFunction
                            << "Check : "
                            << " triangle " << facei << " vertices " << f
                            << " fc:" << f.centre(surf.points())
                            << " has the same vertices as triangle " << nbrFacei
                            << " vertices " << nbrF
                            << " fc:" << nbrF.centre(surf.points())
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    return triSurface(tris, geometricSurfacePatchList(), newPoints, true);
}


void Foam::isoSurfacePoint::trimToPlanes
(
    const PtrList<plane>& planes,
    const triPointRef& tri,
    DynamicList<point>& newTriPoints
)
{
    // Buffer for generated triangles
    DynamicList<triPoints> insideTrisA;
    storeOp insideOpA(insideTrisA);

    // Buffer for generated triangles
    DynamicList<triPoints> insideTrisB;
    storeOp insideOpB(insideTrisB);

    triPointRef::dummyOp dop;

    // Store starting triangle in insideTrisA
    insideOpA(triPoints(tri.a(), tri.b(), tri.c()));


    bool useA = true;

    forAll(planes, faceI)
    {
        const plane& pln = planes[faceI];

        if (useA)
        {
            insideTrisB.clear();
            for (const triPoints& tri : insideTrisA)
            {
                triPointRef(tri).sliceWithPlane(pln, insideOpB, dop);
            }
        }
        else
        {
            insideTrisA.clear();
            for (const triPoints& tri : insideTrisB)
            {
                triPointRef(tri).sliceWithPlane(pln, insideOpA, dop);
            }
        }
        useA = !useA;
    }


    DynamicList<triPoints>& newTris = (useA ? insideTrisA : insideTrisB);

    newTriPoints.reserve(newTriPoints.size() + 3*newTris.size());

    // Transfer
    for (const triPoints& tri : newTris)
    {
        newTriPoints.append(tri[0]);
        newTriPoints.append(tri[1]);
        newTriPoints.append(tri[2]);
    }
}


void Foam::isoSurfacePoint::trimToBox
(
    const treeBoundBox& bb,
    DynamicList<point>& triPoints,
    DynamicList<label>& triMap      // trimmed to original tris
)
{
    if (debug)
    {
        Pout<< "isoSurfacePoint : trimming to " << bb << endl;
    }

    // Generate inwards pointing planes
    PtrList<plane> planes(boundBox::nFaces());
    forAll(treeBoundBox::faceNormals, faceI)
    {
        const vector& n = treeBoundBox::faceNormals[faceI];
        planes.set(faceI, new plane(bb.faceCentre(faceI), -n));
    }

    label nTris = triPoints.size()/3;

    DynamicList<point> newTriPoints(triPoints.size()/16);
    triMap.setCapacity(nTris/16);

    label vertI = 0;
    for (label triI = 0; triI < nTris; triI++)
    {
        const point& p0 = triPoints[vertI++];
        const point& p1 = triPoints[vertI++];
        const point& p2 = triPoints[vertI++];

        label oldNPoints = newTriPoints.size();
        trimToPlanes
        (
            planes,
            triPointRef(p0, p1, p2),
            newTriPoints
        );

        label nCells = (newTriPoints.size()-oldNPoints)/3;
        for (label i = 0; i < nCells; i++)
        {
            triMap.append(triI);
        }
    }

    if (debug)
    {
        Pout<< "isoSurfacePoint : trimmed from " << nTris
            << " down to " << triMap.size()
            << " triangles." << endl;
    }

    triPoints.transfer(newTriPoints);
}


void Foam::isoSurfacePoint::trimToBox
(
    const treeBoundBox& bb,
    DynamicList<point>& triPoints,  // new points
    DynamicList<label>& triMap,     // map from (new) triangle to original
    labelList& triPointMap,         // map from (new) point to original
    labelList& interpolatedPoints,  // labels of newly introduced points
    List<FixedList<label, 3>>& interpolatedOldPoints,// and their interpolation
    List<FixedList<scalar, 3>>& interpolationWeights
)
{
    const List<point> oldTriPoints(triPoints);

    // Trim triPoints, return map
    trimToBox(bb, triPoints, triMap);


    // Find point correspondence:
    // - one-to-one for preserved original points (triPointMap)
    // - interpolation for newly introduced points
    //   (interpolatedOldPoints)
    label sz = oldTriPoints.size()/100;
    DynamicList<label> dynInterpolatedPoints(sz);
    DynamicList<FixedList<label, 3>> dynInterpolatedOldPoints(sz);
    DynamicList<FixedList<scalar, 3>> dynInterpolationWeights(sz);


    triPointMap.setSize(triPoints.size());
    forAll(triMap, triI)
    {
        label oldTriI = triMap[triI];

        // Find point correspondence. Assumes coordinate is bit-copy.
        for (label i = 0; i < 3; i++)
        {
            label pointI = 3*triI+i;
            const point& pt = triPoints[pointI];

            // Compare to old-triangle's points
            label matchPointI = -1;
            for (label j = 0; j < 3; j++)
            {
                label oldPointI = 3*oldTriI+j;
                if (pt == oldTriPoints[oldPointI])
                {
                    matchPointI = oldPointI;
                    break;
                }
            }

            triPointMap[pointI] = matchPointI;

            // If new point: calculate and store interpolation
            if (matchPointI == -1)
            {
                dynInterpolatedPoints.append(pointI);

                FixedList<label, 3> oldPoints
                (
                    {3*oldTriI, 3*oldTriI+1, 3*oldTriI+2}
                );
                dynInterpolatedOldPoints.append(oldPoints);

                triPointRef tri(oldTriPoints, oldPoints);
                barycentric2D bary = tri.pointToBarycentric(pt);
                FixedList<scalar, 3> weights({bary.a(), bary.b(), bary.c()});

                dynInterpolationWeights.append(weights);
            }
        }
    }

    interpolatedPoints.transfer(dynInterpolatedPoints);
    interpolatedOldPoints.transfer(dynInterpolatedOldPoints);
    interpolationWeights.transfer(dynInterpolationWeights);
}


Foam::triSurface Foam::isoSurfacePoint::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldFaces,
    labelList& oldToNewPoints,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        ListOps::createWithValue<bool>(s.size(), newToOldFaces, true, false)
    );

    newToOldPoints.setSize(s.points().size());
    oldToNewPoints.setSize(s.points().size());
    oldToNewPoints = -1;
    {
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Renumber labels for triangle
                const labelledTri& tri = s[oldFacei];

                forAll(tri, fp)
                {
                    label oldPointi = tri[fp];

                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }

    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = s.points()[newToOldPoints[i]];
    }
    // Extract faces
    List<labelledTri> newTriangles(newToOldFaces.size());

    forAll(newToOldFaces, i)
    {
        // Get old vertex labels
        const labelledTri& tri = s[newToOldFaces[i]];

        newTriangles[i][0] = oldToNewPoints[tri[0]];
        newTriangles[i][1] = oldToNewPoints[tri[1]];
        newTriangles[i][2] = oldToNewPoints[tri[2]];
        newTriangles[i].region() = tri.region();
    }

    // Reuse storage.
    return triSurface(newTriangles, s.patches(), newPoints, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfacePoint::isoSurfacePoint
(
    const volScalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const isoSurfaceParams& params,
    const bitSet& /*unused*/
)
:
    isoSurfaceBase
    (
        cellValues.mesh(),
        cellValues.primitiveField(),
        pointValues,
        iso,
        params
    ),
    cValsPtr_(nullptr),
    mergeDistance_(params.mergeTol()*mesh_.bounds().mag())
{
    const bool regularise = (params.filter() != filterType::NONE);
    const fvMesh& fvmesh = cellValues.mesh();

    if (debug)
    {
        Pout<< "isoSurfacePoint:" << nl
            << "    isoField      : " << cellValues.name() << nl
            << "    cell min/max  : " << minMax(cVals_) << nl
            << "    point min/max : " << minMax(pVals_) << nl
            << "    isoValue      : " << iso << nl
            << "    filter        : " << Switch(regularise) << nl
            << "    mergeTol      : " << params.mergeTol() << nl
            << endl;
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Rewrite input field
    // ~~~~~~~~~~~~~~~~~~~
    // Rewrite input volScalarField to have interpolated values
    // on separated patches.

    cValsPtr_.reset(adaptPatchFields(cellValues).ptr());


    // Construct cell centres field consistent with cVals
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Generate field to interpolate. This is identical to the mesh.C()
    // except on separated coupled patches and on empty patches.

    slicedVolVectorField meshC
    (
        IOobject
        (
            "C",
            fvmesh.pointsInstance(),
            polyMesh::meshSubDir,
            fvmesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        fvmesh,
        dimLength,
        fvmesh.cellCentres(),
        fvmesh.faceCentres()
    );
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Adapt separated coupled (proc and cyclic) patches
        if (pp.coupled())
        {
            auto& pfld = const_cast<fvPatchVectorField&>
            (
                meshC.boundaryField()[patchi]
            );

            bitSet isCollocated
            (
                collocatedFaces(refCast<const coupledPolyPatch>(pp))
            );

            forAll(isCollocated, i)
            {
                if (!isCollocated[i])
                {
                    pfld[i] = mesh_.faceCentres()[pp.start()+i];
                }
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {
            auto& bfld = const_cast<slicedVolVectorField::Boundary&>
            (
                meshC.boundaryField()
            );

            // Clear old value. Cannot resize it since is a slice.
            // Set new value we can change

            bfld.release(patchi);
            bfld.set
            (
                patchi,
                new calculatedFvPatchField<vector>
                (
                    fvmesh.boundary()[patchi],
                    meshC
                )
            );

            // Change to face centres
            bfld[patchi] = pp.patchSlice(mesh_.faceCentres());
        }
    }


    // Pre-calculate patch-per-face to avoid whichPatch call.
    labelList boundaryRegion(mesh_.nBoundaryFaces());

    for (const polyPatch& pp : patches)
    {
        SubList<label>(boundaryRegion, pp.size(), pp.offset()) = pp.index();
    }


    // Determine if any cut through face/cell
    calcCutTypes(boundaryRegion, meshC, cValsPtr_(), pVals_);

    if (debug)
    {
        volScalarField debugField
        (
            IOobject
            (
                "isoSurfacePoint.cutType",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvmesh,
            dimensionedScalar(dimless, Zero)
        );

        auto& debugFld = debugField.primitiveFieldRef();

        forAll(cellCutType_, celli)
        {
            debugFld[celli] = cellCutType_[celli];
        }

        Pout<< "Writing cut types:"
            << debugField.objectPath() << endl;

        debugField.write();
    }


    DynamicList<point> snappedPoints(nCutCells_);

    // Per cc -1 or a point inside snappedPoints.
    labelList snappedCc;
    if (regularise)
    {
        calcSnappedCc
        (
            boundaryRegion,
            meshC,
            cValsPtr_(),
            pVals_,

            snappedPoints,
            snappedCc
        );
    }
    else
    {
        snappedCc.setSize(mesh_.nCells());
        snappedCc = -1;
    }



    if (debug)
    {
        Pout<< "isoSurfacePoint : shifted " << snappedPoints.size()
            << " cell centres to intersection." << endl;
    }

    label nCellSnaps = snappedPoints.size();


    // Per point -1 or a point inside snappedPoints.
    labelList snappedPoint;
    if (regularise)
    {
        // Determine if point is on boundary.
        bitSet isBoundaryPoint(mesh_.nPoints());

        for (const polyPatch& pp : patches)
        {
            // Mark all boundary points that are not physically coupled
            // (so anything but collocated coupled patches)

            if (pp.coupled())
            {
                const coupledPolyPatch& cpp =
                    refCast<const coupledPolyPatch>(pp);

                bitSet isCollocated(collocatedFaces(cpp));

                forAll(isCollocated, i)
                {
                    if (!isCollocated[i])
                    {
                        const face& f = mesh_.faces()[cpp.start()+i];

                        isBoundaryPoint.set(f);
                    }
                }
            }
            else
            {
                forAll(pp, i)
                {
                    const face& f = mesh_.faces()[pp.start()+i];

                    isBoundaryPoint.set(f);
                }
            }
        }

        calcSnappedPoint
        (
            isBoundaryPoint,
            boundaryRegion,
            meshC,
            cValsPtr_(),
            pVals_,

            snappedPoints,
            snappedPoint
        );
    }
    else
    {
        snappedPoint.setSize(mesh_.nPoints());
        snappedPoint = -1;
    }

    if (debug)
    {
        Pout<< "isoSurfacePoint : shifted " << snappedPoints.size()-nCellSnaps
            << " vertices to intersection." << endl;
    }


    // Use a triSurface as a temporary for various operations
    triSurface tmpsurf;

    {
        DynamicList<point> triPoints(3*nCutCells_);
        DynamicList<label> triMeshCells(nCutCells_);

        generateTriPoints
        (
            cValsPtr_(),
            pVals_,

            meshC,
            mesh_.points(),

            snappedPoints,
            snappedCc,
            snappedPoint,

            triPoints,      // 3 points of the triangle
            triMeshCells    // per triangle the originating cell
        );

        if (debug)
        {
            Pout<< "isoSurfacePoint : generated " << triMeshCells.size()
                << " unmerged triangles from " << triPoints.size()
                << " unmerged points." << endl;
        }

        label nOldPoints = triPoints.size();

        // Trimmed to original triangle
        DynamicList<label> trimTriMap;
        // Trimmed to original point
        labelList trimTriPointMap;
        if (getClipBounds().good())
        {
            trimToBox
            (
                treeBoundBox(getClipBounds()),
                triPoints,              // new points
                trimTriMap,             // map from (new) triangle to original
                trimTriPointMap,        // map from (new) point to original
                interpolatedPoints_,    // labels of newly introduced points
                interpolatedOldPoints_, // and their interpolation
                interpolationWeights_
            );
            triMeshCells = labelField(triMeshCells, trimTriMap);
        }


        // Merge points and compact out non-valid triangles
        labelList triMap;           // merged to unmerged triangle
        tmpsurf = stitchTriPoints
        (
            true,               // check for duplicate tris
            triPoints,
            triPointMergeMap_,  // unmerged to merged point
            triMap
        );

        if (debug)
        {
            Pout<< "isoSurfacePoint : generated " << triMap.size()
                << " merged triangles." << endl;
        }


        if (getClipBounds().good())
        {
            // Adjust interpolatedPoints_
            inplaceRenumber(triPointMergeMap_, interpolatedPoints_);

            // Adjust triPointMergeMap_
            labelList newTriPointMergeMap(nOldPoints, -1);
            forAll(trimTriPointMap, trimPointI)
            {
                label oldPointI = trimTriPointMap[trimPointI];
                if (oldPointI >= 0)
                {
                    label pointI = triPointMergeMap_[trimPointI];
                    if (pointI >= 0)
                    {
                        newTriPointMergeMap[oldPointI] = pointI;
                    }
                }
            }
            triPointMergeMap_.transfer(newTriPointMergeMap);
        }

        meshCells_.setSize(triMap.size());
        forAll(triMap, i)
        {
            meshCells_[i] = triMeshCells[triMap[i]];
        }
    }

    if (debug)
    {
        Pout<< "isoSurfacePoint : checking " << tmpsurf.size()
            << " triangles for validity." << endl;

        forAll(tmpsurf, facei)
        {
            triSurfaceTools::validTri(tmpsurf, facei);
        }

        fileName stlFile = mesh_.time().path() + ".stl";
        Pout<< "Dumping surface to " << stlFile << endl;
        tmpsurf.write(stlFile);
    }


    // Transfer to mesh storage. Note, an iso-surface has no zones
    {
        // Recover the pointField
        pointField pts;
        tmpsurf.swapPoints(pts);

        // Transcribe from triFace to face
        faceList faces;
        tmpsurf.triFaceFaces(faces);

        tmpsurf.clearOut();

        Mesh updated(std::move(pts), std::move(faces), surfZoneList());

        this->Mesh::transfer(updated);
    }
}


// ************************************************************************* //
