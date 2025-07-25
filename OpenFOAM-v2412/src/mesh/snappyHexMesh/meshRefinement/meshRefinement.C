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

#include "meshRefinement.H"
#include "volMesh.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "decompositionMethod.H"
#include "regionSplit.H"
#include "fvMeshDistribute.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "mapDistributePolyMesh.H"
#include "localPointRegion.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "calculatedPointPatchFields.H"
#include "cyclicSlipPointPatchFields.H"
#include "processorPointPatch.H"
#include "globalIndex.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "searchableSurfaces.H"
#include "treeBoundBox.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMeshTools.H"
#include "motionSmoother.H"
#include "faceSet.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"
#include "PackedBoolList.H"

// Leak path
#include "shortestPathSet.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshRefinement, 0);
}


const Foam::Enum
<
    Foam::meshRefinement::MeshType
>
Foam::meshRefinement::MeshTypeNames
({
    { MeshType::CASTELLATED, "castellated" },
    { MeshType::CASTELLATEDBUFFERLAYER, "castellatedBufferLayer" },
    { MeshType::CASTELLATEDBUFFERLAYER2, "castellatedBufferLayer2" }
});


const Foam::Enum
<
    Foam::meshRefinement::debugType
>
Foam::meshRefinement::debugTypeNames
({
    { debugType::MESH, "mesh" },
    { debugType::OBJINTERSECTIONS, "intersections" },
    { debugType::FEATURESEEDS, "featureSeeds" },
    { debugType::ATTRACTION, "attraction" },
    { debugType::LAYERINFO, "layerInfo" },
});


//const Foam::Enum
//<
//    Foam::meshRefinement::outputType
//>
//Foam::meshRefinement::outputTypeNames
//({
//    { outputType::OUTPUTLAYERINFO, "layerInfo" }
//});


const Foam::Enum
<
    Foam::meshRefinement::writeType
>
Foam::meshRefinement::writeTypeNames
({
    { writeType::WRITEMESH, "mesh" },
    { writeType::NOWRITEREFINEMENT, "noRefinement" },
    { writeType::WRITELEVELS, "scalarLevels" },
    { writeType::WRITELAYERSETS, "layerSets" },
    { writeType::WRITELAYERFIELDS, "layerFields" },
});


Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel_;

//Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel_;

// Inside/outside test for polyMesh:.findCell()
//   2.4.x : default = polyMesh::FACE_DIAG_TRIS
//   1712  : default = polyMesh::CELL_TETS
//
// - CELL_TETS is better with concave cells, but much slower.
// - use faster method (FACE_DIAG_TRIS) here

static const Foam::polyMesh::cellDecomposition
    findCellMode(Foam::polyMesh::FACE_DIAG_TRIS);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::globalFaceCount(const labelList& elems) const
{
    //- Check for duplicates
    const bitSet isElem(mesh_.nFaces(), elems);
    if (label(isElem.count()) != elems.size())
    {
        FatalErrorInFunction << "Problem Duplicates:"
            << " isElem:" << isElem.count()
            << " elems:" << elems.size()
            << exit(FatalError);
    }

    //- Check for same entries on coupled faces
    {
        bitSet isElem2(isElem);
        syncTools::swapFaceList(mesh_, isElem2);

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (isElem2[facei] != isElem[facei])
            {
                FatalErrorInFunction
                    << "at face:" << facei
                    << " at:" << mesh_.faceCentres()[facei]
                    << " patch:" << mesh_.boundaryMesh().whichPatch(facei)
                    << " isElem:" << isElem[facei]
                    << " isElem2:" << isElem2[facei]
                    << exit(FatalError);
            }
        }
    }

    //- Count number of master elements
    const bitSet isMaster(syncTools::getMasterFaces(mesh_));
    label count = 0;
    for (const label i : isElem)
    {
        if (isMaster[i])
        {
            count++;
        }
    }
    return returnReduce(count, sumOp<label>());
}


void Foam::meshRefinement::calcNeighbourData
(
    labelList& neiLevel,
    pointField& neiCc
)  const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    const label nBoundaryFaces = mesh_.nBoundaryFaces();

    if (neiLevel.size() != nBoundaryFaces || neiCc.size() != nBoundaryFaces)
    {
        FatalErrorInFunction
            << nBoundaryFaces << " neiLevel:" << neiLevel.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet addedPatchIDSet(meshedPatches());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelUList& faceCells = pp.faceCells();
        const vectorField::subField faceCentres = pp.faceCentres();
        const vectorField::subField faceAreas = pp.faceAreas();

        label bFacei = pp.start()-mesh_.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = cellCentres[faceCells[i]];
                bFacei++;
            }
        }
        else if (addedPatchIDSet.found(patchi))
        {
            // Face was introduced from cell-cell intersection. Try to
            // reconstruct other side cell(centre). Three possibilities:
            // - cells same size.
            // - preserved cell smaller. Not handled.
            // - preserved cell larger.
            forAll(faceCells, i)
            {
                // Extrapolate the face centre.
                const vector fn = normalised(faceAreas[i]);

                label own = faceCells[i];
                label ownLevel = cellLevel[own];
                label faceLevel = meshCutter_.faceLevel(pp.start()+i);
                if (faceLevel < 0)
                {
                    // Due to e.g. face merging no longer a consistent
                    // refinementlevel of face. Assume same as cell
                    faceLevel = ownLevel;
                }

                // Normal distance from face centre to cell centre
                scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                if (faceLevel > ownLevel)
                {
                    // Other cell more refined. Adjust normal distance
                    d *= 0.5;
                }
                neiLevel[bFacei] = faceLevel;
                // Calculate other cell centre by extrapolation
                neiCc[bFacei] = faceCentres[i] + d*fn;
                bFacei++;
            }
        }
        else
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = faceCentres[i];
                bFacei++;
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh_, neiCc);
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);
}


void Foam::meshRefinement::calcCellCellRays
(
    const pointField& neiCc,
    const labelList& neiLevel,
    const labelList& testFaces,
    pointField& start,
    pointField& end,
    labelList& minLevel
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();


    // Mark all non-coupled or coupled+master faces. Leaves only slave of
    // coupled unset.
    bitSet isMaster(mesh_.nBoundaryFaces(), true);
    {
        for (const polyPatch& pp : mesh_.boundaryMesh())
        {
            if (pp.coupled() && !refCast<const coupledPolyPatch>(pp).owner())
            {
                isMaster.unset(labelRange(pp.offset(), pp.size()));
            }
        }
    }


    start.setSize(testFaces.size());
    end.setSize(testFaces.size());
    minLevel.setSize(testFaces.size());

    forAll(testFaces, i)
    {
        const label facei = testFaces[i];
        const label own = mesh_.faceOwner()[facei];

        if (mesh_.isInternalFace(facei))
        {
            const label nei = mesh_.faceNeighbour()[facei];

            start[i] = cellCentres[own];
            end[i] = cellCentres[nei];
            minLevel[i] = min(cellLevel[own], cellLevel[nei]);
        }
        else
        {
            const label bFacei = facei - mesh_.nInternalFaces();

            if (isMaster[bFacei])
            {
                start[i] = cellCentres[own];
                end[i] = neiCc[bFacei];
            }
            else
            {
                // Slave face
                start[i] = neiCc[bFacei];
                end[i] = cellCentres[own];
            }
            minLevel[i] = min(cellLevel[own], neiLevel[bFacei]);
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }
}


void Foam::meshRefinement::updateIntersections(const labelList& changedFaces)
{
    // Stats on edges to test. Count proc faces only once.
    bitSet isMasterFace(syncTools::getMasterFaces(mesh_));

    {
        const label nMasterFaces =
            returnReduce(isMasterFace.count(), sumOp<label>());

        label nChangedFaces = 0;
        forAll(changedFaces, i)
        {
            if (isMasterFace.test(changedFaces[i]))
            {
                ++nChangedFaces;
            }
        }
        reduce(nChangedFaces, sumOp<label>());

        if (!dryRun_)
        {
            Info<< "Edge intersection testing:" << nl
                << "    Number of edges             : " << nMasterFaces << nl
                << "    Number of edges to retest   : " << nChangedFaces
                << endl;
        }
    }


    // Get boundary face centre and level. Coupled aware.
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);

    // Collect segments we want to test for
    pointField start(changedFaces.size());
    pointField end(changedFaces.size());
    {
        labelList minLevel;
        calcCellCellRays
        (
            neiCc,
            neiLevel,
            changedFaces,
            start,
            end,
            minLevel
        );
    }


    // Do tests in one go
    labelList surfaceHit;
    {
        labelList surfaceLevel;
        surfaces_.findHigherIntersection
        (
            shells_,
            start,
            end,
            labelList(start.size(), -1),    // accept any intersection
            surfaceHit,
            surfaceLevel
        );
    }

    // Keep just surface hit
    forAll(surfaceHit, i)
    {
        surfaceIndex_[changedFaces[i]] = surfaceHit[i];
    }

    // Make sure both sides have same information. This should be
    // case in general since same vectors but just to make sure.
    syncTools::syncFaceList(mesh_, surfaceIndex_, maxEqOp<label>());

    label nTotHits = returnReduce(countHits(), sumOp<label>());

    if (!dryRun_)
    {
        Info<< "    Number of intersected edges : " << nTotHits << endl;
    }

    // Set files to same time as mesh
    setInstance(mesh_.facesInstance());
}


void Foam::meshRefinement::nearestFace
(
    const labelUList& startFaces,
    const bitSet& isBlockedFace,

    autoPtr<mapDistribute>& mapPtr,
    labelList& faceToStart,
    const label nIter
) const
{
    // From startFaces walk out (but not through isBlockedFace). Returns
    // faceToStart which is the index into startFaces (or rather distributed
    // version of it). E.g.
    // pointField startFc(mesh.faceCentres(), startFaces);
    // mapPtr().distribute(startFc);
    // forAll(faceToStart, facei)
    // {
    //    const label sloti = faceToStart[facei];
    //    if (sloti != -1)
    //    {
    //        Pout<< "face:" << mesh.faceCentres()[facei]
    //            << " nearest:" << startFc[sloti]
    //            << endl;
    //    }
    // }


    // Consecutive global numbering for start elements
    const globalIndex globalStart(startFaces.size());

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh_.nCells());
    List<topoDistanceData<label>> faceData(mesh_.nFaces());

    // Mark blocked faces to there not visited
    for (const label facei : isBlockedFace)
    {
        faceData[facei] = topoDistanceData<label>(0, -1);
    }

    List<topoDistanceData<label>> startData(startFaces.size());
    forAll(startFaces, i)
    {
        const label facei = startFaces[i];
        if (isBlockedFace[facei])
        {
            FatalErrorInFunction << "Start face:" << facei
                << " at:" << mesh_.faceCentres()[facei]
                << " is also blocked" << exit(FatalError);
        }
        startData[i] = topoDistanceData<label>(0, globalStart.toGlobal(i));
    }


    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
    (
        mesh_,
        startFaces,
        startData,
        faceData,
        cellData,
        0
    );
    deltaCalc.iterate(nIter);

    // And extract

    faceToStart.setSize(mesh_.nFaces());
    faceToStart = -1;
    bool haveWarned = false;
    forAll(faceData, facei)
    {
        if (!faceData[facei].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some faces, e.g. face " << facei
                    << " at " << mesh_.faceCentres()[facei]
                    << endl;
                haveWarned = true;
            }
        }
        else
        {
            faceToStart[facei] = faceData[facei].data();
        }
    }

    // Compact
    List<Map<label>> compactMap;
    mapPtr.reset(new mapDistribute(globalStart, faceToStart, compactMap));
}


void Foam::meshRefinement::nearestPatch
(
    const labelList& adaptPatchIDs,
    labelList& nearestPatch,
    labelList& nearestZone
) const
{
    // Determine nearest patch for all mesh faces. Used when removing cells
    // to give some reasonable patch to exposed faces.

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    nearestZone.setSize(mesh_.nFaces(), -1);

    if (adaptPatchIDs.size())
    {
        nearestPatch.setSize(mesh_.nFaces(), adaptPatchIDs[0]);


        // Get per-face the zone or -1
        labelList faceToZone(mesh_.nFaces(), -1);
        {
            for (const faceZone& fz : mesh_.faceZones())
            {
                UIndirectList<label>(faceToZone, fz) = fz.index();
            }
        }


        // Count number of faces in adaptPatchIDs
        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = patches[adaptPatchIDs[i]];
            nFaces += pp.size();
        }

        // Field on cells and faces.
        List<topoDistanceData<labelPair>> cellData(mesh_.nCells());
        List<topoDistanceData<labelPair>> faceData(mesh_.nFaces());

        // Start of changes
        labelList patchFaces(nFaces);
        List<topoDistanceData<labelPair>> patchData(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            label patchi = adaptPatchIDs[i];
            const polyPatch& pp = patches[patchi];

            forAll(pp, i)
            {
                patchFaces[nFaces] = pp.start()+i;
                patchData[nFaces] = topoDistanceData<labelPair>
                (
                    0,
                    labelPair
                    (
                        patchi,
                        faceToZone[pp.start()+i]
                    )
                );
                nFaces++;
            }
        }

        // Propagate information inwards
        FaceCellWave<topoDistanceData<labelPair>> deltaCalc
        (
            mesh_,
            patchFaces,
            patchData,
            faceData,
            cellData,
            mesh_.globalData().nTotalCells()+1
        );

        // And extract

        bool haveWarned = false;
        forAll(faceData, facei)
        {
            if (!faceData[facei].valid(deltaCalc.data()))
            {
                if (!haveWarned)
                {
                    WarningInFunction
                        << "Did not visit some faces, e.g. face " << facei
                        << " at " << mesh_.faceCentres()[facei] << endl
                        << "Assigning these faces to patch "
                        << adaptPatchIDs[0]
                        << endl;
                    haveWarned = true;
                }
            }
            else
            {
                const labelPair& data = faceData[facei].data();
                nearestPatch[facei] = data.first();
                nearestZone[facei] = data.second();
            }
        }
    }
    else
    {
        // Use patch 0
        nearestPatch.setSize(mesh_.nFaces(), 0);
    }
}


Foam::labelList Foam::meshRefinement::nearestPatch
(
    const labelList& adaptPatchIDs
) const
{
    labelList nearestAdaptPatch;
    labelList nearestAdaptZone;
    nearestPatch(adaptPatchIDs, nearestAdaptPatch, nearestAdaptZone);
    return nearestAdaptPatch;
}


void Foam::meshRefinement::nearestIntersection
(
    const labelList& surfacesToTest,
    const labelList& testFaces,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2
) const
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nBoundaryFaces());
    pointField neiCc(mesh_.nBoundaryFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Collect segments

    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    calcCellCellRays
    (
        neiCc,
        neiLevel,
        testFaces,
        start,
        end,
        minLevel
    );

    // Do tests in one go

    surfaces_.findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1,
        hit1,
        region1,
        surface2,
        hit2,
        region2
    );
}


Foam::labelList Foam::meshRefinement::nearestIntersection
(
    const labelList& surfacesToTest,
    const label defaultRegion
) const
{
    // Determine nearest intersection for all mesh faces. Used when removing
    // cells to give some reasonable patch to exposed faces. Use this
    // function instead of nearestPatch if you don't have patches yet.


    // All faces with any intersection
    const labelList testFaces(intersectedFaces());

    // Find intersection (if any)
    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    nearestIntersection
    (
        surfacesToTest,
        testFaces,

        surface1,
        hit1,
        region1,
        surface2,
        hit2,
        region2
    );


    labelList nearestRegion(mesh_.nFaces(), defaultRegion);

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh_.nCells());
    List<topoDistanceData<label>> faceData(mesh_.nFaces());

    // Start walking from all intersected faces
    DynamicList<label> patchFaces(testFaces.size());
    DynamicList<topoDistanceData<label>> patchData(testFaces.size());
    forAll(testFaces, i)
    {
        label facei = testFaces[i];
        if (surface1[i] != -1)
        {
            patchFaces.append(facei);
            label regioni = surfaces_.globalRegion(surface1[i], region1[i]);
            patchData.append(topoDistanceData<label>(0, regioni));
        }
        else if (surface2[i] != -1)
        {
            patchFaces.append(facei);
            label regioni = surfaces_.globalRegion(surface2[i], region2[i]);
            patchData.append(topoDistanceData<label>(0, regioni));
        }
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
    (
        mesh_,
        patchFaces,
        patchData,
        faceData,
        cellData,
        mesh_.globalData().nTotalCells()+1
    );

    // And extract

    bool haveWarned = false;
    forAll(faceData, facei)
    {
        if (!faceData[facei].valid(deltaCalc.data()))
        {
            if (!haveWarned)
            {
                WarningInFunction
                    << "Did not visit some faces, e.g. face " << facei
                    << " at " << mesh_.faceCentres()[facei] << endl
                    << "Assigning these faces to global region "
                    << defaultRegion << endl;
                haveWarned = true;
            }
        }
        else
        {
            nearestRegion[facei] = faceData[facei].data();
        }
    }

    return nearestRegion;
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<scalar>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    scalarField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minEqOp<scalar>(),
        GREAT
    );
    scalarField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxEqOp<scalar>(),
        -GREAT
    );
    forAll(minFld, pointi)
    {
        const scalar& minVal = minFld[pointi];
        const scalar& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<point>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    pointField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minMagSqrEqOp<point>(),
        point(GREAT, GREAT, GREAT)
    );
    pointField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    forAll(minFld, pointi)
    {
        const point& minVal = minFld[pointi];
        const point& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::checkData()
{
    Pout<< "meshRefinement::checkData() : Checking refinement structure."
        << endl;
    meshCutter_.checkMesh();

    Pout<< "meshRefinement::checkData() : Checking refinement levels."
        << endl;
    meshCutter_.checkRefinementLevels(1, labelList(0));


    const label nBnd = mesh_.nBoundaryFaces();

    Pout<< "meshRefinement::checkData() : Checking synchronization."
        << endl;

    // Check face centres
    {
        // Boundary face centres
        pointField::subList boundaryFc
        (
            mesh_.faceCentres(),
            nBnd,
            mesh_.nInternalFaces()
        );

        // Get neighbouring face centres
        pointField neiBoundaryFc(boundaryFc);
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            neiBoundaryFc,
            eqOp<point>()
        );

        // Compare
        testSyncBoundaryFaceList
        (
            mergeDistance_,
            "testing faceCentres : ",
            boundaryFc,
            neiBoundaryFc
        );
    }

    // Check meshRefinement
    const labelList& surfIndex = surfaceIndex();


    {
        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(nBnd);
        pointField neiCc(nBnd);
        calcNeighbourData(neiLevel, neiCc);

        // Collect segments we want to test for
        pointField start(mesh_.nFaces());
        pointField end(mesh_.nFaces());

        forAll(start, facei)
        {
            start[facei] = mesh_.cellCentres()[mesh_.faceOwner()[facei]];

            if (mesh_.isInternalFace(facei))
            {
                end[facei] = mesh_.cellCentres()[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                end[facei] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do tests in one go
        labelList surfaceHit;
        {
            labelList surfaceLevel;
            surfaces_.findHigherIntersection
            (
                shells_,
                start,
                end,
                labelList(start.size(), -1),    // accept any intersection
                surfaceHit,
                surfaceLevel
            );
        }
        // Get the coupled hit
        labelList neiHit
        (
            SubList<label>
            (
                surfaceHit,
                nBnd,
                mesh_.nInternalFaces()
            )
        );
        syncTools::swapBoundaryFaceList(mesh_, neiHit);

        // Check
        forAll(surfaceHit, facei)
        {
            if (surfIndex[facei] != surfaceHit[facei])
            {
                if (mesh_.isInternalFace(facei))
                {
                    WarningInFunction
                        << "Internal face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfIndex[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " neiCc:"
                        << mesh_.cellCentres()[mesh_.faceNeighbour()[facei]]
                        << endl;
                }
                else if
                (
                    surfIndex[facei]
                 != neiHit[facei-mesh_.nInternalFaces()]
                )
                {
                    WarningInFunction
                        << "Boundary face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfIndex[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " neiCc:"
                        << neiCc[facei-mesh_.nInternalFaces()]
                        << " end:" << end[facei]
                        << " ownLevel:"
                        << meshCutter_.cellLevel()[mesh_.faceOwner()[facei]]
                        << " faceLevel:"
                        << meshCutter_.faceLevel(facei)
                        << endl;
                }
            }
        }
    }
    {
        labelList::subList boundarySurface
        (
            surfaceIndex_,
            mesh_.nBoundaryFaces(),
            mesh_.nInternalFaces()
        );

        labelList neiBoundarySurface(boundarySurface);
        syncTools::swapBoundaryFaceList
        (
            mesh_,
            neiBoundarySurface
        );

        // Compare
        testSyncBoundaryFaceList
        (
            0,                              // tolerance
            "testing surfaceIndex() : ",
            boundarySurface,
            neiBoundarySurface
        );
    }


    // Find duplicate faces
    Pout<< "meshRefinement::checkData() : Counting duplicate faces."
        << endl;

    labelList duplicateFace
    (
        localPointRegion::findDuplicateFaces
        (
            mesh_,
            identity(mesh_.nBoundaryFaces(), mesh_.nInternalFaces())
        )
    );

    // Count
    {
        label nDup = 0;

        forAll(duplicateFace, i)
        {
            if (duplicateFace[i] != -1)
            {
                nDup++;
            }
        }
        nDup /= 2;  // will have counted both faces of duplicate
        Pout<< "meshRefinement::checkData() : Found " << nDup
            << " duplicate pairs of faces." << endl;
    }
}


void Foam::meshRefinement::setInstance(const fileName& inst)
{
    meshCutter_.setInstance(inst);
    surfaceIndex_.instance() = inst;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::doRemoveCells
(
    const labelList& cellsToRemove,
    const labelList& exposedFaces,
    const labelList& exposedPatchIDs,
    removeCells& cellRemover
)
{
    polyTopoChange meshMod(mesh_);

    // Arbitrary: put exposed faces into last patch.
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    // Remove any unnecessary fields
    mesh_.clearOut();
    mesh_.moving(false);

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false, true);
    mapPolyMesh& map = *mapPtr;

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map.hasMotionPoints())
    {
        mesh_.movePoints(map.preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());
    setInstance(mesh_.facesInstance());

    // Update local mesh data
    cellRemover.updateMesh(map);

    // Update intersections. Recalculate intersections for exposed faces.
    labelList newExposedFaces = renumber
    (
        map.reverseFaceMap(),
        exposedFaces
    );

    //Pout<< "removeCells : updating intersections for "
    //    << newExposedFaces.size() << " newly exposed faces." << endl;

    updateMesh(map, newExposedFaces);

    return mapPtr;
}


void Foam::meshRefinement::splitFace
(
    const face& f,
    const labelPair& split,

    face& f0,
    face& f1
)
{
    if (split.find(-1) != -1)
    {
        FatalErrorInFunction<< "Illegal split " << split
            << " on face " << f << exit(FatalError);
    }

    label nVerts = split[1]-split[0];
    if (nVerts < 0)
    {
        nVerts += f.size();
    }
    nVerts += 1;


    // Split into f0, f1
    f0.resize_nocopy(nVerts);

    label fp = split[0];
    forAll(f0, i)
    {
        f0[i] = f[fp];
        fp = f.fcIndex(fp);
    }

    f1.resize_nocopy(f.size()-f0.size()+2);
    fp = split[1];
    forAll(f1, i)
    {
        f1[i] = f[fp];
        fp = f.fcIndex(fp);
    }
}


void Foam::meshRefinement::doSplitFaces
(
    const labelList& splitFaces,
    const labelPairList& splits,
    const labelPairList& splitPatches,
    //const List<Pair<point>>& splitPoints,
    polyTopoChange& meshMod
) const
{
    face f0, f1;

    forAll(splitFaces, i)
    {
        const label facei = splitFaces[i];
        const auto& split = splits[i];
        const auto& twoPatches = splitPatches[i];

        const face& f = mesh_.faces()[facei];

        // Split as start and end index in face
        splitFace(f, split, f0, f1);

        // Determine face properties
        label own = mesh_.faceOwner()[facei];
        label nei = -1;
        label patch0 = -1;
        label patch1 = -1;
        if (facei >= mesh_.nInternalFaces())
        {
            patch0 =
            (
                twoPatches[0] != -1
              ? twoPatches[0]
              : mesh_.boundaryMesh().whichPatch(facei)
            );
            patch1 =
            (
                twoPatches[1] != -1
              ? twoPatches[1]
              : mesh_.boundaryMesh().whichPatch(facei)
            );
        }
        else
        {
            nei = mesh_.faceNeighbour()[facei];
        }

        label zonei = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zonei != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zonei];
            zoneFlip = fz.flipMap()[fz.whichFace(facei)];
        }


        if (debug)
        {
            Pout<< "face:" << facei << " verts:" << f
                << " split into f0:" << f0
                << " f1:" << f1 << endl;
        }

        // Change/add faces
        meshMod.modifyFace
        (
            f0,                         // modified face
            facei,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patch0,                     // patch for face
            zonei,                      // zone for face
            zoneFlip                    // face flip in zone
        );

        meshMod.addFace
        (
            f1,                         // modified face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master point
            -1,                         // master edge
            facei,                      // master face
            false,                      // face flip
            patch1,                     // patch for face
            zonei,                      // zone for face
            zoneFlip                    // face flip in zone
        );


        //// Move points
        //meshMod.modifyPoint
        //(
        //    f[split[0]],
        //    splitPoints[i][0],
        //    -1,
        //    true
        //);
        //meshMod.modifyPoint
        //(
        //    f[split[1]],
        //    splitPoints[i][1],
        //    -1,
        //    true
        //);
    }
}


Foam::label Foam::meshRefinement::splitFacesUndo
(
    const labelList& splitFaces,
    const labelPairList& splits,
    const labelPairList& splitPatches,
    const dictionary& motionDict,

    labelList& duplicateFace,
    List<labelPair>& baffles
)
{
    label nSplit = returnReduce(splitFaces.size(), sumOp<label>());
    Info<< nl
        << "Splitting " << nSplit
        << " faces across diagonals" << "." << nl << endl;

    // To undo: store original faces
    faceList originalFaces(UIndirectList<face>(mesh_.faces(), splitFaces));
    labelPairList facePairs(splitFaces.size(), labelPair(-1, -1));


    {
        polyTopoChange meshMod(mesh_);
        meshMod.setCapacity
        (
            meshMod.points().size(),
            meshMod.faces().size()+splitFaces.size(),
            mesh_.nCells()
        );

        // Insert the mesh changes
        doSplitFaces(splitFaces, splits, splitPatches, meshMod);

        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());


        // Update local mesh data
        // ~~~~~~~~~~~~~~~~~~~~~~

        forAll(originalFaces, i)
        {
            inplaceRenumber(map.reversePointMap(), originalFaces[i]);
        }

        {
            Map<label> splitFaceToIndex(invertToMap(splitFaces));

            forAll(map.faceMap(), facei)
            {
                label oldFacei = map.faceMap()[facei];

                const auto oldFaceFnd = splitFaceToIndex.cfind(oldFacei);

                if (oldFaceFnd.good())
                {
                    labelPair& twoFaces = facePairs[oldFaceFnd.val()];
                    if (twoFaces[0] == -1)
                    {
                        twoFaces[0] = facei;
                    }
                    else if (twoFaces[1] == -1)
                    {
                        twoFaces[1] = facei;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "problem: twoFaces:" << twoFaces
                            << exit(FatalError);
                    }
                }
            }
        }


        // Update baffle data
        // ~~~~~~~~~~~~~~~~~~

        if (duplicateFace.size())
        {
            meshRefinement::updateList
            (
                map.faceMap(),
                label(-1),
                duplicateFace
            );
        }

        const labelList& oldToNewFaces = map.reverseFaceMap();
        forAll(baffles, i)
        {
            labelPair& baffle = baffles[i];
            baffle.first() = oldToNewFaces[baffle.first()];
            baffle.second() = oldToNewFaces[baffle.second()];

            if (baffle.first() == -1 || baffle.second() == -1)
            {
                FatalErrorInFunction
                    << "Removed baffle : faces:" << baffle
                    << exit(FatalError);
            }
        }


        // Update intersections
        // ~~~~~~~~~~~~~~~~~~~~

        {
            DynamicList<label> changedFaces(facePairs.size());
            forAll(facePairs, i)
            {
                changedFaces.append(facePairs[i].first());
                changedFaces.append(facePairs[i].second());
            }

            // Update intersections on changed faces
            updateMesh(map, growFaceCellFace(changedFaces));
        }
    }



    // Undo loop
    // ~~~~~~~~~
    // Maintains two bits of information:
    // facePairs     : two faces originating from the same face
    // originalFaces : original face in current vertices


    for (label iteration = 0; iteration < 100; iteration++)
    {
        Info<< nl
            << "Undo iteration " << iteration << nl
            << "----------------" << endl;


        // Check mesh for errors
        // ~~~~~~~~~~~~~~~~~~~~~

        faceSet errorFaces(mesh_, "errorFaces", mesh_.nBoundaryFaces());
        bool hasErrors = motionSmoother::checkMesh
        (
            false,  // report
            mesh_,
            motionDict,
            errorFaces,
            dryRun_
        );
        if (!hasErrors)
        {
            break;
        }

        // Extend faces
        {
            const labelList grownFaces(growFaceCellFace(errorFaces));
            errorFaces.clear();
            errorFaces.insert(grownFaces);
        }


        // Merge face pairs
        // ~~~~~~~~~~~~~~~~
        // (if one of the faces is in the errorFaces set)

        polyTopoChange meshMod(mesh_);

        // Indices (in facePairs) of merged faces
        labelHashSet mergedIndices(facePairs.size());
        forAll(facePairs, index)
        {
            const labelPair& twoFaces = facePairs[index];

            if
            (
                errorFaces.found(twoFaces.first())
             || errorFaces.found(twoFaces.second())
            )
            {
                const face& originalFace = originalFaces[index];


                // Determine face properties
                label own = mesh_.faceOwner()[twoFaces[0]];
                label nei = -1;
                label patchi = -1;
                if (twoFaces[0] >= mesh_.nInternalFaces())
                {
                    patchi = mesh_.boundaryMesh().whichPatch(twoFaces[0]);
                }
                else
                {
                    nei = mesh_.faceNeighbour()[twoFaces[0]];
                }

                label zonei = mesh_.faceZones().whichZone(twoFaces[0]);
                bool zoneFlip = false;
                if (zonei != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zonei];
                    zoneFlip = fz.flipMap()[fz.whichFace(twoFaces[0])];
                }

                // Modify first face
                meshMod.modifyFace
                (
                    originalFace,               // modified face
                    twoFaces[0],                // label of face
                    own,                        // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchi,                     // patch for face
                    zonei,                      // zone for face
                    zoneFlip                    // face flip in zone
                );
                // Merge second face into first
                meshMod.removeFace(twoFaces[1], twoFaces[0]);

                mergedIndices.insert(index);
            }
        }

        label n = returnReduce(mergedIndices.size(), sumOp<label>());

        Info<< "Detected " << n
            << " split faces that will be restored to their original faces."
            << nl << endl;

        if (n == 0)
        {
            // Nothing to be restored
            break;
        }

        nSplit -= n;


        // Remove any unnecessary fields
        mesh_.clearOut();
        mesh_.moving(false);

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false, true);
        mapPolyMesh& map = *mapPtr;

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map.hasMotionPoints())
        {
            mesh_.movePoints(map.preMotionPoints());
        }
        else
        {
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());


        // Update local mesh data
        // ~~~~~~~~~~~~~~~~~~~~~~

        {
            const labelList& oldToNewFaces = map.reverseFaceMap();
            const labelList& oldToNewPoints = map.reversePointMap();

            // Compact out merged faces
            DynamicList<label> changedFaces(mergedIndices.size());

            label newIndex = 0;
            forAll(facePairs, index)
            {
                const labelPair& oldSplit = facePairs[index];
                label new0 = oldToNewFaces[oldSplit[0]];
                label new1 = oldToNewFaces[oldSplit[1]];

                if (!mergedIndices.found(index))
                {
                    // Faces still split
                    if (new0 < 0 || new1 < 0)
                    {
                        FatalErrorInFunction
                            << "Problem: oldFaces:" << oldSplit
                            << " newFaces:" << labelPair(new0, new1)
                            << exit(FatalError);
                    }

                    facePairs[newIndex] = labelPair(new0, new1);
                    originalFaces[newIndex] = renumber
                    (
                        oldToNewPoints,
                        originalFaces[index]
                    );
                    newIndex++;
                }
                else
                {
                    // Merged face. Only new0 kept.
                    if (new0 < 0 || new1 == -1)
                    {
                        FatalErrorInFunction
                            << "Problem: oldFaces:" << oldSplit
                            << " newFace:" << labelPair(new0, new1)
                            << exit(FatalError);
                    }
                    changedFaces.append(new0);
                }
            }

            facePairs.setSize(newIndex);
            originalFaces.setSize(newIndex);


            // Update intersections
            updateMesh(map, growFaceCellFace(changedFaces));
        }

        // Update baffle data
        // ~~~~~~~~~~~~~~~~~~
        {
            if (duplicateFace.size())
            {
                meshRefinement::updateList
                (
                    map.faceMap(),
                    label(-1),
                    duplicateFace
                );
            }

            const labelList& reverseFaceMap = map.reverseFaceMap();
            forAll(baffles, i)
            {
                labelPair& baffle = baffles[i];
                baffle.first() = reverseFaceMap[baffle.first()];
                baffle.second() = reverseFaceMap[baffle.second()];

                if (baffle.first() == -1 || baffle.second() == -1)
                {
                    FatalErrorInFunction
                        << "Removed baffle : faces:" << baffle
                        << exit(FatalError);
                }
            }
        }

    }

    return nSplit;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshRefinement::meshRefinement
(
    fvMesh& mesh,
    const scalar mergeDistance,
    const bool overwrite,
    const refinementSurfaces& surfaces,
    const refinementFeatures& features,
    const shellSurfaces& shells,
    const shellSurfaces& limitShells,
    const labelUList& checkFaces,
    const MeshType meshType,
    const bool dryRun
)
:
    mesh_(mesh),
    mergeDistance_(mergeDistance),
    overwrite_(overwrite),
    oldInstance_(mesh.pointsInstance()),
    surfaces_(surfaces),
    features_(features),
    shells_(shells),
    limitShells_(limitShells),
    meshType_(meshType),
    dryRun_(dryRun),
    meshCutter_
    (
        mesh,
        false   // do not try to read history.
    ),
    surfaceIndex_
    (
        IOobject
        (
            "surfaceIndex",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        labelList(mesh_.nFaces(), -1)
    ),
    userFaceData_(0)
{
    // recalculate intersections for all faces
    updateIntersections(checkFaces);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::meshRefinement::surfaceIndex() const
{
    if (surfaceIndex_.size() != mesh_.nFaces())
    {
        const_cast<meshRefinement&>(*this).updateIntersections
        (
            identity(mesh_.nFaces())
        );
    }
    return surfaceIndex_;
}


Foam::labelList& Foam::meshRefinement::surfaceIndex()
{
    if (surfaceIndex_.size() != mesh_.nFaces())
    {
        updateIntersections(identity(mesh_.nFaces()));
    }
    return surfaceIndex_;
}


Foam::label Foam::meshRefinement::countHits() const
{
    // Stats on edges to test. Count proc faces only once.
    bitSet isMasterFace(syncTools::getMasterFaces(mesh_));

    label nHits = 0;

    const labelList& surfIndex = surfaceIndex();

    forAll(surfIndex, facei)
    {
        if (surfIndex[facei] >= 0 && isMasterFace.test(facei))
        {
            ++nHits;
        }
    }
    return nHits;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const bool keepZoneFaces,
    const bool keepBaffles,
    const labelList& singleProcPoints,
    const scalarField& cellWeights,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    autoPtr<mapDistributePolyMesh> map;

    if (Pstream::parRun())
    {
        // Wanted distribution
        labelList distribution;


        // Faces where owner and neighbour are not 'connected' so can
        // go to different processors.
        boolList blockedFace;
        label nUnblocked = 0;

        // Faces that move as block onto single processor
        PtrList<labelList> specifiedProcessorFaces;
        labelList specifiedProcessor;

        // Pairs of baffles
        List<labelPair> couples;

        // Constraints from decomposeParDict
        decomposer.setConstraints
        (
            mesh_,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples
        );


        if (keepZoneFaces || keepBaffles)
        {
            if (keepZoneFaces)
            {
                // Determine decomposition to keep/move surface zones
                // on one processor. The reason is that snapping will make these
                // into baffles, move and convert them back so if they were
                // proc boundaries after baffling&moving the points might be no
                // longer synchronised so recoupling will fail. To prevent this
                // keep owner&neighbour of such a surface zone on the same
                // processor.

                const faceZoneMesh& fZones = mesh_.faceZones();
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

                // Get faces whose owner and neighbour should stay together,
                // i.e. they are not 'blocked'.

                forAll(fZones, zonei)
                {
                    const faceZone& fZone = fZones[zonei];

                    forAll(fZone, i)
                    {
                        label facei = fZone[i];
                        if (blockedFace[facei])
                        {
                            if
                            (
                                mesh_.isInternalFace(facei)
                             || pbm[pbm.whichPatch(facei)].coupled()
                            )
                            {
                                blockedFace[facei] = false;
                                nUnblocked++;
                            }
                        }
                    }
                }


                // If the faceZones are not synchronised the blockedFace
                // might not be synchronised. If you are sure the faceZones
                // are synchronised remove below check.
                syncTools::syncFaceList
                (
                    mesh_,
                    blockedFace,
                    andEqOp<bool>()     // combine operator
                );
            }
            reduce(nUnblocked, sumOp<label>());

            if (keepZoneFaces)
            {
                Info<< "Found " << nUnblocked
                    << " zoned faces to keep together." << endl;
            }


            // Extend un-blockedFaces with any cyclics
            {
                boolList separatedCoupledFace(mesh_.nFaces(), false);
                selectSeparatedCoupledFaces(separatedCoupledFace);

                label nSeparated = 0;
                forAll(separatedCoupledFace, facei)
                {
                    if (separatedCoupledFace[facei])
                    {
                        if (blockedFace[facei])
                        {
                            blockedFace[facei] = false;
                            nSeparated++;
                        }
                    }
                }
                reduce(nSeparated, sumOp<label>());
                Info<< "Found " << nSeparated
                    << " separated coupled faces to keep together." << endl;

                nUnblocked += nSeparated;
            }


            if (keepBaffles)
            {
                const label nBnd = mesh_.nBoundaryFaces();

                labelList coupledFace(mesh_.nFaces(), -1);
                {
                    // Get boundary baffles that need to stay together
                    List<labelPair> allCouples
                    (
                        localPointRegion::findDuplicateFacePairs(mesh_)
                    );

                    // Merge with any couples from
                    // decompositionMethod::setConstraints
                    forAll(couples, i)
                    {
                        const labelPair& baffle = couples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                    forAll(allCouples, i)
                    {
                        const labelPair& baffle = allCouples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                }

                couples.setSize(nBnd);
                label nCpl = 0;
                forAll(coupledFace, facei)
                {
                    if (coupledFace[facei] != -1 && facei < coupledFace[facei])
                    {
                        couples[nCpl++] = labelPair(facei, coupledFace[facei]);
                    }
                }
                couples.setSize(nCpl);
            }
            label nCouples = returnReduce(couples.size(), sumOp<label>());

            if (keepBaffles)
            {
                Info<< "Found " << nCouples << " baffles to keep together."
                    << endl;
            }
        }


        // Make sure blockedFace not set on couples
        forAll(couples, i)
        {
            const labelPair& baffle = couples[i];
            blockedFace[baffle.first()] = false;
            blockedFace[baffle.second()] = false;
        }

        if (returnReduceOr(singleProcPoints.size()))
        {
            // Modify couples to force all selected points be on the same
            // processor

            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

            label nPointFaces = 0;
            for (const label pointi : singleProcPoints)
            {
                for (const label facei : mesh_.pointFaces()[pointi])
                {
                    if (blockedFace[facei])
                    {
                        if
                        (
                            mesh_.isInternalFace(facei)
                         || pbm[pbm.whichPatch(facei)].coupled()
                        )
                        {
                            blockedFace[facei] = false;
                            nPointFaces++;
                        }
                    }
                }
            }
            reduce(nPointFaces, sumOp<label>());
            Info<< "Found " << nPointFaces
                << " additional point-coupled faces to keep together." << endl;

            nUnblocked += nPointFaces;
        }


        distribution = decomposer.decompose
        (
            mesh_,
            cellWeights,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples                 // explicit connections
        );

        if (debug)
        {
            labelList nProcCells(distributor.countCells(distribution));
            Pout<< "Wanted distribution:" << nProcCells << endl;

            Pstream::listCombineReduce(nProcCells, plusEqOp<label>());

            Pout<< "Wanted resulting decomposition:" << endl;
            forAll(nProcCells, proci)
            {
                Pout<< "    " << proci << '\t' << nProcCells[proci] << endl;
            }
            Pout<< endl;
        }

        // Do actual sending/receiving of mesh
        map = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        distribute(map());

        // Set correct instance (for if overwrite)
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());


        if (debug && keepZoneFaces)
        {
            const faceZoneMesh& fZones = mesh_.faceZones();
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

            // Check that faceZone faces are indeed internal
            forAll(fZones, zonei)
            {
                const faceZone& fZone = fZones[zonei];

                forAll(fZone, i)
                {
                    label facei = fZone[i];
                    label patchi = pbm.whichPatch(facei);

                    if (patchi >= 0 && pbm[patchi].coupled())
                    {
                        WarningInFunction
                            << "Face at " << mesh_.faceCentres()[facei]
                            << " on zone " << fZone.name()
                            << " is on coupled patch " << pbm[patchi].name()
                            << endl;
                    }
                }
            }
        }
    }

    return map;
}


Foam::labelList Foam::meshRefinement::intersectedFaces() const
{
    label nBoundaryFaces = 0;

    const labelList& surfIndex = surfaceIndex();

    forAll(surfIndex, facei)
    {
        if (surfIndex[facei] != -1)
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;

    forAll(surfIndex, facei)
    {
        if (surfIndex[facei] != -1)
        {
            surfaceFaces[nBoundaryFaces++] = facei;
        }
    }
    return surfaceFaces;
}


Foam::labelList Foam::meshRefinement::intersectedPoints() const
{
    const faceList& faces = mesh_.faces();

    // Mark all points on faces that will become baffles
    bitSet isBoundaryPoint(mesh_.nPoints());

    const labelList& surfIndex = surfaceIndex();

    forAll(surfIndex, facei)
    {
        if (surfIndex[facei] != -1)
        {
            isBoundaryPoint.set(faces[facei]);
        }
    }

    //// Insert all meshed patches.
    //labelList adaptPatchIDs(meshedPatches());
    //forAll(adaptPatchIDs, i)
    //{
    //    label patchi = adaptPatchIDs[i];
    //
    //    if (patchi != -1)
    //    {
    //        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
    //        forAll(pp, i)
    //        {
    //            isBoundaryPoint.set(faces[pp.start()+i]);
    //        }
    //    }
    //}

    // Make sure all processors have the same data
    syncTools::syncPointList(mesh_, isBoundaryPoint, orEqOp<unsigned int>(), 0);

    return isBoundaryPoint.sortedToc();
}


//- Create patch from set of patches
Foam::autoPtr<Foam::indirectPrimitivePatch> Foam::meshRefinement::makePatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFacei = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFacei++;
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>(mesh.faces(), addressing),
        mesh.points()
    );
}


Foam::tmp<Foam::pointVectorField> Foam::meshRefinement::makeDisplacementField
(
    const pointMesh& pMesh,
    const labelList& adaptPatchIDs
)
{
    const polyMesh& mesh = pMesh();

    // Construct displacement field.
    const pointBoundaryMesh& pointPatches = pMesh.boundary();

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchVectorField::typeName
    );

    forAll(adaptPatchIDs, i)
    {
        patchFieldTypes[adaptPatchIDs[i]] =
            fixedValuePointPatchVectorField::typeName;
    }

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
        }
    }

    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.
    return tmp<pointVectorField>::New
    (
        IOobject
        (
            "pointDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector(dimLength, Zero),
        patchFieldTypes
    );
}


void Foam::meshRefinement::checkCoupledFaceZones(const polyMesh& mesh)
{
    const faceZoneMesh& fZones = mesh.faceZones();

    // Check any zones are present anywhere and in same order

    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = fZones.names();
        Pstream::allGatherList(zoneNames);

        // All have same data now. Check.
        forAll(zoneNames, proci)
        {
            if
            (
                proci != Pstream::myProcNo()
             && (zoneNames[proci] != zoneNames[Pstream::myProcNo()])
            )
            {
                FatalErrorInFunction
                    << "faceZones are not synchronised on processors." << nl
                    << "Processor " << proci << " has faceZones "
                    << zoneNames[proci] << nl
                    << "Processor " << Pstream::myProcNo()
                    << " has faceZones "
                    << zoneNames[Pstream::myProcNo()] << nl
                    << exit(FatalError);
            }
        }
    }

    // Check that coupled faces are present on both sides.

    labelList faceToZone(mesh.nBoundaryFaces(), -1);

    forAll(fZones, zonei)
    {
        const faceZone& fZone = fZones[zonei];

        forAll(fZone, i)
        {
            label bFacei = fZone[i]-mesh.nInternalFaces();

            if (bFacei >= 0)
            {
                if (faceToZone[bFacei] == -1)
                {
                    faceToZone[bFacei] = zonei;
                }
                else if (faceToZone[bFacei] == zonei)
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is twice in zone!"
                        << abort(FatalError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is also in zone "
                        << fZones[faceToZone[bFacei]].name()
                        << abort(FatalError);
                }
            }
        }
    }

    labelList neiFaceToZone(faceToZone);
    syncTools::swapBoundaryFaceList(mesh, neiFaceToZone);

    forAll(faceToZone, i)
    {
        if (faceToZone[i] != neiFaceToZone[i])
        {
            FatalErrorInFunction
                << "Face " << mesh.nInternalFaces()+i
                << " is in zone " << faceToZone[i]
                << ", its coupled face is in zone " << neiFaceToZone[i]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::calculateEdgeWeights
(
    const polyMesh& mesh,
    const bitSet& isMasterEdge,
    const labelList& meshPoints,
    const edgeList& edges,
    scalarField& edgeWeights,
    scalarField& invSumWeight
)
{
    const pointField& pts = mesh.points();

    // Calculate edgeWeights and inverse sum of edge weights
    edgeWeights.setSize(isMasterEdge.size());
    invSumWeight.setSize(meshPoints.size());

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        scalar eMag = max
        (
            SMALL,
            mag
            (
                pts[meshPoints[e[1]]]
              - pts[meshPoints[e[0]]]
            )
        );
        edgeWeights[edgei] = 1.0/eMag;
    }

    // Sum per point all edge weights
    weightedSum
    (
        mesh,
        isMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        scalarField(meshPoints.size(), 1.0),  // data
        invSumWeight
    );

    // Inplace invert
    forAll(invSumWeight, pointi)
    {
        scalar w = invSumWeight[pointi];

        if (w > 0.0)
        {
            invSumWeight[pointi] = 1.0/w;
        }
    }
}


Foam::label Foam::meshRefinement::appendPatch
(
    fvMesh& mesh,
    const label insertPatchi,
    const word& patchName,
    const dictionary& patchDict
)
{
    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // The new location (at the end)
    const label patchi = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.push_back
    (
        polyPatch::New
        (
            patchName,
            patchDict,
            insertPatchi,
            polyPatches
        )
    );

    fvPatches.push_back
    (
        fvPatch::New
        (
            polyPatches[patchi],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addPatchFields<volScalarField>
    (
        mesh,
        fvPatchFieldBase::calculatedType()
    );
    addPatchFields<volVectorField>
    (
        mesh,
        fvPatchFieldBase::calculatedType()
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        fvPatchFieldBase::calculatedType()
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        fvPatchFieldBase::calculatedType()
    );
    addPatchFields<volTensorField>
    (
        mesh,
        fvPatchFieldBase::calculatedType()
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        fvsPatchFieldBase::calculatedType()
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        fvsPatchFieldBase::calculatedType()
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        fvsPatchFieldBase::calculatedType()
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        fvsPatchFieldBase::calculatedType()
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        fvsPatchFieldBase::calculatedType()
    );
    return patchi;
}


Foam::label Foam::meshRefinement::addPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchInfo
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    const label patchi = polyPatches.findPatchID(patchName);
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    label insertPatchi = polyPatches.size();
    label startFacei = mesh.nFaces();

    forAll(polyPatches, patchi)
    {
        const polyPatch& pp = polyPatches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchi = patchi;
            startFacei = pp.start();
            break;
        }
    }

    dictionary patchDict(patchInfo);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFacei);

    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    label addedPatchi = appendPatch(mesh, insertPatchi, patchName, patchDict);


    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(addedPatchi+1);
    for (label i = 0; i < insertPatchi; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchi; i < addedPatchi; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[addedPatchi] = insertPatchi;

    // Shuffle into place
    polyPatches.reorder(oldToNew, true);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    return insertPatchi;
}


Foam::label Foam::meshRefinement::addMeshedPatch
(
    const word& name,
    const dictionary& patchInfo
)
{
    label meshedi = meshedPatches_.find(name);

    if (meshedi != -1)
    {
        // Already there. Get corresponding polypatch
        return mesh_.boundaryMesh().findPatchID(name);
    }
    else
    {
        // Add patch
        label patchi = addPatch(mesh_, name, patchInfo);

//        dictionary patchDict(patchInfo);
//        patchDict.set("nFaces", 0);
//        patchDict.set("startFace", 0);
//        autoPtr<polyPatch> ppPtr
//        (
//            polyPatch::New
//            (
//                name,
//                patchDict,
//                0,
//                mesh_.boundaryMesh()
//            )
//        );
//        label patchi = fvMeshTools::addPatch
//        (
//            mesh_,
//            ppPtr(),
//            dictionary(),       // optional field values
//            fvPatchFieldBase::calculatedType(),
//            true
//        );

        // Store
        meshedPatches_.append(name);

        // Clear patch based addressing
        faceToCoupledPatch_.clear();

        return patchi;
    }
}


Foam::labelList Foam::meshRefinement::meshedPatches() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<label> patchIDs(meshedPatches_.size());
    forAll(meshedPatches_, i)
    {
        label patchi = patches.findPatchID(meshedPatches_[i]);

        if (patchi == -1)
        {
            WarningInFunction
                << "Problem : did not find patch " << meshedPatches_[i]
                << endl << "Valid patches are " << patches.names()
                << endl;    //abort(FatalError);
        }
        else if (!polyPatch::constraintType(patches[patchi].type()))
        {
            patchIDs.append(patchi);
        }
    }

    return patchIDs;
}


Foam::label Foam::meshRefinement::addFaceZone
(
    const word& fzName,
    const word& masterPatch,
    const word& slavePatch,
    const surfaceZonesInfo::faceZoneType& fzType
)
{
    label zonei = surfaceZonesInfo::addFaceZone
    (
        fzName,   //name
        labelList(0),   //addressing
        boolList(0),    //flipmap
        mesh_
    );

    faceZoneToMasterPatch_.insert(fzName, masterPatch);
    faceZoneToSlavePatch_.insert(fzName, slavePatch);
    faceZoneToType_.insert(fzName, fzType);

    return zonei;
}


bool Foam::meshRefinement::getFaceZoneInfo
(
    const word& fzName,
    label& masterPatchID,
    label& slavePatchID,
    surfaceZonesInfo::faceZoneType& fzType
) const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!faceZoneToMasterPatch_.found(fzName))
    {
        return false;
    }
    else
    {
        const word& masterName = faceZoneToMasterPatch_[fzName];
        masterPatchID = pbm.findPatchID(masterName);

        const word& slaveName = faceZoneToSlavePatch_[fzName];
        slavePatchID = pbm.findPatchID(slaveName);

        fzType = faceZoneToType_[fzName];

        return true;
    }
}


Foam::label Foam::meshRefinement::addPointZone(const word& name)
{
    pointZoneMesh& pointZones = mesh_.pointZones();

    label zoneI = pointZones.findZoneID(name);

    if (zoneI == -1)
    {
        zoneI = pointZones.size();
        pointZones.clearAddressing();

        pointZones.emplace_back
        (
            name,           // name
            zoneI,          // index
            pointZones      // pointZoneMesh
        );
    }
    return zoneI;
}


void Foam::meshRefinement::selectSeparatedCoupledFaces(boolList& selected) const
{
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMI.
        const auto* cpp = isA<coupledPolyPatch>(pp);

        if (cpp && (cpp->separated() || !cpp->parallel()))
        {
            SubList<bool>(selected, pp.size(), pp.start()) = true;
        }
    }
}


Foam::labelList Foam::meshRefinement::countEdgeFaces
(
    const uindirectPrimitivePatch& pp
) const
{
    // Count number of faces per edge. Parallel consistent.

    const labelListList& edgeFaces = pp.edgeFaces();
    labelList nEdgeFaces(edgeFaces.size());
    forAll(edgeFaces, edgei)
    {
        nEdgeFaces[edgei] = edgeFaces[edgei].size();
    }

    // Sync across processor patches
    if (Pstream::parRun())
    {
        const globalMeshData& globalData = mesh_.globalData();
        const mapDistribute& map = globalData.globalEdgeSlavesMap();
        const indirectPrimitivePatch& cpp = globalData.coupledPatch();

        // Match pp edges to coupled edges
        labelList patchEdges;
        labelList coupledEdges;
        PackedBoolList sameEdgeOrientation;
        PatchTools::matchEdges
        (
            pp,
            cpp,
            patchEdges,
            coupledEdges,
            sameEdgeOrientation
        );

        // Convert patch-edge data into cpp-edge data
        labelList coupledNEdgeFaces(map.constructSize(), Zero);
        UIndirectList<label>(coupledNEdgeFaces, coupledEdges) =
            UIndirectList<label>(nEdgeFaces, patchEdges);

        // Synchronise
        globalData.syncData
        (
            coupledNEdgeFaces,
            globalData.globalEdgeSlaves(),
            globalData.globalEdgeTransformedSlaves(),
            map,
            plusEqOp<label>()
        );

        // Convert back from cpp-edge to patch-edge
        UIndirectList<label>(nEdgeFaces, patchEdges) =
            UIndirectList<label>(coupledNEdgeFaces, coupledEdges);
    }
    return nEdgeFaces;
}


Foam::label Foam::meshRefinement::findCell
(
    const polyMesh& mesh,
    const vector& perturbVec,
    const point& p
)
{
    // Force calculation of base points (needs to be synchronised)
    (void)mesh.tetBasePtIs();

    label celli = mesh.findCell(p, findCellMode);

    if (returnReduceAnd(celli < 0))
    {
        // See if we can perturb a bit
        celli = mesh.findCell(p+perturbVec, findCellMode);
    }
    return celli;
}


Foam::label Foam::meshRefinement::findRegion
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const vector& perturbVec,
    const point& p
)
{
    label regioni = -1;

    // Force calculation of base points (needs to be synchronised)
    (void)mesh.tetBasePtIs();

    label celli = mesh.findCell(p, findCellMode);
    if (celli != -1)
    {
        regioni = cellToRegion[celli];
    }
    reduce(regioni, maxOp<label>());

    if (regioni == -1)
    {
        // See if we can perturb a bit
        celli = mesh.findCell(p+perturbVec, findCellMode);
        if (celli != -1)
        {
            regioni = cellToRegion[celli];
        }
        reduce(regioni, maxOp<label>());
    }
    return regioni;
}


Foam::fileName Foam::meshRefinement::writeLeakPath
(
    const polyMesh& mesh,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const boolList& blockedFace,
    coordSetWriter& writer
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    fileName outputDir
    (
        mesh.time().globalPath()
      / functionObject::outputPrefix
      / mesh.pointsInstance()
    );
    outputDir.clean();  // Remove unneeded ".."

    // Write the leak path

    meshSearch searchEngine(mesh);
    shortestPathSet leakPath
    (
        "leakPath",
        mesh,
        searchEngine,
        coordSet::coordFormatNames[coordSet::coordFormat::DISTANCE],
        false,  //true,
        50,     // tbd. Number of iterations
        pbm.groupPatchIDs()["wall"],
        locationsInMesh,
        locationsOutsideMesh,
        blockedFace
    );

    // Split leak path according to segment. Note: segment index
    // is global (= index in locationsInsideMesh)
    List<pointList> segmentPoints;
    List<scalarList> segmentDist;
    {
        label nSegments = 0;
        if (leakPath.segments().size())
        {
            nSegments = max(leakPath.segments())+1;
        }
        reduce(nSegments, maxOp<label>());

        labelList nElemsPerSegment(nSegments, Zero);
        for (label segmenti : leakPath.segments())
        {
            nElemsPerSegment[segmenti]++;
        }
        segmentPoints.setSize(nElemsPerSegment.size());
        segmentDist.setSize(nElemsPerSegment.size());
        forAll(nElemsPerSegment, i)
        {
            segmentPoints[i].setSize(nElemsPerSegment[i]);
            segmentDist[i].setSize(nElemsPerSegment[i]);
        }
        nElemsPerSegment = 0;

        forAll(leakPath, elemi)
        {
            label segmenti = leakPath.segments()[elemi];
            pointList& points = segmentPoints[segmenti];
            scalarList& dist = segmentDist[segmenti];
            label& n = nElemsPerSegment[segmenti];

            points[n] = leakPath[elemi];
            dist[n] = leakPath.distance()[elemi];
            n++;
        }
    }

    PtrList<coordSet> allLeakPaths(segmentPoints.size());
    forAll(allLeakPaths, segmenti)
    {
        // Collect data from all processors
        List<pointList> gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = std::move(segmentPoints[segmenti]);
        Pstream::gatherList(gatheredPts);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = std::move(segmentDist[segmenti]);
        Pstream::gatherList(gatheredDist);

        // Combine processor lists into one big list.
        pointList allPts
        (
            ListListOps::combine<pointList>
            (
                gatheredPts, accessOp<pointList>()
            )
        );
        scalarList allDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );

        // Sort according to distance
        labelList indexSet(Foam::sortedOrder(allDist));

        allLeakPaths.set
        (
            segmenti,
            new coordSet
            (
                leakPath.name(),
                leakPath.axis(),
                pointList(allPts, indexSet),
                //scalarList(allDist, indexSet)
                scalarList(allPts.size(), scalar(segmenti))
            )
        );
    }

    fileName fName;
    if (Pstream::master())
    {
        List<scalarField> allLeakData(allLeakPaths.size());
        forAll(allLeakPaths, segmenti)
        {
            allLeakData[segmenti] = allLeakPaths[segmenti].distance();
        }

        writer.nFields(1);

        writer.open
        (
            allLeakPaths,
            (outputDir / allLeakPaths[0].name())
        );

        fName = writer.write("leakPath", allLeakData);

        // Force writing to finish before FatalError exit
        writer.close(true);
    }

    // Probably do not need to broadcast name (only written on master anyhow)
    Pstream::broadcast(fName);

    return fName;
}


// Modify cellRegion to be consistent with locationsInMesh.
// - all regions not in locationsInMesh are set to -1
// - check that all regions inside locationsOutsideMesh are set to -1
Foam::label Foam::meshRefinement::findRegions
(
    const polyMesh& mesh,
    const vector& perturbVec,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const label nRegions,
    labelList& cellRegion,
    const boolList& blockedFace,
    // Leak-path
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter
)
{
    bitSet insideCell(mesh.nCells());

    // Mark all cells reachable from locationsInMesh
    labelList insideRegions(locationsInMesh.size());
    forAll(insideRegions, i)
    {
        // Find the region containing the point
        label regioni = findRegion
        (
            mesh,
            cellRegion,
            perturbVec,
            locationsInMesh[i]
        );

        insideRegions[i] = regioni;

        // Mark all cells in the region as being inside
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == regioni)
            {
                insideCell.set(celli);
            }
        }
    }



    // Check that all the locations outside the
    // mesh do not conflict with those inside
    forAll(locationsOutsideMesh, i)
    {
        // Find the region containing the point,
        // and the corresponding inside region index

        label indexi;
        label regioni = findRegion
        (
            mesh,
            cellRegion,
            perturbVec,
            locationsOutsideMesh[i]
        );

        if (regioni == -1 && (indexi = insideRegions.find(regioni)) != -1)
        {
            if (leakPathFormatter)
            {
                const fileName fName
                (
                    writeLeakPath
                    (
                        mesh,
                        locationsInMesh,
                        locationsOutsideMesh,
                        blockedFace,
                        leakPathFormatter.constCast()
                    )
                );
                Info<< "Dumped leak path to " << fName << endl;
            }

            auto& err =
            (
                exitIfLeakPath
              ? FatalErrorInFunction
              : WarningInFunction
            );

            err << "Location in mesh " << locationsInMesh[indexi]
                << " is inside same mesh region " << regioni
                << " as one of the locations outside mesh "
                << locationsOutsideMesh << endl;

            if (exitIfLeakPath)
            {
                FatalError << exit(FatalError);
            }
        }
    }


    label nRemove = 0;

    // Now update cellRegion to -1 for unreachable cells
    forAll(insideCell, celli)
    {
        if (!insideCell.test(celli))
        {
            cellRegion[celli] = -1;
            ++nRemove;
        }
        else if (cellRegion[celli] == -1)
        {
            ++nRemove;
        }
    }

    return nRemove;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMeshRegions
(
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const bool exitIfLeakPath,
    const refPtr<coordSetWriter>& leakPathFormatter
)
{
    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.

    boolList blockedFace(mesh_.nFaces(), false);
    selectSeparatedCoupledFaces(blockedFace);

    regionSplit cellRegion(mesh_, blockedFace);

    label nRemove = findRegions
    (
        mesh_,
        vector::uniform(mergeDistance_),   // perturbVec
        locationsInMesh,
        locationsOutsideMesh,
        cellRegion.nRegions(),
        cellRegion,
        blockedFace,
        // Leak-path
        exitIfLeakPath,
        leakPathFormatter
    );

    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(nRemove);
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] == -1)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nTotCellsToRemove = returnReduce
    (
        cellsToRemove.size(),
        sumOp<label>()
    );


    autoPtr<mapPolyMesh> mapPtr;
    if (nTotCellsToRemove > 0)
    {
        label nCellsToKeep =
            mesh_.globalData().nTotalCells()
          - nTotCellsToRemove;

        Info<< "Keeping all cells containing points " << locationsInMesh << endl
            << "Selected for keeping : "
            << nCellsToKeep
            << " cells." << endl;


        // Remove cells
        removeCells cellRemover(mesh_);

        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatch;

        label nExposedFaces = returnReduce(exposedFaces.size(), sumOp<label>());
        if (nExposedFaces)
        {
            // FatalErrorInFunction
            //    << "Removing non-reachable cells should only expose"
            //    << " boundary faces" << nl
            //    << "ExposedFaces:" << exposedFaces << abort(FatalError);

            // Patch for exposed faces for lack of anything sensible.
            label defaultPatch = 0;
            if (globalToMasterPatch.size())
            {
                defaultPatch = globalToMasterPatch[0];
            }

            WarningInFunction
                << "Removing non-reachable cells exposes "
                << nExposedFaces << " internal or coupled faces." << endl
                << "    These get put into patch " << defaultPatch << endl;
            exposedPatch.setSize(exposedFaces.size(), defaultPatch);
        }

        mapPtr = doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatch,
            cellRemover
        );
    }
    return mapPtr;
}


void Foam::meshRefinement::distribute(const mapDistributePolyMesh& map)
{
    // mesh_ already distributed; distribute my member data

    // surfaceQueries_ ok.

    // refinement
    meshCutter_.distribute(map);

    // surfaceIndex is face data.
    map.distributeFaceData(surfaceIndex_);

    // faceToPatch (baffles that were on coupled faces) is not maintained
    // (since baffling also disconnects points)
    faceToCoupledPatch_.clear();

    // maintainedFaces are indices of faces.
    forAll(userFaceData_, i)
    {
        map.distributeFaceData(userFaceData_[i].second());
    }

    // Redistribute surface and any fields on it.
    {
        Random rndGen(653213);

        // Get local mesh bounding box. Single box for now.
        List<treeBoundBox> meshBb
        (
            1,
            treeBoundBox(mesh_.points()).extend(rndGen, 1e-4)
        );

        // Distribute all geometry (so refinementSurfaces and shellSurfaces)
        searchableSurfaces& geometry =
            const_cast<searchableSurfaces&>(surfaces_.geometry());

        forAll(geometry, i)
        {
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            geometry[i].distribute
            (
                meshBb,
                false,          // do not keep outside triangles
                faceMap,
                pointMap
            );

            if (faceMap)
            {
                // (ab)use the instance() to signal current modification time
                geometry[i].instance() = geometry[i].time().timeName();
            }

            faceMap.clear();
            pointMap.clear();
        }
    }
}


void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces
)
{
    Map<label> dummyMap(0);

    updateMesh(map, changedFaces, dummyMap, dummyMap, dummyMap);
}


void Foam::meshRefinement::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    // For now only meshCutter has storable/retrievable data.
    meshCutter_.storeData
    (
        pointsToStore,
        facesToStore,
        cellsToStore
    );
}


void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // For now only meshCutter has storable/retrievable data.

    // Update numbering of cells/vertices.
    meshCutter_.updateMesh
    (
        map,
        pointsToRestore,
        facesToRestore,
        cellsToRestore
    );

    // Update surfaceIndex
    updateList(map.faceMap(), label(-1), surfaceIndex_);

    // Update faceToCoupledPatch_
    {
        Map<label> newFaceToPatch(faceToCoupledPatch_.size());
        forAllConstIters(faceToCoupledPatch_, iter)
        {
            const label newFacei = map.reverseFaceMap()[iter.key()];

            if (newFacei >= 0)
            {
                newFaceToPatch.insert(newFacei, iter.val());
            }
        }
        faceToCoupledPatch_.transfer(newFaceToPatch);
    }


    // Update cached intersection information
    updateIntersections(changedFaces);

    // Update maintained faces
    forAll(userFaceData_, i)
    {
        labelList& data = userFaceData_[i].second();

        if (userFaceData_[i].first() == KEEPALL)
        {
            // extend list with face-from-face data
            updateList(map.faceMap(), label(-1), data);
        }
        else if (userFaceData_[i].first() == MASTERONLY)
        {
            // keep master only
            labelList newFaceData(map.faceMap().size(), -1);

            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0 && map.reverseFaceMap()[oldFacei] == facei)
                {
                    newFaceData[facei] = data[oldFacei];
                }
            }
            data.transfer(newFaceData);
        }
        else
        {
            // remove any face that has been refined i.e. referenced more than
            // once.

            // 1. Determine all old faces that get referenced more than once.
            // These get marked with -1 in reverseFaceMap
            labelList reverseFaceMap(map.reverseFaceMap());

            forAll(map.faceMap(), facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] != facei)
                    {
                        // facei is slave face. Mark old face.
                        reverseFaceMap[oldFacei] = -1;
                    }
                }
            }

            // 2. Map only faces with intact reverseFaceMap
            labelList newFaceData(map.faceMap().size(), -1);
            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] == facei)
                    {
                        newFaceData[facei] = data[oldFacei];
                    }
                }
            }
            data.transfer(newFaceData);
        }
    }
}


bool Foam::meshRefinement::write() const
{
    bool writeOk = mesh_.write();

    // Make sure that any distributed surfaces (so ones which probably have
    // been changed) get written as well.
    // Note: should ideally have some 'modified' flag to say whether it
    // has been changed or not.
    searchableSurfaces& geometry =
        const_cast<searchableSurfaces&>(surfaces_.geometry());

    forAll(geometry, i)
    {
        searchableSurface& s = geometry[i];

        // Check if instance() of surface is not constant or system.
        // Is good hint that surface is distributed.
        if
        (
            s.instance() != s.time().system()
         && s.instance() != s.time().caseSystem()
         && s.instance() != s.time().constant()
         && s.instance() != s.time().caseConstant()
        )
        {
            // Make sure it gets written to current time, not constant.
            s.instance() = s.time().timeName();
            writeOk = writeOk && s.write();
        }
    }

    return writeOk;
}


Foam::bitSet Foam::meshRefinement::getMasterPoints
(
    const polyMesh& mesh,
    const labelList& meshPoints
)
{
    const label myProci = UPstream::myProcNo();
    const globalIndex globalPoints(meshPoints.size());

    labelList myPoints
    (
        Foam::identity(globalPoints.range(myProci))
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        myPoints,
        minEqOp<label>(),
        labelMax
    );


    bitSet isPatchMasterPoint(meshPoints.size());
    forAll(meshPoints, pointi)
    {
        if (myPoints[pointi] == globalPoints.toGlobal(myProci, pointi))
        {
            isPatchMasterPoint.set(pointi);
        }
    }

    return isPatchMasterPoint;
}


Foam::bitSet Foam::meshRefinement::getMasterEdges
(
    const polyMesh& mesh,
    const labelList& meshEdges
)
{
    const label myProci = UPstream::myProcNo();
    const globalIndex globalEdges(meshEdges.size());

    labelList myEdges
    (
        Foam::identity(globalEdges.range(myProci))
    );

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        myEdges,
        minEqOp<label>(),
        labelMax
    );


    bitSet isMasterEdge(meshEdges.size());
    forAll(meshEdges, edgei)
    {
        if (myEdges[edgei] == globalEdges.toGlobal(myProci, edgei))
        {
            isMasterEdge.set(edgei);
        }
    }

    return isMasterEdge;
}


void Foam::meshRefinement::printMeshInfo
(
    const bool debug,
    const string& msg,
    const bool printCellLevel
)
const
{
    const globalMeshData& pData = mesh_.globalData();

    if (debug)
    {
        Pout<< msg.c_str()
            << " : cells(local):" << mesh_.nCells()
            << "  faces(local):" << mesh_.nFaces()
            << "  points(local):" << mesh_.nPoints()
            << endl;
    }

    {
        bitSet isMasterFace(syncTools::getMasterFaces(mesh_));
        label nMasterFaces = isMasterFace.count();

        bitSet isMeshMasterPoint(syncTools::getMasterPoints(mesh_));
        label nMasterPoints = isMeshMasterPoint.count();

        Info<< msg.c_str()
            << " : cells:" << pData.nTotalCells()
            << "  faces:" << returnReduce(nMasterFaces, sumOp<label>())
            << "  points:" << returnReduce(nMasterPoints, sumOp<label>());

        if (UPstream::parRun())
        {
            const scalar nIdealCells =
                scalar(pData.nTotalCells())/Pstream::nProcs();
            const scalar unbalance = returnReduce
            (
                mag(1.0-mesh_.nCells()/nIdealCells),
                maxOp<scalar>()
            );
            Info<< "  unbalance:" << unbalance;
        }
        Info<< endl;
    }

    if (printCellLevel)
    {
        const labelList& cellLevel = meshCutter_.cellLevel();

        labelList nCells(gMax(cellLevel)+1, Zero);

        forAll(cellLevel, celli)
        {
            nCells[cellLevel[celli]]++;
        }

        Pstream::listCombineGather(nCells, plusEqOp<label>());

        /// Pstream::broadcast(nCells);
        if (Pstream::master())
        {
            Info<< "Cells per refinement level:" << endl;
            forAll(nCells, leveli)
            {
                Info<< "    " << leveli << '\t' << nCells[leveli]
                    << endl;
            }
        }
    }
}


Foam::word Foam::meshRefinement::timeName() const
{
    if (overwrite_ && mesh_.time().timeIndex() == 0)
    {
        return oldInstance_;
    }

    return mesh_.time().timeName();
}


void Foam::meshRefinement::dumpRefinementLevel() const
{
    // Note: use time().timeName(), not meshRefinement::timeName()
    // so as to dump the fields to 0, not to constant.
    {
        volScalarField volRefLevel
        (
            IOobject
            (
                "cellLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(volRefLevel, celli)
        {
            volRefLevel[celli] = cellLevel[celli];
        }

        volRefLevel.write();
    }

    // Dump pointLevel
    {
        const pointMesh& pMesh = pointMesh::New(mesh_);

        pointScalarField pointRefLevel
        (
            IOobject
            (
                "pointLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            pMesh,
            dimensionedScalar(dimless, Zero)
        );

        const labelList& pointLevel = meshCutter_.pointLevel();

        forAll(pointRefLevel, pointi)
        {
            pointRefLevel[pointi] = pointLevel[pointi];
        }

        pointRefLevel.write();
    }
}


void Foam::meshRefinement::dumpIntersections(const fileName& prefix) const
{
    {
        OFstream str(prefix + "_edges.obj");
        label verti = 0;
        Pout<< "meshRefinement::dumpIntersections :"
            << " Writing cellcentre-cellcentre intersections to file "
            << str.name() << endl;


        // Redo all intersections
        // ~~~~~~~~~~~~~~~~~~~~~~

        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(mesh_.nBoundaryFaces());
        pointField neiCc(mesh_.nBoundaryFaces());
        calcNeighbourData(neiLevel, neiCc);

        labelList intersectionFaces(intersectedFaces());

        // Collect segments we want to test for
        pointField start(intersectionFaces.size());
        pointField end(intersectionFaces.size());
        {
            labelList minLevel;
            calcCellCellRays
            (
                neiCc,
                labelList(neiCc.size(), -1),
                intersectionFaces,
                start,
                end,
                minLevel
            );
        }


        // Do tests in one go
        labelList surfaceHit;
        List<pointIndexHit> surfaceHitInfo;
        surfaces_.findAnyIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceHitInfo
        );

        forAll(intersectionFaces, i)
        {
            if (surfaceHit[i] != -1)
            {
                meshTools::writeOBJ(str, start[i]);
                verti++;
                meshTools::writeOBJ(str, surfaceHitInfo[i].point());
                verti++;
                meshTools::writeOBJ(str, end[i]);
                verti++;
                str << "l " << verti-2 << ' ' << verti-1 << nl
                    << "l " << verti-1 << ' ' << verti << nl;
            }
        }
    }

    Pout<< endl;
}


void Foam::meshRefinement::write
(
    const debugType debugFlags,
    const writeType writeFlags,
    const fileName& prefix
) const
{
    if (writeFlags & WRITEMESH)
    {
        write();
    }

    if (writeFlags && !(writeFlags & NOWRITEREFINEMENT))
    {
        meshCutter_.write();

        // Force calculation before writing
        (void)surfaceIndex();
        surfaceIndex_.write();
    }

    if (writeFlags & WRITELEVELS)
    {
        dumpRefinementLevel();
    }

    if ((debugFlags & OBJINTERSECTIONS) && prefix.size())
    {
        dumpIntersections(prefix);
    }
}


void Foam::meshRefinement::removeFiles(const polyMesh& mesh)
{
    IOobject io
    (
        "dummy",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh
    );
    fileName setsDir(io.path());

    if (topoSet::debug) DebugVar(setsDir);

    // Remove local files
    if (exists(setsDir/"surfaceIndex"))
    {
        rm(setsDir/"surfaceIndex");
    }

    // Remove other files
    hexRef8::removeFiles(mesh);
}


Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel()
{
    return writeLevel_;
}


void Foam::meshRefinement::writeLevel(const writeType flags)
{
    writeLevel_ = flags;
}


//Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel()
//{
//    return outputLevel_;
//}
//
//
//void Foam::meshRefinement::outputLevel(const outputType flags)
//{
//    outputLevel_ = flags;
//}


const Foam::dictionary& Foam::meshRefinement::subDict
(
    const dictionary& dict,
    const word& keyword,
    const bool noExit,
    enum keyType::option matchOpt
)
{
    const dictionary* dictptr = dict.findDict(keyword, matchOpt);

    if (!dictptr)
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << keyword
            << "' not found (or not a dictionary) in dictionary "
            << dict.relativeName() << nl;

        if (noExit)
        {
            // Dummy return
            return dictionary::null;
        }
        else
        {
            FatalIOError << exit(FatalIOError);
        }
    }

    return *dictptr;
}


Foam::ITstream& Foam::meshRefinement::lookup
(
    const dictionary& dict,
    const word& keyword,
    const bool noExit,
    enum keyType::option matchOpt
)
{
    const entry* eptr = dict.findEntry(keyword, matchOpt);

    if (!eptr)
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << keyword << "' not found in dictionary "
            << dict.relativeName() << nl;

        if (noExit)
        {
            // Dummy return
            return ITstream::empty_stream();
        }
        else
        {
            FatalIOError << exit(FatalIOError);
        }
    }

    return eptr->stream();
}


// ************************************************************************* //
