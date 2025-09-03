/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
#include "trackedParticle.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "faceSet.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyTopoChange.H"
#include "mapDistributePolyMesh.H"
#include "Cloud.H"
//#include "globalIndex.H"
#include "cellSet.H"
#include "treeDataCell.H"

#include "cellCuts.H"
#include "refineCell.H"
#include "hexCellLooper.H"
#include "meshCutter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- To compare normals
    static bool less(const vector& x, const vector& y)
    {
        for (direction i = 0; i < vector::nComponents; i++)
        {
            if (x[i] < y[i])
            {
                return true;
            }
            else if (x[i] > y[i])
            {
                return false;
            }
        }
        // All components the same
        return false;
    }


    //- To compare normals
    class normalLess
    {
        const vectorList& values_;

    public:

        normalLess(const vectorList& values)
        :
            values_(values)
        {}

        bool operator()(const label a, const label b) const
        {
            return less(values_[a], values_[b]);
        }
    };


    //- Template specialization for pTraits<labelList> so we can have fields
    template<>
    class pTraits<labelList>
    {

    public:

        //- Component type
        typedef labelList cmptType;
    };

    //- Template specialization for pTraits<labelList> so we can have fields
    template<>
    class pTraits<vectorList>
    {

    public:

        //- Component type
        typedef vectorList cmptType;
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get faces (on the new mesh) that have in some way been affected by the
// mesh change. Picks up all faces but those that are between two
// unrefined faces. (Note that of an unchanged face the edge still might be
// split but that does not change any face centre or cell centre.
Foam::labelList Foam::meshRefinement::getChangedFaces
(
    const mapPolyMesh& map,
    const labelList& oldCellsToRefine
)
{
    const polyMesh& mesh = map.mesh();

    labelList changedFaces;

    // For reporting: number of masterFaces changed
    label nMasterChanged = 0;
    bitSet isMasterFace(syncTools::getMasterFaces(mesh));

    {
        // Mark any face on a cell which has been added or changed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that refining a face changes the face centre (for a warped face)
        // which changes the cell centre. This again changes the cellcentre-
        // cellcentre edge across all faces of the cell.
        // Note: this does not happen for unwarped faces but unfortunately
        // we don't have this information.

        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const cellList& cells = mesh.cells();
        const label nInternalFaces = mesh.nInternalFaces();

        // Mark refined cells on old mesh
        bitSet oldRefineCell(map.nOldCells(), oldCellsToRefine);

        // Mark refined faces
        bitSet refinedInternalFace(nInternalFaces);

        // 1. Internal faces

        for (label faceI = 0; faceI < nInternalFaces; faceI++)
        {
            label oldOwn = map.cellMap()[faceOwner[faceI]];
            label oldNei = map.cellMap()[faceNeighbour[faceI]];

            if
            (
                oldOwn >= 0 && !oldRefineCell.test(oldOwn)
             && oldNei >= 0 && !oldRefineCell.test(oldNei)
            )
            {
                // Unaffected face since both neighbours were not refined.
            }
            else
            {
                refinedInternalFace.set(faceI);
            }
        }


        // 2. Boundary faces

        boolList refinedBoundaryFace(mesh.nBoundaryFaces(), false);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            label faceI = pp.start();

            forAll(pp, i)
            {
                label oldOwn = map.cellMap()[faceOwner[faceI]];

                if (oldOwn >= 0 && !oldRefineCell.test(oldOwn))
                {
                    // owner did exist and wasn't refined.
                }
                else
                {
                    refinedBoundaryFace[faceI-nInternalFaces] = true;
                }
                faceI++;
            }
        }

        // Synchronise refined face status
        syncTools::syncBoundaryFaceList
        (
            mesh,
            refinedBoundaryFace,
            orEqOp<bool>()
        );


        // 3. Mark all faces affected by refinement. Refined faces are in
        //    - refinedInternalFace
        //    - refinedBoundaryFace
        boolList changedFace(mesh.nFaces(), false);

        forAll(refinedInternalFace, faceI)
        {
            if (refinedInternalFace.test(faceI))
            {
                const cell& ownFaces = cells[faceOwner[faceI]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
                const cell& neiFaces = cells[faceNeighbour[faceI]];
                forAll(neiFaces, neiI)
                {
                    changedFace[neiFaces[neiI]] = true;
                }
            }
        }

        forAll(refinedBoundaryFace, i)
        {
            if (refinedBoundaryFace[i])
            {
                const cell& ownFaces = cells[faceOwner[i+nInternalFaces]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            changedFace,
            orEqOp<bool>()
        );


        // Now we have in changedFace marked all affected faces. Pack.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        changedFaces = findIndices(changedFace, true);

        // Count changed master faces.
        nMasterChanged = 0;

        forAll(changedFace, faceI)
        {
            if (changedFace[faceI] && isMasterFace.test(faceI))
            {
                nMasterChanged++;
            }
        }

    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "getChangedFaces : Detected "
            << " local:" << changedFaces.size()
            << " global:" << returnReduce(nMasterChanged, sumOp<label>())
            << " changed faces out of " << mesh.globalData().nTotalFaces()
            << endl;

        faceSet changedFacesSet(mesh, "changedFaces", changedFaces);
        Pout<< "getChangedFaces : Writing " << changedFaces.size()
            << " changed faces to faceSet " << changedFacesSet.name()
            << endl;
        changedFacesSet.write();
    }

    return changedFaces;
}


// Mark cell for refinement (if not already marked). Return false if
// refinelimit hit. Keeps running count (in nRefine) of cells marked for
// refinement
bool Foam::meshRefinement::markForRefine
(
    const label markValue,
    const label nAllowRefine,

    label& cellValue,
    label& nRefine
)
{
    if (cellValue == -1)
    {
        cellValue = markValue;
        nRefine++;
    }

    return nRefine <= nAllowRefine;
}


void Foam::meshRefinement::markFeatureCellLevel
(
    const pointField& keepPoints,
    labelList& maxFeatureLevel
) const
{
    // We want to refine all cells containing a feature edge.
    // - don't want to search over whole mesh
    // - don't want to build octree for whole mesh
    // - so use tracking to follow the feature edges
    //
    // 1. find non-manifold points on feature edges (i.e. start of feature edge
    //    or 'knot')
    // 2. seed particle starting at keepPoint going to this non-manifold point
    // 3. track particles to their non-manifold point
    //
    // 4. track particles across their connected feature edges, marking all
    //    visited cells with their level (through trackingData)
    // 5. do 4 until all edges have been visited.


    // Find all start cells of features. Is done by tracking from keepPoint.
    Cloud<trackedParticle> startPointCloud
    (
        mesh_,
        Foam::zero{},
        "startPointCloud"
    );


    // Features are identical on all processors. Number them so we know
    // what to seed. Do this on only the processor that
    // holds the keepPoint.

    for (const point& keepPoint : keepPoints)
    {
        const label celli = mesh_.cellTree().findInside(keepPoint);

        if (celli != -1)
        {
            // I am the processor that holds the keepPoint

            forAll(features_, feati)
            {
                const edgeMesh& featureMesh = features_[feati];
                const label featureLevel = features_.levels()[feati][0];
                const labelListList& pointEdges = featureMesh.pointEdges();

                // Find regions on edgeMesh
                labelList edgeRegion;
                label nRegions = featureMesh.regions(edgeRegion);


                bitSet regionVisited(nRegions);


                // 1. Seed all 'knots' in edgeMesh


                forAll(pointEdges, pointi)
                {
                    if (pointEdges[pointi].size() != 2)
                    {
                        if (debug&meshRefinement::FEATURESEEDS)
                        {
                            Pout<< "Adding particle from point:" << pointi
                                << " coord:" << featureMesh.points()[pointi]
                                << " since number of emanating edges:"
                                << pointEdges[pointi].size()
                                << endl;
                        }

                        // Non-manifold point. Create particle.
                        startPointCloud.addParticle
                        (
                            new trackedParticle
                            (
                                mesh_,
                                keepPoint,
                                celli,
                                featureMesh.points()[pointi],   // endpos
                                featureLevel,                   // level
                                feati,                          // featureMesh
                                pointi,                         // end point
                                -1                              // feature edge
                            )
                        );

                        // Mark
                        if (pointEdges[pointi].size() > 0)
                        {
                            label e0 = pointEdges[pointi][0];
                            label regioni = edgeRegion[e0];
                            regionVisited.set(regioni);
                        }
                    }
                }


                // 2. Any regions that have not been visited at all? These can
                //    only be circular regions!
                forAll(featureMesh.edges(), edgei)
                {
                    if (regionVisited.set(edgeRegion[edgei]))
                    {
                        const edge& e = featureMesh.edges()[edgei];
                        label pointi = e.start();
                        if (debug&meshRefinement::FEATURESEEDS)
                        {
                            Pout<< "Adding particle from point:" << pointi
                                << " coord:" << featureMesh.points()[pointi]
                                << " on circular region:" << edgeRegion[edgei]
                                << endl;
                        }

                        // Non-manifold point. Create particle.
                        startPointCloud.addParticle
                        (
                            new trackedParticle
                            (
                                mesh_,
                                keepPoint,
                                celli,
                                featureMesh.points()[pointi],   // endpos
                                featureLevel,                   // level
                                feati,                          // featureMesh
                                pointi,                         // end point
                                -1                              // feature edge
                            )
                        );
                    }
                }
            }
        }
    }


    // Largest refinement level of any feature passed through
    maxFeatureLevel = labelList(mesh_.nCells(), -1);

    // Whether edge has been visited.
    List<bitSet> featureEdgeVisited(features_.size());

    forAll(features_, featI)
    {
        featureEdgeVisited[featI].setSize(features_[featI].edges().size());
        featureEdgeVisited[featI] = false;
    }

    // Database to pass into trackedParticle::move
    trackedParticle::trackingData td
    (
        startPointCloud,
        maxFeatureLevel,
        featureEdgeVisited
    );


    // Track all particles to their end position (= starting feature point)
    // Note that the particle might have started on a different processor
    // so this will transfer across nicely until we can start tracking proper.
    scalar maxTrackLen = 2.0*mesh_.bounds().mag();

    if (debug&meshRefinement::FEATURESEEDS)
    {
        Pout<< "Tracking " << startPointCloud.size()
            << " particles over distance " << maxTrackLen
            << " to find the starting cell" << endl;
    }
    startPointCloud.move(startPointCloud, td, maxTrackLen);


    // Reset levels
    maxFeatureLevel = -1;
    forAll(features_, featI)
    {
        featureEdgeVisited[featI] = false;
    }


    // Start with empty cloud
    Cloud<trackedParticle> cloud(mesh_, Foam::zero{}, "featureCloud");

    if (debug&meshRefinement::FEATURESEEDS)
    {
        Pout<< "Constructing cloud for cell marking" << endl;
    }

    for (const trackedParticle& startTp : startPointCloud)
    {
        const label featI = startTp.i();
        const label pointI = startTp.j();

        const edgeMesh& featureMesh = features_[featI];
        const labelList& pEdges = featureMesh.pointEdges()[pointI];

        // Now shoot particles down all pEdges.
        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];

            if (featureEdgeVisited[featI].set(edgeI))
            {
                // Unvisited edge. Make the particle go to the other point
                // on the edge.

                const edge& e = featureMesh.edges()[edgeI];
                label otherPointi = e.otherVertex(pointI);

                trackedParticle* tp(new trackedParticle(startTp));
                tp->start() = tp->position();
                tp->end() = featureMesh.points()[otherPointi];
                tp->j() = otherPointi;
                tp->k() = edgeI;

                if (debug&meshRefinement::FEATURESEEDS)
                {
                    Pout<< "Adding particle for point:" << pointI
                        << " coord:" << tp->position()
                        << " feature:" << featI
                        << " to track to:" << tp->end()
                        << endl;
                }

                cloud.addParticle(tp);
            }
        }
    }

    startPointCloud.clear();


    while (true)
    {
        // Track all particles to their end position.
        if (debug&meshRefinement::FEATURESEEDS)
        {
            Pout<< "Tracking " << cloud.size()
                << " particles over distance " << maxTrackLen
                << " to mark cells" << endl;
        }
        cloud.move(cloud, td, maxTrackLen);

        // Make particle follow edge.
        for (trackedParticle& tp : cloud)
        {
            const label featI = tp.i();
            const label pointI = tp.j();

            const edgeMesh& featureMesh = features_[featI];
            const labelList& pEdges = featureMesh.pointEdges()[pointI];

            // Particle now at pointI. Check connected edges to see which one
            // we have to visit now.

            bool keepParticle = false;

            forAll(pEdges, i)
            {
                label edgeI = pEdges[i];

                if (featureEdgeVisited[featI].set(edgeI))
                {
                    // Unvisited edge. Make the particle go to the other point
                    // on the edge.

                    const edge& e = featureMesh.edges()[edgeI];
                    label otherPointi = e.otherVertex(pointI);

                    tp.start() = tp.position();
                    tp.end() = featureMesh.points()[otherPointi];
                    tp.j() = otherPointi;
                    tp.k() = edgeI;
                    keepParticle = true;
                    break;
                }
            }

            if (!keepParticle)
            {
                // Particle at 'knot' where another particle already has been
                // seeded. Delete particle.
                cloud.deleteParticle(tp);
            }
        }


        if (debug&meshRefinement::FEATURESEEDS)
        {
            Pout<< "Remaining particles " << cloud.size() << endl;
        }

        if (!returnReduceOr(cloud.size()))
        {
            break;
        }
    }



    //if (debug&meshRefinement::FEATURESEEDS)
    //{
    //    forAll(maxFeatureLevel, cellI)
    //    {
    //        if (maxFeatureLevel[cellI] != -1)
    //        {
    //            Pout<< "Feature went through cell:" << cellI
    //                << " coord:" << mesh_.cellCentres()[cellI]
    //                << " level:" << maxFeatureLevel[cellI]
    //                << endl;
    //        }
    //    }
    //}
}


// Calculates list of cells to refine based on intersection with feature edge.
Foam::label Foam::meshRefinement::markFeatureRefinement
(
    const pointField& keepPoints,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    // Largest refinement level of any feature passed through
    labelList maxFeatureLevel;
    markFeatureCellLevel(keepPoints, maxFeatureLevel);

    // See which cells to refine. maxFeatureLevel will hold highest level
    // of any feature edge that passed through.

    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    forAll(maxFeatureLevel, cellI)
    {
        if (maxFeatureLevel[cellI] > cellLevel[cellI])
        {
            // Mark
            if
            (
               !markForRefine
                (
                    0,                      // surface (n/a)
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                )
            )
            {
                // Reached limit
                break;
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine,  sumOp<label>());
}


// Mark cells for distance-to-feature based refinement.
Foam::label Foam::meshRefinement::markInternalDistanceToFeatureRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Detect if there are any distance shells
    if (features_.maxDistance() <= 0.0)
    {
        return 0;
    }

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList maxLevel;
    features_.findHigherLevel(testCc, testLevels, maxLevel);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Mark cells for non-surface intersection based refinement.
Foam::label Foam::meshRefinement::markInternalRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList maxLevel;
    shells_.findHigherLevel(testCc, testLevels, maxLevel);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::label Foam::meshRefinement::unmarkInternalRefinement
(
    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(nRefine);
    labelList testLevels(nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] >= 0)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList levelShell;
    limitShells_.findLevel(testCc, testLevels, levelShell);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] >= 0)
        {
            if (levelShell[testI] != -1)
            {
                refineCell[cellI] = -1;
                nRefine--;
            }
            testI++;
        }
    }

    return returnReduce(oldNRefine-nRefine, sumOp<label>());
}


// Collect faces that are intersected and whose neighbours aren't yet marked
// for refinement.
Foam::labelList Foam::meshRefinement::getRefineCandidateFaces
(
    const labelList& refineCell
) const
{
    labelList testFaces(mesh_.nFaces());

    label nTest = 0;

    const labelList& surfIndex = surfaceIndex();

    labelList boundaryRefineCell;
    syncTools::swapBoundaryCellList(mesh_, refineCell, boundaryRefineCell);

    forAll(surfIndex, faceI)
    {
        if (surfIndex[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                if (refineCell[own] == -1 || refineCell[nei] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
            else
            {
                const label bFacei = faceI - mesh_.nInternalFaces();
                if (refineCell[own] == -1 || boundaryRefineCell[bFacei] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
        }
    }
    testFaces.setSize(nTest);

    return testFaces;
}


// Mark cells for surface intersection based refinement.
Foam::label Foam::meshRefinement::markSurfaceRefinement
(
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    // Use cached surfaceIndex_ to detect if any intersection. If so
    // re-intersect to determine level wanted.

    // Collect candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~

    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    // ~~~~~~~~~~~~~~~~

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


    // Do test for higher intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surfaceHit;
    labelList surfaceMinLevel;
    surfaces_.findHigherIntersection
    (
        shells_,
        start,
        end,
        minLevel,

        surfaceHit,
        surfaceMinLevel
    );


    // Mark cells for refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label surfI = surfaceHit[i];

        if (surfI != -1)
        {
            // Found intersection with surface with higher wanted
            // refinement. Check if the level field on that surface
            // specifies an even higher level. Note:this is weird. Should
            // do the check with the surfaceMinLevel whilst intersecting the
            // surfaces?

            label own = mesh_.faceOwner()[faceI];

            if (surfaceMinLevel[i] > cellLevel[own])
            {
                // Owner needs refining
                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    break;
                }
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];
                if (surfaceMinLevel[i] > cellLevel[nei])
                {
                    // Neighbour needs refining
                    if
                    (
                       !markForRefine
                        (
                            surfI,
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Count number of matches of first argument in second argument
Foam::label Foam::meshRefinement::countMatches
(
    const List<point>& normals1,
    const List<point>& normals2,
    const scalar tol
) const
{
    label nMatches = 0;

    forAll(normals1, i)
    {
        const vector& n1 = normals1[i];

        forAll(normals2, j)
        {
            const vector& n2 = normals2[j];

            if (magSqr(n1-n2) < tol)
            {
                nMatches++;
                break;
            }
        }
    }

    return nMatches;
}


//bool Foam::meshRefinement::highCurvature
//(
//    const scalar minCosAngle,
//    const point& p0,
//    const vector& n0,
//    const point& p1,
//    const vector& n1
//) const
//{
//    return ((n0&n1) < minCosAngle);
//}
bool Foam::meshRefinement::highCurvature
(
    const scalar minCosAngle,
    const scalar lengthScale,
    const point& p0,
    const vector& n0,
    const point& p1,
    const vector& n1
) const
{
    const scalar cosAngle = (n0&n1);

    if (cosAngle < minCosAngle)
    {
        // Sharp feature, independent of intersection points
        return true;
    }
    else if (cosAngle > 1-1e-6)
    {
        // Co-planar
        return false;
    }
    else if (lengthScale > SMALL)
    {
        // Calculate radius of curvature

        vector axis(n0 ^ n1);
        const plane pl0(p0, (n0 ^ axis));
        const scalar r1 = pl0.normalIntersect(p1, n1);

        //const point radiusPoint(p1+r1*n1);
        //DebugVar(radiusPoint);
        const plane pl1(p1, (n1 ^ axis));
        const scalar r0 = pl1.normalIntersect(p0, n0);

        //const point radiusPoint(p0+r0*n0);
        //DebugVar(radiusPoint);
        //- Take average radius. Bit ad hoc
        const scalar radius = 0.5*(mag(r1)+mag(r0));

        if (radius < lengthScale)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}
//XXXXX


// Mark cells for surface curvature based refinement.
Foam::label Foam::meshRefinement::markSurfaceCurvatureRefinement
(
    const scalar curvature,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // 1. local test: any cell on more than one surface gets refined
    // (if its current level is < max of the surface max level)

    // 2. 'global' test: any cell on only one surface with a neighbour
    // on a different surface gets refined (if its current level etc.)


    const bitSet isMasterFace(syncTools::getMasterFaces(mesh_));


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined.
    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    // Note: uses isMasterFace otherwise could be call to calcCellCellRays
    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            label nei = mesh_.faceNeighbour()[faceI];

            start[i] = cellCentres[own];
            end[i] = cellCentres[nei];
            minLevel[i] = min(cellLevel[own], cellLevel[nei]);
        }
        else
        {
            label bFaceI = faceI - mesh_.nInternalFaces();

            start[i] = cellCentres[own];
            end[i] = neiCc[bFaceI];

            if (!isMasterFace[faceI])
            {
                std::swap(start[i], end[i]);
            }

            minLevel[i] = min(cellLevel[own], neiLevel[bFaceI]);
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // If no curvature supplied behave as before
    const bool hasCurvatureLevels = (max(surfaces_.maxCurvatureLevel()) > 0);


    // Test for all intersections (with surfaces of higher max level than
    // minLevel) and cache per cell the interesting inter
    labelListList cellSurfLevels(mesh_.nCells());
    List<pointList> cellSurfLocations(mesh_.nCells());
    List<vectorList> cellSurfNormals(mesh_.nCells());

    {
        // Per segment the intersection point of the surfaces
        List<pointList> surfaceLocation;
        // Per segment the normals of the surfaces
        List<vectorList> surfaceNormal;
        // Per segment the list of levels of the surfaces
        labelListList surfaceLevel;

        surfaces_.findAllIntersections
        (
            start,
            end,
            minLevel,           // max level of surface has to be bigger
                                // than min level of neighbouring cells

            labelList(surfaces_.maxLevel().size(), 0),  // min level
            surfaces_.maxLevel(),

            surfaceLocation,
            surfaceNormal,
            surfaceLevel
        );


        // Sort the data according to surface location. This will guarantee
        // that on coupled faces both sides visit the intersections in
        // the same order so will decide the same
        labelList visitOrder;
        forAll(surfaceNormal, pointI)
        {
            pointList& pLocations = surfaceLocation[pointI];
            vectorList& pNormals = surfaceNormal[pointI];
            labelList& pLevel = surfaceLevel[pointI];

            sortedOrder(pNormals, visitOrder, normalLess(pLocations));

            pLocations = List<point>(pLocations, visitOrder);
            pNormals = List<vector>(pNormals, visitOrder);
            pLevel = labelUIndList(pLevel, visitOrder);
        }


        //- At some point could just return the intersected surface+region
        //  and derive all the surface information (maxLevel, curvatureLevel)
        //  from that - we're now doing it inside refinementSurfaces itself.
        //// Per segment the surfaces hit
        //List<labelList> hitSurface;
        //List<pointList> hitLocation;
        //List<labelList> hitRegion;
        //List<vectorField> hitNormal;
        //surfaces_.findAllIntersections
        //(
        //    identity(surfaces_.surfaces()), // all refinement geometries
        //    start,
        //    end,
        //
        //    hitSurface,
        //    hitLocation,
        //    hitRegion,
        //    hitNormal
        //);
        //
        //// Filter out levels. minLevel = (mesh) cellLevel (on inbetween face).
        //// Ignore any surface with higher level
        //const auto& maxLevel = surfaces_.maxLevel();
        //labelList visitOrder;
        //DynamicList<label> valid;
        //forAll(hitSurface, segmenti)
        //{
        //    const label meshLevel = minLevel[segmenti];
        //
        //    auto& fSurface = hitSurface[segmenti];
        //    auto& fLocation = hitLocation[segmenti];
        //    auto& fRegion = hitRegion[segmenti];
        //    auto& fNormal = hitNormal[segmenti];
        //
        //    // Sort the data according to intersection location. This will
        //    // guarantee
        //    // that on coupled faces both sides visit the intersections in
        //    // the same order so will decide the same
        //    sortedOrder(fLocation, visitOrder, normalLess(hfLocation));
        //    fLocation = List<point>(fLocation, visitOrder);
        //    fSurface = labelUIndList(fSurface, visitOrder);
        //    fRegion = labelUIndList(fRegion, visitOrder);
        //    fNormal = List<vector>(fNormal, visitOrder);
        //
        //    // Filter out any intersections with surfaces outside cell level.
        //    // Note that min refinement level of surfaces is ignored.
        //    valid.clear();
        //    forAll(fSurface, hiti)
        //    {
        //        const label regioni =
        //            surfaces_.globalRegion(fSurface[hiti], fRegion[hiti]);
        //        if (meshLevel < maxLevel[regioni]) //&& >= minLevel(regioni)
        //        {
        //            valid.append(hiti);
        //        }
        //    }
        //    fLocation = List<point>(fLocation, valid);
        //    fSurface = labelUIndList(fSurface, valid);
        //    fRegion = labelUIndList(fRegion, valid);
        //    fNormal = List<vector>(fNormal, valid);
        //}

        // Clear out unnecessary data
        start.clear();
        end.clear();
        minLevel.clear();

        // Convert face-wise data to cell.
        forAll(testFaces, i)
        {
            label faceI = testFaces[i];
            label own = mesh_.faceOwner()[faceI];

            const pointList& fPoints = surfaceLocation[i];
            const vectorList& fNormals = surfaceNormal[i];
            const labelList& fLevels = surfaceLevel[i];

            forAll(fNormals, hitI)
            {
                if (fLevels[hitI] > cellLevel[own])
                {
                    cellSurfLevels[own].append(fLevels[hitI]);
                    cellSurfLocations[own].append(fPoints[hitI]);
                    cellSurfNormals[own].append(fNormals[hitI]);
                }

                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (fLevels[hitI] > cellLevel[nei])
                    {
                        cellSurfLevels[nei].append(fLevels[hitI]);
                        cellSurfLocations[nei].append(fPoints[hitI]);
                        cellSurfNormals[nei].append(fNormals[hitI]);
                    }
                }
            }
        }
    }



    // Bit of statistics
    if (debug)
    {
        label nSet = 0;
        label nNormals = 0;
        forAll(cellSurfNormals, cellI)
        {
            const vectorList& normals = cellSurfNormals[cellI];
            if (normals.size())
            {
                nSet++;
                nNormals += normals.size();
            }
        }
        reduce(nSet, sumOp<label>());
        reduce(nNormals, sumOp<label>());
        Info<< "markSurfaceCurvatureRefinement :"
            << " cells:" << mesh_.globalData().nTotalCells()
            << " of which with normals:" << nSet
            << " ; total normals stored:" << nNormals
            << endl;
    }



    bool reachedLimit = false;


    // 1. Check normals per cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    for
    (
        label cellI = 0;
        !reachedLimit && cellI < cellSurfLocations.size();
        cellI++
    )
    {
        const pointList& points = cellSurfLocations[cellI];
        const vectorList& normals = cellSurfNormals[cellI];
        const labelList& levels = cellSurfLevels[cellI];

        // Current cell size
        const scalar cellSize =
            meshCutter_.level0EdgeLength()/pow(2.0, cellLevel[cellI]);

        // n^2 comparison of all normals in a cell
        for (label i = 0; !reachedLimit && i < normals.size(); i++)
        {
            for (label j = i+1; !reachedLimit && j < normals.size(); j++)
            {
                // TBD: calculate curvature size (if curvatureLevel specified)
                //      and pass in instead of cellSize

                //if ((normals[i] & normals[j]) < curvature)
                if
                (
                    highCurvature
                    (
                        curvature,
                        (hasCurvatureLevels ? cellSize : scalar(0)),
                        points[i],
                        normals[i],
                        points[j],
                        normals[j]
                    )
                )
                {
                    const label maxLevel = max(levels[i], levels[j]);

                    if (cellLevel[cellI] < maxLevel)
                    {
                        if
                        (
                            !markForRefine
                            (
                                maxLevel,
                                nAllowRefine,
                                refineCell[cellI],
                                nRefine
                            )
                        )
                        {
                            if (debug)
                            {
                                Pout<< "Stopped refining since reaching my cell"
                                    << " limit of " << mesh_.nCells()+7*nRefine
                                    << endl;
                            }
                            reachedLimit = true;
                            break;
                        }
                    }
                }
            }
        }
    }



    // 2. Find out a measure of surface curvature
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Look at normals between neighbouring surfaces
    // Loop over all faces. Could only be checkFaces, except if they're coupled

    // Internal faces
    for
    (
        label faceI = 0;
        !reachedLimit && faceI < mesh_.nInternalFaces();
        faceI++
    )
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        const pointList& ownPoints = cellSurfLocations[own];
        const vectorList& ownNormals = cellSurfNormals[own];
        const labelList& ownLevels = cellSurfLevels[own];
        const pointList& neiPoints = cellSurfLocations[nei];
        const vectorList& neiNormals = cellSurfNormals[nei];
        const labelList& neiLevels = cellSurfLevels[nei];


        // Special case: owner normals set is a subset of the neighbour
        // normals. Do not do curvature refinement since other cell's normals
        // do not add any information. Avoids weird corner extensions of
        // refinement regions.
        bool ownIsSubset =
            countMatches(ownNormals, neiNormals)
         == ownNormals.size();

        bool neiIsSubset =
            countMatches(neiNormals, ownNormals)
         == neiNormals.size();


        if (!ownIsSubset && !neiIsSubset)
        {
            // Current average cell size. Is min? max? average?
            const scalar cellSize =
                meshCutter_.level0EdgeLength()
              / pow(2.0, min(cellLevel[own], cellLevel[nei]));

            // n^2 comparison of between ownNormals and neiNormals
            for (label i = 0; !reachedLimit &&  i < ownNormals.size(); i++)
            {
                for (label j = 0; !reachedLimit && j < neiNormals.size(); j++)
                {
                    // Have valid data on both sides. Check curvature.
                    //if ((ownNormals[i] & neiNormals[j]) < curvature)
                    if
                    (
                        highCurvature
                        (
                            curvature,
                            (hasCurvatureLevels ? cellSize : scalar(0)),
                            ownPoints[i],
                            ownNormals[i],
                            neiPoints[j],
                            neiNormals[j]
                        )
                    )
                    {
                        // See which side to refine.
                        if (cellLevel[own] < ownLevels[i])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    ownLevels[i],
                                    nAllowRefine,
                                    refineCell[own],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                        if (cellLevel[nei] < neiLevels[j])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    neiLevels[j],
                                    nAllowRefine,
                                    refineCell[nei],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }


    // Send over surface point/normal to neighbour cell.
//    labelListList neiSurfaceLevel;
//    syncTools::swapBoundaryCellList(mesh_, cellSurfLevels, neiSurfaceLevel);
    List<vectorList> neiSurfacePoints;
    syncTools::swapBoundaryCellList(mesh_, cellSurfLocations, neiSurfacePoints);
    List<vectorList> neiSurfaceNormals;
    syncTools::swapBoundaryCellList(mesh_, cellSurfNormals, neiSurfaceNormals);

    // Boundary faces
    for
    (
        label faceI = mesh_.nInternalFaces();
        !reachedLimit && faceI < mesh_.nFaces();
        faceI++
    )
    {
        label own = mesh_.faceOwner()[faceI];
        label bFaceI = faceI - mesh_.nInternalFaces();

        const pointList& ownPoints = cellSurfLocations[own];
        const vectorList& ownNormals = cellSurfNormals[own];
        const labelList& ownLevels = cellSurfLevels[own];

        const pointList& neiPoints = neiSurfacePoints[bFaceI];
        const vectorList& neiNormals = neiSurfaceNormals[bFaceI];
        //const labelList& neiLevels = neiSurfaceLevel[bFacei];

        // Special case: owner normals set is a subset of the neighbour
        // normals. Do not do curvature refinement since other cell's normals
        // do not add any information. Avoids weird corner extensions of
        // refinement regions.
        bool ownIsSubset =
            countMatches(ownNormals, neiNormals)
         == ownNormals.size();

        bool neiIsSubset =
            countMatches(neiNormals, ownNormals)
         == neiNormals.size();


        if (!ownIsSubset && !neiIsSubset)
        {
            // Current average cell size. Is min? max? average?
            const scalar cellSize =
                meshCutter_.level0EdgeLength()
              / pow(2.0, min(cellLevel[own], neiLevel[bFaceI]));

            // n^2 comparison of between ownNormals and neiNormals
            for (label i = 0; !reachedLimit &&  i < ownNormals.size(); i++)
            {
                for (label j = 0; !reachedLimit && j < neiNormals.size(); j++)
                {
                    // Have valid data on both sides. Check curvature.
                    //if ((ownNormals[i] & neiNormals[j]) < curvature)
                    if
                    (
                        highCurvature
                        (
                            curvature,
                            (hasCurvatureLevels ? cellSize : scalar(0)),
                            ownPoints[i],
                            ownNormals[i],
                            neiPoints[j],
                            neiNormals[j]
                        )
                    )
                    {
                        if (cellLevel[own] < ownLevels[i])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    ownLevels[i],
                                    nAllowRefine,
                                    refineCell[own],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine
                                        << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }


    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


bool Foam::meshRefinement::isGap
(
    const scalar planarCos,
    const vector& point0,
    const vector& normal0,

    const vector& point1,
    const vector& normal1
) const
{
    //- Hits differ and angles are oppositeish and
    //  hits have a normal distance
    vector d = point1-point0;
    scalar magD = mag(d);

    if (magD > mergeDistance())
    {
        scalar cosAngle = (normal0 & normal1);

        vector avg = Zero;
        if (cosAngle < (-1+planarCos))
        {
            // Opposite normals
            avg = 0.5*(normal0-normal1);
        }
        else if (cosAngle > (1-planarCos))
        {
            avg = 0.5*(normal0+normal1);
        }

        if (avg != vector::zero)
        {
            avg /= mag(avg);

            // Check normal distance of intersection locations
            if (mag(avg&d) > mergeDistance())
            {
                return true;
            }
        }
    }

    return false;
}


// Mark small gaps
bool Foam::meshRefinement::isNormalGap
(
    const scalar planarCos,
    const label level0,
    const vector& point0,
    const vector& normal0,

    const label level1,
    const vector& point1,
    const vector& normal1
) const
{
    //const scalar edge0Len = meshCutter_.level0EdgeLength();

    //- Hits differ and angles are oppositeish and
    //  hits have a normal distance
    vector d = point1-point0;
    scalar magD = mag(d);

    if (magD > mergeDistance())
    {
        scalar cosAngle = (normal0 & normal1);

        vector avgNormal = Zero;
        if (cosAngle < (-1+planarCos))
        {
            // Opposite normals
            avgNormal = 0.5*(normal0-normal1);
        }
        else if (cosAngle > (1-planarCos))
        {
            avgNormal = 0.5*(normal0+normal1);
        }

        if (avgNormal != vector::zero)
        {
            avgNormal /= mag(avgNormal);

            //scalar normalDist = mag(d&avgNormal);
            //const scalar avgCellSize = edge0Len/pow(2.0, 0.5*(level0+level1));
            //if (normalDist > 1*avgCellSize)
            //{
            //    Pout<< "** skipping large distance " << endl;
            //    return false;
            //}

            // Check average normal with respect to intersection locations
            if (mag(avgNormal&d/magD) > Foam::cos(degToRad(45.0)))
            {
                return true;
            }
        }
    }

    return false;
}


bool Foam::meshRefinement::checkProximity
(
    const scalar planarCos,
    const label nAllowRefine,

    const label surfaceLevel,       // current intersection max level
    const vector& surfaceLocation,  // current intersection location
    const vector& surfaceNormal,    // current intersection normal

    const label cellI,

    label& cellMaxLevel,        // cached max surface level for this cell
    vector& cellMaxLocation,    // cached surface location for this cell
    vector& cellMaxNormal,      // cached surface normal for this cell

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    // Test if surface applicable
    if (surfaceLevel > cellLevel[cellI])
    {
        if (cellMaxLevel == -1)
        {
            // First visit of cell. Store
            cellMaxLevel = surfaceLevel;
            cellMaxLocation = surfaceLocation;
            cellMaxNormal = surfaceNormal;
        }
        else
        {
            // Second or more visit.
            // Check if
            //  - different location
            //  - opposite surface

            bool closeSurfaces = isNormalGap
            (
                planarCos,
                cellLevel[cellI],   //cellMaxLevel,
                cellMaxLocation,
                cellMaxNormal,
                cellLevel[cellI],   //surfaceLevel,
                surfaceLocation,
                surfaceNormal
            );

            // Set normal to that of highest surface. Not really necessary
            // over here but we reuse cellMax info when doing coupled faces.
            if (surfaceLevel > cellMaxLevel)
            {
                cellMaxLevel = surfaceLevel;
                cellMaxLocation = surfaceLocation;
                cellMaxNormal = surfaceNormal;
            }


            if (closeSurfaces)
            {
                //Pout<< "Found gap:" << nl
                //    << "    location:" << surfaceLocation
                //    << "\tnormal  :" << surfaceNormal << nl
                ///    << "    location:" << cellMaxLocation
                //    << "\tnormal  :" << cellMaxNormal << nl
                //    << "\tcos     :" << (surfaceNormal&cellMaxNormal) << nl
                //    << endl;

                return markForRefine
                (
                    surfaceLevel,   // mark with any non-neg number.
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );
            }
        }
    }

    // Did not reach refinement limit.
    return true;
}


Foam::label Foam::meshRefinement::markProximityRefinement
(
    const scalar planarCos,

    const labelList& surfaceMinLevel,
    const labelList& surfaceMaxLevel,

    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    // 1. local test: any cell on more than one surface gets refined
    // (if its current level is < max of the surface max level)

    // 2. 'global' test: any cell on only one surface with a neighbour
    // on a different surface gets refined (if its current level etc.)


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined.
    const labelList testFaces(getRefineCandidateFaces(refineCell));

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


    // Test for all intersections (with surfaces of higher gap level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.
    labelList cellMaxLevel(mesh_.nCells(), -1);
    vectorField cellMaxNormal(mesh_.nCells(), Zero);
    pointField cellMaxLocation(mesh_.nCells(), Zero);

    {
        // Per segment the normals of the surfaces
        List<vectorList> surfaceLocation;
        List<vectorList> surfaceNormal;
        // Per segment the list of levels of the surfaces
        labelListList surfaceLevel;

        surfaces_.findAllIntersections
        (
            start,
            end,
            minLevel,           // gap level of surface has to be bigger
                                // than min level of neighbouring cells

            surfaceMinLevel,
            surfaceMaxLevel,

            surfaceLocation,
            surfaceNormal,
            surfaceLevel
        );
        // Clear out unnecessary data
        start.clear();
        end.clear();
        minLevel.clear();


        //// Extract per cell information on the surface with the highest max

        forAll(testFaces, i)
        {
            label faceI = testFaces[i];
            label own = mesh_.faceOwner()[faceI];

            const labelList& fLevels = surfaceLevel[i];
            const vectorList& fPoints = surfaceLocation[i];
            const vectorList& fNormals = surfaceNormal[i];

            forAll(fLevels, hitI)
            {
                checkProximity
                (
                    planarCos,
                    nAllowRefine,

                    fLevels[hitI],
                    fPoints[hitI],
                    fNormals[hitI],

                    own,
                    cellMaxLevel[own],
                    cellMaxLocation[own],
                    cellMaxNormal[own],

                    refineCell,
                    nRefine
                );
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                forAll(fLevels, hitI)
                {
                    checkProximity
                    (
                        planarCos,
                        nAllowRefine,

                        fLevels[hitI],
                        fPoints[hitI],
                        fNormals[hitI],

                        nei,
                        cellMaxLevel[nei],
                        cellMaxLocation[nei],
                        cellMaxNormal[nei],

                        refineCell,
                        nRefine
                    );
                }
            }
        }
    }

    // 2. Find out a measure of surface curvature and region edges.
    // Send over surface region and surface normal to neighbour cell.

    labelList neiBndMaxLevel(mesh_.nBoundaryFaces());
    pointField neiBndMaxLocation(mesh_.nBoundaryFaces());
    vectorField neiBndMaxNormal(mesh_.nBoundaryFaces());

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label bFaceI = faceI-mesh_.nInternalFaces();
        label own = mesh_.faceOwner()[faceI];

        neiBndMaxLevel[bFaceI] = cellMaxLevel[own];
        neiBndMaxLocation[bFaceI] = cellMaxLocation[own];
        neiBndMaxNormal[bFaceI] = cellMaxNormal[own];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxLevel);
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxLocation);
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxNormal);

    // Loop over all faces. Could only be checkFaces.. except if they're coupled

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        if (cellMaxLevel[own] != -1 && cellMaxLevel[nei] != -1)
        {
            // Have valid data on both sides. Check planarCos.
            if
            (
                isNormalGap
                (
                    planarCos,
                    cellLevel[own], //cellMaxLevel[own],
                    cellMaxLocation[own],
                    cellMaxNormal[own],
                    cellLevel[nei], //cellMaxLevel[nei],
                    cellMaxLocation[nei],
                    cellMaxNormal[nei]
                )
            )
            {
                // See which side to refine
                if (cellLevel[own] < cellMaxLevel[own])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[own],
                            nAllowRefine,
                            refineCell[own],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }

                if (cellLevel[nei] < cellMaxLevel[nei])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[nei],
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }
            }
        }
    }
    // Boundary faces
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label bFaceI = faceI - mesh_.nInternalFaces();

        if (cellLevel[own] < cellMaxLevel[own] && neiBndMaxLevel[bFaceI] != -1)
        {
            // Have valid data on both sides. Check planarCos.
            if
            (
                isNormalGap
                (
                    planarCos,
                    cellLevel[own], //cellMaxLevel[own],
                    cellMaxLocation[own],
                    cellMaxNormal[own],
                    neiLevel[bFaceI], //neiBndMaxLevel[bFaceI],
                    neiBndMaxLocation[bFaceI],
                    neiBndMaxNormal[bFaceI]
                )
            )
            {
                if
                (
                    !markForRefine
                    (
                        cellMaxLevel[own],
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining since reaching my cell"
                            << " limit of " << mesh_.nCells()+7*nRefine
                            << endl;
                    }
                    break;
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate list of cells to refine. Gets for any edge (start - end)
// whether it intersects the surface. Does more accurate test and checks
// the wanted level on the surface intersected.
// Does approximate precalculation of how many cells can be refined before
// hitting overall limit maxGlobalCells.
Foam::labelList Foam::meshRefinement::refineCandidates
(
    const pointField& keepPoints,
    const scalar curvature,
    const scalar planarAngle,

    const bool featureRefinement,
    const bool featureDistanceRefinement,
    const bool internalRefinement,
    const bool surfaceRefinement,
    const bool curvatureRefinement,
    const bool smallFeatureRefinement,
    const bool gapRefinement,
    const bool bigGapRefinement,
    const bool spreadGapSize,
    const label maxGlobalCells,
    const label maxLocalCells
) const
{
    label totNCells = mesh_.globalData().nTotalCells();

    labelList cellsToRefine;

    if (totNCells >= maxGlobalCells)
    {
        Info<< "No cells marked for refinement since reached limit "
            << maxGlobalCells << '.' << endl;
    }
    else
    {
        // Every cell I refine adds 7 cells. Estimate the number of cells
        // I am allowed to refine.
        // Assume perfect distribution so can only refine as much the fraction
        // of the mesh I hold. This prediction step prevents us having to do
        // lots of reduces to keep count of the total number of cells selected
        // for refinement.

        //scalar fraction = scalar(mesh_.nCells())/totNCells;
        //label myMaxCells = label(maxGlobalCells*fraction);
        //label nAllowRefine = (myMaxCells - mesh_.nCells())/7;
        ////label nAllowRefine = (maxLocalCells - mesh_.nCells())/7;
        //
        //Pout<< "refineCandidates:" << nl
        //    << "    total cells:" << totNCells << nl
        //    << "    local cells:" << mesh_.nCells() << nl
        //    << "    local fraction:" << fraction << nl
        //    << "    local allowable cells:" << myMaxCells << nl
        //    << "    local allowable refinement:" << nAllowRefine << nl
        //    << endl;

        //- Disable refinement shortcut. nAllowRefine is per processor limit.
        label nAllowRefine = labelMax / Pstream::nProcs();

        // Marked for refinement (>= 0) or not (-1). Actual value is the
        // index of the surface it intersects / shell it is inside.
        labelList refineCell(mesh_.nCells(), -1);
        label nRefine = 0;


        // Swap neighbouring cell centres and cell level
        labelList neiLevel(mesh_.nBoundaryFaces());
        pointField neiCc(mesh_.nBoundaryFaces());
        calcNeighbourData(neiLevel, neiCc);


        const scalar planarCos = Foam::cos(degToRad(planarAngle));


        // Cells pierced by feature lines
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureRefinement)
        {
            label nFeatures = markFeatureRefinement
            (
                keepPoints,
                nAllowRefine,

                refineCell,
                nRefine
            );

            Info<< "Marked for refinement due to explicit features             "
                << ": " << nFeatures << " cells."  << endl;
        }

        // Inside distance-to-feature shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureDistanceRefinement)
        {
            label nShell = markInternalDistanceToFeatureRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to distance to explicit features "
                ": " << nShell << " cells."  << endl;
        }

        // Inside refinement shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        if (internalRefinement)
        {
            label nShell = markInternalRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to refinement shells             "
                << ": " << nShell << " cells."  << endl;
        }

        // Refinement based on intersection of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (surfaceRefinement)
        {
            label nSurf = markSurfaceRefinement
            (
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to surface intersection          "
                << ": " << nSurf << " cells."  << endl;

            // Refine intersected-cells only inside gaps. See
            // markInternalGapRefinement to refine all cells inside gaps.
            if
            (
                planarCos >= -1
             && planarCos <= 1
             && max(shells_.maxGapLevel()) > 0
            )
            {
                label nGapSurf = markSurfaceGapRefinement
                (
                    planarCos,
                    nAllowRefine,
                    neiLevel,
                    neiCc,

                    refineCell,
                    nRefine
                );
                Info<< "Marked for refinement due to surface intersection"
                    << " (at gaps)"
                    << ": " << nGapSurf << " cells."  << endl;
            }
        }

        // Refinement based on curvature of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            curvatureRefinement
         && (curvature >= -1 && curvature <= 1)
         && (surfaces_.minLevel() != surfaces_.maxLevel())
        )
        {
            label nCurv = markSurfaceCurvatureRefinement
            (
                curvature,
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to curvature/regions             "
                << ": " << nCurv << " cells."  << endl;
        }


        // Refinement based on features smaller than cell size
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (smallFeatureRefinement)
        {
            label nGap = 0;
            if
            (
                planarCos >= -1
             && planarCos <= 1
             && max(shells_.maxGapLevel()) > 0
            )
            {
                nGap = markSmallFeatureRefinement
                (
                    planarCos,
                    nAllowRefine,
                    neiLevel,
                    neiCc,

                    refineCell,
                    nRefine
                );
            }
            Info<< "Marked for refinement due to close opposite surfaces       "
                << ": " << nGap << " cells."  << endl;

            label nCurv = 0;
            if (max(surfaces_.maxCurvatureLevel()) > 0)
            {
                nCurv = markSurfaceFieldRefinement
                (
                    nAllowRefine,
                    neiLevel,
                    neiCc,

                    refineCell,
                    nRefine
                );
                Info<< "Marked for refinement"
                    << " due to curvature                     "
                    << ": " << nCurv << " cells."  << endl;
            }
        }


        // Refinement based on gap (only neighbouring cells)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const labelList& surfaceGapLevel = surfaces_.gapLevel();

        if
        (
            gapRefinement
         && (planarCos >= -1 && planarCos <= 1)
         && (max(surfaceGapLevel) > -1)
        )
        {
            Info<< "Specified gap level : " << max(surfaceGapLevel)
                << ", planar angle " << planarAngle << endl;

            label nGap = markProximityRefinement
            (
                planarCos,

                labelList(surfaceGapLevel.size(), -1),  // surface min level
                surfaceGapLevel,                        // surface max level

                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to close opposite surfaces       "
                << ": " << nGap << " cells."  << endl;
        }


        // Refinement based on gaps larger than cell size
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            bigGapRefinement
         && (planarCos >= -1 && planarCos <= 1)
         && max(shells_.maxGapLevel()) > 0
        )
        {
            // Refine based on gap information provided by shell and nearest
            // surface
            labelList numGapCells;
            scalarField gapSize;
            label nGap = markInternalGapRefinement
            (
                planarCos,
                spreadGapSize,
                nAllowRefine,

                refineCell,
                nRefine,
                numGapCells,
                gapSize
            );
            Info<< "Marked for refinement due to opposite surfaces"
                << "             "
                << ": " << nGap << " cells."  << endl;
        }


        // Limit refinement
        // ~~~~~~~~~~~~~~~~

        if (limitShells_.shells().size())
        {
            label nUnmarked = unmarkInternalRefinement(refineCell, nRefine);
            if (nUnmarked > 0)
            {
                Info<< "Unmarked for refinement due to limit shells"
                    << "                : " << nUnmarked << " cells."  << endl;
            }
        }



        // Pack cells-to-refine
        // ~~~~~~~~~~~~~~~~~~~~

        cellsToRefine.setSize(nRefine);
        nRefine = 0;

        forAll(refineCell, cellI)
        {
            if (refineCell[cellI] != -1)
            {
                cellsToRefine[nRefine++] = cellI;
            }
        }
    }

    return cellsToRefine;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Remove any unnecessary fields
    mesh_.clearOut();
    mesh_.moving(false);

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false);
    mapPolyMesh& map = *mapPtr;

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map.hasMotionPoints())
    {
        mesh_.movePoints(map.preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, cellsToRefine));

    return mapPtr;
}


// Load balancing
Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    labelList& cellsToRefine,
    const scalar maxLoadUnbalance,
    const label maxCellUnbalance
)
{
    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        // First check if we need to balance at all. Precalculate number of
        // cells after refinement and see what maximum difference is.
        const scalar nNewCells =
            scalar(mesh_.nCells() + 7*cellsToRefine.size());
        const scalar nNewCellsAll = returnReduce(nNewCells, sumOp<scalar>());
        const scalar nIdealNewCells = nNewCellsAll / Pstream::nProcs();
        const scalar unbalance = returnReduce
        (
            mag(1.0-nNewCells/nIdealNewCells),
            maxOp<scalar>()
        );

        // Trigger the balancing to avoid too early balancing for better
        // scaling performance.
        const scalar nNewCellsOnly = scalar(7*cellsToRefine.size());

        const label maxNewCells =
            label(returnReduce(nNewCellsOnly, maxOp<scalar>()));

        const label maxDeltaCells =
            label(mag(returnReduce(nNewCells, maxOp<scalar>())-nIdealNewCells));

        // New trigger to avoid too early balancing
        // 1. Check if globally one proc exceeds the maxCellUnbalance value
        //    related to the new added cells at the refinement loop
        // 2. Check if globally one proc exceeds the maxCellUnbalance based on
        //    the average cell count a proc should have
        if
        (
            (maxNewCells <= maxCellUnbalance)
         && (maxDeltaCells <= maxCellUnbalance)
        )
        {
            Info<< "Skipping balancing since trigger value not reached:" << "\n"
                << "    Trigger cell count: " << maxCellUnbalance << nl
                << "    Max new cell count in proc: " << maxNewCells << nl
                << "    Max difference between new cells and balanced: "
                << maxDeltaCells << nl
                << "    Max load unbalance " << maxLoadUnbalance
                << nl <<endl;
        }
        else
        {
            if (unbalance <= maxLoadUnbalance)
            {
                Info<< "Skipping balancing since max unbalance " << unbalance
                    << " is less than allowable " << maxLoadUnbalance
                    << endl;
            }
            else
            {
                Info<< "Balancing since max unbalance " << unbalance
                    << " is larger than allowable " << maxLoadUnbalance
                    << endl;

                scalarField cellWeights(mesh_.nCells(), 1);
                forAll(cellsToRefine, i)
                {
                    cellWeights[cellsToRefine[i]] += 7;
                }

                distMap = balance
                (
                    false,  //keepZoneFaces
                    false,  //keepBaffles
                    labelList::null(),  //singleProcPoints
                    cellWeights,
                    decomposer,
                    distributor
                );

                // Update cells to refine
                distMap().distributeCellIndices(cellsToRefine);

                Info<< "Balanced mesh in = "
                    << mesh_.time().cpuTimeIncrement() << " s" << endl;

                printMeshInfo(debug, "After balancing " + msg, true);


                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Writing balanced " << msg
                        << " mesh to time " << timeName() << endl;
                    write
                    (
                        debugType(debug),
                        writeType(writeLevel() | WRITEMESH),
                        mesh_.time().path()/timeName()
                    );
                    Pout<< "Dumped debug data in = "
                        << mesh_.time().cpuTimeIncrement() << " s" << endl;

                    // test all is still synced across proc patches
                    checkData();
                }
            }
        }
    }

    return distMap;
}


// Do refinement of consistent set of cells followed by truncation and
// load balancing.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::refineAndBalance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& cellsToRefine,
    const scalar maxLoadUnbalance,
    const label maxCellUnbalance
)
{
    // Refinement
    // ~~~~~~~~~~

    refine(cellsToRefine);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined but unbalanced " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    printMeshInfo(debug, "After refinement " + msg, true);


    // Load balancing
    // ~~~~~~~~~~~~~~

    labelList noCellsToRefine;

    auto distMap = balance
    (
        msg,
        decomposer,
        distributor,
        noCellsToRefine,    // mesh is already refined; no need to predict
        maxLoadUnbalance,
        maxCellUnbalance
    );

    return distMap;
}


// Do load balancing followed by refinement of consistent set of cells.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::balanceAndRefine
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& initCellsToRefine,
    const scalar maxLoadUnbalance,
    const label maxCellUnbalance
)
{
    labelList cellsToRefine(initCellsToRefine);

    //{
    //    globalIndex globalCells(mesh_.nCells());
    //
    //    Info<< "** Distribution before balancing/refining:" << endl;
    //    for (const int procI : Pstream::allProcs())
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}
    //{
    //    globalIndex globalCells(cellsToRefine.size());
    //
    //    Info<< "** Cells to be refined:" << endl;
    //    for (const int procI : Pstream::allProcs())
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}



    // Load balancing
    // ~~~~~~~~~~~~~~

    auto distMap = balance
    (
        msg,
        decomposer,
        distributor,
        cellsToRefine,
        maxLoadUnbalance,
        maxCellUnbalance
    );


    // Refinement
    // ~~~~~~~~~~
    // Note: uses updated cellsToRefine

    refine(cellsToRefine);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;

    //{
    //    globalIndex globalCells(mesh_.nCells());
    //
    //    Info<< "** After refinement distribution:" << endl;
    //    for (const int procI : Pstream::allProcs())
    //    {
    //        Info<< "    " << procI << '\t'
    //            << globalCells.localSize(procI) << endl;
    //    }
    //    Info<< endl;
    //}

    printMeshInfo(debug, "After refinement " + msg, true);

    return distMap;
}


Foam::labelList Foam::meshRefinement::directionalRefineCandidates
(
    const label maxGlobalCells,
    const label maxLocalCells,
    const labelList& currentLevel,
    const direction dir
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label totNCells = mesh_.globalData().nTotalCells();

    labelList cellsToRefine;

    if (totNCells >= maxGlobalCells)
    {
        Info<< "No cells marked for refinement since reached limit "
            << maxGlobalCells << '.' << endl;
    }
    else
    {
        // Disable refinement shortcut. nAllowRefine is per processor limit.
        label nAllowRefine = labelMax / Pstream::nProcs();

        // Marked for refinement (>= 0) or not (-1). Actual value is the
        // index of the surface it intersects / shell it is inside
        labelList refineCell(mesh_.nCells(), -1);
        label nRefine = 0;

        // Find cells inside the shells with directional levels
        labelList insideShell;
        shells_.findDirectionalLevel
        (
            cellCentres,
            cellLevel,
            currentLevel,  // current directional level
            dir,
            insideShell
        );

        // Mark for refinement
        forAll(insideShell, celli)
        {
            if (insideShell[celli] >= 0)
            {
                bool reachedLimit = !markForRefine
                (
                    insideShell[celli],    // mark with any positive value
                    nAllowRefine,
                    refineCell[celli],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+nAllowRefine << endl;
                    }
                    break;
                }
            }
        }

        // Limit refinement
        // ~~~~~~~~~~~~~~~~

        {
            label nUnmarked = unmarkInternalRefinement(refineCell, nRefine);
            if (nUnmarked > 0)
            {
                Info<< "Unmarked for refinement due to limit shells"
                    << "                : " << nUnmarked << " cells."  << endl;
            }
        }



        // Pack cells-to-refine
        // ~~~~~~~~~~~~~~~~~~~~

        cellsToRefine.setSize(nRefine);
        nRefine = 0;

        forAll(refineCell, cellI)
        {
            if (refineCell[cellI] != -1)
            {
                cellsToRefine[nRefine++] = cellI;
            }
        }
    }

    return cellsToRefine;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::directionalRefine
(
    const string& msg,
    const direction cmpt,
    const labelList& cellsToRefine
)
{
    // Set splitting direction
    vector refDir(Zero);
    refDir[cmpt] = 1;
    List<refineCell> refCells(cellsToRefine.size());
    forAll(cellsToRefine, i)
    {
        refCells[i] = refineCell(cellsToRefine[i], refDir);
    }

    // How to walk circumference of cells
    hexCellLooper cellWalker(mesh_);

    // Analyse cuts
    cellCuts cuts(mesh_, cellWalker, refCells);

    // Cell cutter
    Foam::meshCutter meshRefiner(mesh_);

    polyTopoChange meshMod(mesh_);

    // Insert mesh refinement into polyTopoChange.
    meshRefiner.setRefinement(cuts, meshMod);

    // Remove any unnecessary fields
    mesh_.clearOut();
    mesh_.moving(false);

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(*morphMap);

    // Move mesh (since morphing does not do this)
    if (morphMap().hasMotionPoints())
    {
        mesh_.movePoints(morphMap().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update stored refinement pattern
    meshRefiner.updateMesh(*morphMap);

    // Update intersection info
    updateMesh(*morphMap, getChangedFaces(*morphMap, cellsToRefine));

    return morphMap;
}


// ************************************************************************* //
