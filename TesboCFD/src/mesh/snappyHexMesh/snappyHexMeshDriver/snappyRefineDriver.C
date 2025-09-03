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

\*---------------------------------------------------------------------------*/

#include "snappyRefineDriver.H"
#include "meshRefinement.H"
#include "fvMesh.H"
#include "Time.H"
#include "cellSet.H"
#include "syncTools.H"
#include "refinementParameters.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "mapDistributePolyMesh.H"
#include "unitConversion.H"
#include "snapParameters.H"
#include "localPointRegion.H"
#include "IOmanip.H"
#include "labelVector.H"
#include "profiling.H"
#include "searchableSurfaces.H"
#include "fvMeshSubset.H"
#include "interpolationTable.H"
#include "snappyVoxelMeshDriver.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "addPatchCellLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(snappyRefineDriver, 0);
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyRefineDriver::snappyRefineDriver
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelUList& globalToMasterPatch,
    const labelUList& globalToSlavePatch,
    coordSetWriter& setFormatter,
    refPtr<surfaceWriter>& surfFormatter,
    const bool dryRun
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch),
    setFormatter_(setFormatter),
    surfFormatter_(surfFormatter),
    dryRun_(dryRun)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::snappyRefineDriver::featureEdgeRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine
)
{
    if (dryRun_)
    {
        return 0;
    }

    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(edge, "snappyHexMesh::refine::edge");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (meshRefiner_.features().size() && maxIter > 0)
    {
        for (; iter < maxIter; iter++)
        {
            Info<< nl
                << "Feature refinement iteration " << iter << nl
                << "------------------------------" << nl
                << endl;

            labelList candidateCells
            (
                meshRefiner_.refineCandidates
                (
                    refineParams.locationsInMesh(),
                    refineParams.curvature(),
                    refineParams.planarAngle(),

                    true,               // featureRefinement
                    false,              // featureDistanceRefinement
                    false,              // internalRefinement
                    false,              // surfaceRefinement
                    false,              // curvatureRefinement
                    false,              // smallFeatureRefinement
                    false,              // gapRefinement
                    false,              // bigGapRefinement
                    false,              // spreadGapSize
                    refineParams.maxGlobalCells(),
                    refineParams.maxLocalCells()
                )
            );
            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;



            const label nCellsToRefine =
                returnReduce(cellsToRefine.size(), sumOp<label>());

            Info<< "Selected for feature refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            if (nCellsToRefine <= minRefine)
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug > 0)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
            {
                meshRefiner_.balanceAndRefine
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::smallFeatureRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (dryRun_)
    {
        return 0;
    }

    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(feature, "snappyHexMesh::refine::smallFeature");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    // See if any surface has an extendedGapLevel
    const labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    const labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());
    const labelList curvMaxLevel(meshRefiner_.surfaces().maxCurvatureLevel());

    if
    (
        max(surfaceMaxLevel) == 0
     && max(shellMaxLevel) == 0
     && max(curvMaxLevel) == 0
    )
    {
        return iter;
    }

    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Small surface feature refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                true,               // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if (nCellsToRefine == 0)
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::surfaceOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label leakBlockageIter
)
{
    if (dryRun_)
    {
        return 0;
    }

    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(surface, "snappyHexMesh::refine::surface");
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minimum number of surface refinement iterations.
    label overallMaxLevel = max(meshRefiner_.surfaces().maxLevel());

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Surface refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Do optional leak closing (by removing cells)
        if (iter >= leakBlockageIter)
        {
            // Block off intersections with unzoned surfaces with specified
            // leakLevel < iter
            const labelList unnamedSurfaces
            (
                surfaceZonesInfo::getUnnamedSurfaces
                (
                    surfaces.surfZones()
                )
            );

            DynamicList<label> selectedSurfaces(unnamedSurfaces.size());
            for (const label surfi : unnamedSurfaces)
            {
                const label regioni = surfaces.globalRegion(surfi, 0);

                // Take shortcut: assume all cells on surface are refined to
                // its refinement level at iteration iter. So just use the
                // iteration to see if the surface is a candidate.
                if (iter > surfaces.leakLevel()[regioni])
                {
                    selectedSurfaces.append(surfi);
                }
            }

            if
            (
                selectedSurfaces.size()
             && refineParams.locationsOutsideMesh().size()
            )
            {
               meshRefiner_.blockLeakFaces
                (
                    globalToMasterPatch_,
                    globalToSlavePatch_,
                    refineParams.locationsInMesh(),
                    refineParams.zonesInMesh(),
                    refineParams.locationsOutsideMesh(),
                    selectedSurfaces
                );
            }
        }


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                true,               // surfaceRefinement
                true,               // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );
        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::gapOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (dryRun_)
    {
        return 0;
    }

    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minimum number of surface refinement iterations.

    label maxIncrement = 0;
    const labelList& maxLevel = meshRefiner_.surfaces().maxLevel();
    const labelList& gapLevel = meshRefiner_.surfaces().gapLevel();

    forAll(maxLevel, i)
    {
        maxIncrement = max(maxIncrement, gapLevel[i]-maxLevel[i]);
    }

    label iter = 0;

    if (maxIncrement == 0)
    {
        return iter;
    }

    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Gap refinement iteration " << iter << nl
            << "--------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                true,               // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing current mesh to time "
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
            Pout<< "Dumped mesh in = "
                << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGap." << endl;
            cellSet c(mesh, "candidateCellsFromGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Grow by one layer to make sure we're covering the gap
        {
            boolList isCandidateCell(mesh.nCells(), false);
            forAll(candidateCells, i)
            {
                isCandidateCell[candidateCells[i]] = true;
            }

            for (label i=0; i<1; i++)
            {
                boolList newIsCandidateCell(isCandidateCell);

                // Internal faces
                for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
                {
                    label own = mesh.faceOwner()[facei];
                    label nei = mesh.faceNeighbour()[facei];

                    if (isCandidateCell[own] != isCandidateCell[nei])
                    {
                        newIsCandidateCell[own] = true;
                        newIsCandidateCell[nei] = true;
                    }
                }

                // Get coupled boundary condition values
                boolList neiIsCandidateCell;
                syncTools::swapBoundaryCellList
                (
                    mesh,
                    isCandidateCell,
                    neiIsCandidateCell
                );

                // Boundary faces
                for
                (
                    label facei = mesh.nInternalFaces();
                    facei < mesh.nFaces();
                    facei++
                )
                {
                    label own = mesh.faceOwner()[facei];
                    label bFacei = facei-mesh.nInternalFaces();

                    if (isCandidateCell[own] != neiIsCandidateCell[bFacei])
                    {
                        newIsCandidateCell[own] = true;
                    }
                }

                isCandidateCell.transfer(newIsCandidateCell);
            }

            label n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    n++;
                }
            }
            candidateCells.setSize(n);
            n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    candidateCells[n++] = celli;
                }
            }
        }


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGapPlusBuffer." << endl;
            cellSet c(mesh, "candidateCellsFromGapPlusBuffer", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }


        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= maxIncrement
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::surfaceProximityBlock
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    fvMesh& mesh = meshRefiner_.mesh();

    if (min(meshRefiner_.surfaces().blockLevel()) == labelMax)
    {
        return 0;
    }

    label iter = 0;

    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Gap blocking iteration " << iter << nl
            << "------------------------" << nl
            << endl;


        // Determine cells to block
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        meshRefiner_.removeGapCells
        (
            refineParams.planarAngle(),
            meshRefiner_.surfaces().blockLevel(),
            globalToMasterPatch_,
            refineParams.nFilterIter()
        );

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
            Pout<< "Writing gap blocking iteration "
                << iter << " mesh to time " << meshRefiner_.timeName()
                << endl;
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
    return iter;
}


Foam::label Foam::snappyRefineDriver::bigGapOnlyRefine
(
    const refinementParameters& refineParams,
    const bool spreadGapSize,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    if (dryRun_)
    {
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    // See if any surface has an extendedGapLevel
    labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());

    label overallMaxLevel(max(max(surfaceMaxLevel), max(shellMaxLevel)));

    if (overallMaxLevel == 0)
    {
        return iter;
    }


    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Big gap refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                true,               // bigGapRefinement
                spreadGapSize,      // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing current mesh to time "
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
            Pout<< "Dumped mesh in = "
                << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromBigGap." << endl;
            cellSet c(mesh, "candidateCellsFromBigGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::danglingCellRefine
(
    const refinementParameters& refineParams,
    const label nFaces,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    if (dryRun_)
    {
        return 0;
    }

    addProfiling(dangling, "snappyHexMesh::refine::danglingCell");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Dangling coarse cells refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        const cellList& cells = mesh.cells();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        labelList candidateCells;
        {
            cellSet candidateCellSet(mesh, "candidateCells", cells.size()/1000);

            forAll(cells, celli)
            {
                label nIntFaces = 0;
                for (const label meshFacei : cells[celli])
                {
                    const label patchi = pbm.patchID(meshFacei);

                    if (patchi < 0 || pbm[patchi].coupled())
                    {
                        ++nIntFaces;
                    }
                }

                if (nIntFaces == nFaces)
                {
                    candidateCellSet.insert(celli);
                }
            }

            if (debug&meshRefinement::MESH)
            {
                Pout<< "Dumping " << candidateCellSet.size()
                    << " cells to cellSet candidateCellSet." << endl;
                candidateCellSet.instance() = meshRefiner_.timeName();
                candidateCellSet.write();
            }
            candidateCells = candidateCellSet.toc();
        }



        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine. After a few iterations check if too
        // few cells
        if
        (
            nCellsToRefine == 0
         || (
                iter >= 1
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;

            if (refineParams.balanceAtEnd())
            {
                Info<< "Final mesh balancing" << endl;

                meshRefiner_.balance
                (
                    "",
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    0,
                    -1
                );
            }

            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    return iter;
}


// Detect cells with opposing intersected faces of differing refinement
// level and refine them.
Foam::label Foam::snappyRefineDriver::refinementInterfaceRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    if (dryRun_)
    {
        return 0;
    }

    addProfiling(interface, "snappyHexMesh::refine::transition");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (refineParams.interfaceRefine())
    {
        for (;iter < maxIter; iter++)
        {
            Info<< nl
                << "Refinement transition refinement iteration " << iter << nl
                << "--------------------------------------------" << nl
                << endl;

            const labelList& surfaceIndex = meshRefiner_.surfaceIndex();
            const hexRef8& cutter = meshRefiner_.meshCutter();
            const vectorField& fA = mesh.faceAreas();
            const labelList& faceOwner = mesh.faceOwner();


            // Determine cells to refine
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            const cellList& cells = mesh.cells();

            labelList candidateCells;
            {
                // Pass1: pick up cells with differing face level

                cellSet transitionCells
                (
                    mesh,
                    "transitionCells",
                    cells.size()/100
                );

                forAll(cells, celli)
                {
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if (surfaceIndex[facei] != -1)
                        {
                            label fLevel = cutter.faceLevel(facei);
                            if (fLevel != cLevel)
                            {
                                transitionCells.insert(celli);
                            }
                        }
                    }
                }


                cellSet candidateCellSet
                (
                    mesh,
                    "candidateCells",
                    cells.size()/1000
                );

                // Pass2: check for oppositeness

                //for (const label celli : transitionCells)
                //{
                //    const cell& cFaces = cells[celli];
                //    const point& cc = cellCentres[celli];
                //    const scalar rCVol = pow(cellVolumes[celli], -5.0/3.0);
                //
                //    // Determine principal axes of cell
                //    symmTensor R(Zero);
                //
                //    forAll(cFaces, i)
                //    {
                //        label facei = cFaces[i];
                //
                //        const point& fc = faceCentres[facei];
                //
                //        // Calculate face-pyramid volume
                //        scalar pyrVol = 1.0/3.0 * fA[facei] & (fc-cc);
                //
                //        if (faceOwner[facei] != celli)
                //        {
                //            pyrVol = -pyrVol;
                //        }
                //
                //        // Calculate face-pyramid centre
                //        vector pc = (3.0/4.0)*fc + (1.0/4.0)*cc;
                //
                //        R += pyrVol*sqr(pc-cc)*rCVol;
                //    }
                //
                //    //- MEJ: Problem: truncation errors cause complex evs
                //    vector lambdas(eigenValues(R));
                //    const tensor axes(eigenVectors(R, lambdas));
                //
                //
                //    // Check if this cell has
                //    // - opposing sides intersected
                //    // - which are of different refinement level
                //    // - plus the inbetween face
                //
                //    labelVector plusFaceLevel(labelVector(-1, -1, -1));
                //    labelVector minFaceLevel(labelVector(-1, -1, -1));
                //
                //    forAll(cFaces, cFacei)
                //    {
                //        label facei = cFaces[cFacei];
                //
                //        if (surfaceIndex[facei] != -1)
                //        {
                //            label fLevel = cutter.faceLevel(facei);
                //
                //            // Get outwards pointing normal
                //            vector n = fA[facei]/mag(fA[facei]);
                //            if (faceOwner[facei] != celli)
                //            {
                //                n = -n;
                //            }
                //
                //            // What is major direction and sign
                //            direction cmpt = vector::X;
                //            scalar maxComp = (n&axes.x());
                //
                //            scalar yComp = (n&axes.y());
                //            scalar zComp = (n&axes.z());
                //
                //            if (mag(yComp) > mag(maxComp))
                //            {
                //                maxComp = yComp;
                //                cmpt = vector::Y;
                //            }
                //
                //            if (mag(zComp) > mag(maxComp))
                //            {
                //                maxComp = zComp;
                //                cmpt = vector::Z;
                //            }
                //
                //            if (maxComp > 0)
                //            {
                //                plusFaceLevel[cmpt] = max
                //                (
                //                    plusFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //            else
                //            {
                //                minFaceLevel[cmpt] = max
                //                (
                //                    minFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //        }
                //    }
                //
                //    // Check if we picked up any opposite differing level
                //    for (direction dir = 0; dir < vector::nComponents; dir++)
                //    {
                //        if
                //        (
                //            plusFaceLevel[dir] != -1
                //         && minFaceLevel[dir] != -1
                //         && plusFaceLevel[dir] != minFaceLevel[dir]
                //        )
                //        {
                //            candidateCellSet.insert(celli);
                //        }
                //    }
                //}

                const scalar oppositeCos = Foam::cos(degToRad(135.0));

                for (const label celli : transitionCells)
                {
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    // Detect opposite intersection
                    bool foundOpposite = false;

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if
                        (
                            surfaceIndex[facei] != -1
                         && cutter.faceLevel(facei) > cLevel
                        )
                        {
                            // Get outwards pointing normal
                            vector n = fA[facei]/mag(fA[facei]);
                            if (faceOwner[facei] != celli)
                            {
                                n = -n;
                            }

                            // Check for any opposite intersection
                            forAll(cFaces, cFaceI2)
                            {
                                label face2i = cFaces[cFaceI2];

                                if
                                (
                                    face2i != facei
                                 && surfaceIndex[face2i] != -1
                                 && cutter.faceLevel(face2i) > cLevel
                                )
                                {
                                    // Get outwards pointing normal
                                    vector n2 = fA[face2i]/mag(fA[face2i]);
                                    if (faceOwner[face2i] != celli)
                                    {
                                        n2 = -n2;
                                    }


                                    if ((n&n2) < oppositeCos)
                                    {
                                        foundOpposite = true;
                                        break;
                                    }
                                }
                            }

                            if (foundOpposite)
                            {
                                break;
                            }
                        }
                    }


                    if (foundOpposite)
                    {
                        candidateCellSet.insert(celli);
                    }
                }

                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Dumping " << candidateCellSet.size()
                        << " cells to cellSet candidateCellSet." << endl;
                    candidateCellSet.instance() = meshRefiner_.timeName();
                    candidateCellSet.write();
                }
                candidateCells = candidateCellSet.toc();
            }



            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;


            const label nCellsToRefine =
                returnReduce(cellsToRefine.size(), sumOp<label>());

            Info<< "Selected for refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            // Stop when no cells to refine. After a few iterations check if too
            // few cells
            if
            (
                nCellsToRefine == 0
             || (
                    iter >= 1
                 && nCellsToRefine <= refineParams.minRefineCells()
                )
            )
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
            {
                meshRefiner_.balanceAndRefine
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
        }
    }
    return iter;
}

bool Foam::snappyRefineDriver::usesHigherLevel
(
    const labelUList& boundaryPointLevel,
    const labelUList& f,
    const label cLevel
) const
{
    for (const label pointi : f)
    {
        if (boundaryPointLevel[pointi] > cLevel)
        {
            return true;
        }
    }
    return false;
}


Foam::label Foam::snappyRefineDriver::boundaryRefinementInterfaceRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    if (dryRun_)
    {
        return 0;
    }

    addProfiling(interface, "snappyHexMesh::refine::transition");
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (refineParams.interfaceRefine())
    {
        for (;iter < maxIter; iter++)
        {
            Info<< nl
                << "Boundary refinement iteration " << iter << nl
                << "-------------------------------" << nl
                << endl;

            const labelList& surfaceIndex = meshRefiner_.surfaceIndex();
            const hexRef8& cutter = meshRefiner_.meshCutter();
            const labelList& cellLevel = cutter.cellLevel();
            const faceList& faces = mesh.faces();
            const cellList& cells = mesh.cells();


            // Determine cells to refine
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            // Point/face on boundary
            bitSet isBoundaryPoint(mesh.nPoints());
            bitSet isBoundaryFace(mesh.nFaces());
            {
                forAll(surfaceIndex, facei)
                {
                    if (surfaceIndex[facei] != -1)
                    {
                        isBoundaryFace.set(facei);
                        isBoundaryPoint.set(faces[facei]);
                    }
                }
                const labelList meshPatchIDs(meshRefiner_.meshedPatches());
                for (const label patchi : meshPatchIDs)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchi];
                    forAll(pp, i)
                    {
                        isBoundaryFace.set(pp.start()+i);
                        isBoundaryPoint.set(pp[i]);
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    isBoundaryPoint,
                    orEqOp<unsigned int>(),
                    0
                );
            }

            // Mark max boundary face level onto boundary points. All points
            // not on a boundary face stay 0.
            labelList boundaryPointLevel(mesh.nPoints(), 0);
            {
                forAll(cells, celli)
                {
                    const cell& cFaces = cells[celli];
                    const label cLevel = cellLevel[celli];

                    for (const label facei : cFaces)
                    {
                        if (isBoundaryFace(facei))
                        {
                            const face& f = faces[facei];
                            for (const label pointi : f)
                            {
                                boundaryPointLevel[pointi] =
                                max
                                (
                                    boundaryPointLevel[pointi],
                                    cLevel
                                );
                            }
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    boundaryPointLevel,
                    maxEqOp<label>(),
                    label(0)
                );
            }


            // Detect cells with a point but not face on the boundary
            labelList candidateCells;
            {
                const cellList& cells = mesh.cells();

                cellSet candidateCellSet
                (
                    mesh,
                    "candidateCells",
                    cells.size()/100
                );

                forAll(cells, celli)
                {
                    const cell& cFaces = cells[celli];
                    const label cLevel = cellLevel[celli];

                    bool isBoundaryCell = false;
                    for (const label facei : cFaces)
                    {
                        if (isBoundaryFace(facei))
                        {
                            isBoundaryCell = true;
                            break;
                        }
                    }

                    if (!isBoundaryCell)
                    {
                        for (const label facei : cFaces)
                        {
                            const face& f = mesh.faces()[facei];
                            if (usesHigherLevel(boundaryPointLevel, f, cLevel))
                            {
                                candidateCellSet.insert(celli);
                                break;
                            }
                        }
                    }
                }
                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Dumping " << candidateCellSet.size()
                        << " cells to cellSet candidateCellSet." << endl;
                    candidateCellSet.instance() = meshRefiner_.timeName();
                    candidateCellSet.write();
                }
                candidateCells = candidateCellSet.toc();
            }

            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;


            const label nCellsToRefine =
                returnReduce(cellsToRefine.size(), sumOp<label>());

            Info<< "Selected for refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            // Stop when no cells to refine. After a few iterations check if too
            // few cells
            if
            (
                nCellsToRefine == 0
                // || (
                //        iter >= 1
                //     && nCellsToRefine <= refineParams.minRefineCells()
                //    )
            )
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
            {
                meshRefiner_.balanceAndRefine
                (
                    "boundary cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "boundary cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams.maxLoadUnbalance(),
                    refineParams.maxCellUnbalance()
                );
            }
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::removeInsideCells
(
    const refinementParameters& refineParams,
    const label nBufferLayers
)
{
    // Skip if no limitRegion and zero bufferLayers
    if (meshRefiner_.limitShells().shells().size() == 0 && nBufferLayers == 0)
    {
        return;
    }

    if (dryRun_)
    {
        return;
    }

    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    // Remove any cells inside limitShells with level -1
    if (meshRefiner_.limitShells().shells().size())
    {
        meshRefiner_.removeLimitShells
        (
            nBufferLayers,
            1,
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.zonesInMesh(),
            refineParams.locationsOutsideMesh()
        );
    }


    // Fix any additional (e.g. locationsOutsideMesh). Note: probably not
    // necessary.
    meshRefiner_.splitMesh
    (
        nBufferLayers,                  // nBufferLayers
        refineParams.nErodeCellZone(),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh(),
        !refineParams.useLeakClosure(),
        setFormatter_,
        surfFormatter_
    );

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;

        Pout<< "Writing subsetted mesh to time "
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
        Pout<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::label Foam::snappyRefineDriver::shellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (dryRun_)
    {
        return 0;
    }

    if (refineParams.minRefineCells() == -1)
    {
        // Special setting to be able to restart shm on meshes with inconsistent
        // cellLevel/pointLevel
        return 0;
    }

    addProfiling(shell, "snappyHexMesh::refine::shell");
    const fvMesh& mesh = meshRefiner_.mesh();

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = ListOps::createWithValue<label>
    (
        mesh.nFaces(),
        meshRefiner_.intersectedFaces(),
        0, // set value
        -1 // default value
    );

    // Determine the maximum refinement level over all volume refinement
    // regions. This determines the minimum number of shell refinement
    // iterations.
    label overallMaxShellLevel = meshRefiner_.shells().maxLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Shell refinement iteration " << iter << nl
            << "----------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                refineParams.locationsInMesh(),
                refineParams.curvature(),
                refineParams.planarAngle(),

                false,              // featureRefinement
                true,               // featureDistanceRefinement
                true,               // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false,              // spreadGapSize
                refineParams.maxGlobalCells(),
                refineParams.maxLocalCells()
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromShells." << endl;

            cellSet c(mesh, "candidateCellsFromShells", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Problem choosing starting faces for bufferlayers (bFaces)
        //  - we can't use the current intersected boundary faces
        //    (intersectedFaces) since this grows indefinitely
        //  - if we use 0 faces we don't satisfy bufferLayers from the
        //    surface.
        //  - possibly we want to have bFaces only the initial set of faces
        //    and maintain the list while doing the refinement.
        labelList bFaces
        (
            findIndices(meshRefiner_.userFaceData()[0].second(), 0)
        );

        //Info<< "Collected boundary faces : "
        //    << returnReduce(bFaces.size(), sumOp<label>()) << endl;

        labelList cellsToRefine;

        if (refineParams.nBufferLayers() <= 2)
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement
            (
                refineParams.nBufferLayers(),
                candidateCells,                     // cells to refine
                bFaces,                             // faces for nBufferLayers
                1,                                  // point difference
                meshRefiner_.intersectedPoints()    // points to check
            );
        }
        else
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement2
            (
                refineParams.nBufferLayers(),
                candidateCells,                 // cells to refine
                bFaces                          // faces for nBufferLayers
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        const label nCellsToRefine =
            returnReduce(cellsToRefine.size(), sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxShellLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if (returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells()))
        {
            meshRefiner_.balanceAndRefine
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams.maxLoadUnbalance(),
                refineParams.maxCellUnbalance()
            );
        }
    }
    meshRefiner_.userFaceData().clear();

    return iter;
}


Foam::label Foam::snappyRefineDriver::directionalShellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (dryRun_)
    {
        return 0;
    }

    addProfiling(shell, "snappyHexMesh::refine::directionalShell");
    const fvMesh& mesh = meshRefiner_.mesh();
    const shellSurfaces& shells = meshRefiner_.shells();

    labelList& cellLevel =
        const_cast<labelIOList&>(meshRefiner_.meshCutter().cellLevel());
    labelList& pointLevel =
        const_cast<labelIOList&>(meshRefiner_.meshCutter().pointLevel());


    // Determine the minimum and maximum cell levels that are candidates for
    // directional refinement
    const labelPairList dirSelect(shells.directionalSelectLevel());
    label overallMinLevel = labelMax;
    label overallMaxLevel = labelMin;
    forAll(dirSelect, shelli)
    {
        overallMinLevel = min(dirSelect[shelli].first(), overallMinLevel);
        overallMaxLevel = max(dirSelect[shelli].second(), overallMaxLevel);
    }

    if (overallMinLevel > overallMaxLevel)
    {
        return 0;
    }

    // Maintain directional refinement levels
    List<labelVector> dirCellLevel(cellLevel.size());
    forAll(cellLevel, celli)
    {
        dirCellLevel[celli] = labelVector::uniform(cellLevel[celli]);
    }

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Directional shell refinement iteration " << iter << nl
            << "----------------------------------------" << nl
            << endl;

        label nAllRefine = 0;

        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            // Select the cells that need to be refined in certain direction:
            // - cell inside/outside shell
            // - original cellLevel (using mapping) mentioned in levelIncrement
            // - dirCellLevel not yet up to cellLevel+levelIncrement


            // Extract component of directional level
            labelList currentLevel(dirCellLevel.size());
            forAll(dirCellLevel, celli)
            {
                currentLevel[celli] = dirCellLevel[celli][dir];
            }

            labelList candidateCells
            (
                meshRefiner_.directionalRefineCandidates
                (
                    refineParams.maxGlobalCells(),
                    refineParams.maxLocalCells(),
                    currentLevel,
                    dir
                )
            );

            // Extend to keep 2:1 ratio
            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    currentLevel,
                    candidateCells,
                    true
                )
            );

            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;

            const label nCellsToRefine =
                returnReduce(cellsToRefine.size(), sumOp<label>());

            Info<< "Selected for direction " << vector::componentNames[dir]
                << " refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            nAllRefine += nCellsToRefine;

            // Stop when no cells to refine or have done minimum necessary
            // iterations and not enough cells to refine.
            if (nCellsToRefine > 0)
            {
                if (debug)
                {
                    const_cast<Time&>(mesh.time())++;
                }

                const bitSet isRefineCell(mesh.nCells(), cellsToRefine);

                autoPtr<mapPolyMesh> map
                (
                    meshRefiner_.directionalRefine
                    (
                        "directional refinement iteration " + name(iter),
                        dir,
                        cellsToRefine
                    )
                );

                Info<< "Refined mesh in = "
                    << mesh.time().cpuTimeIncrement() << " s" << endl;

                meshRefinement::updateList
                (
                    map().cellMap(),
                    labelVector(0, 0, 0),
                    dirCellLevel
                );

                // Note: edges will have been split. The points might have
                // inherited pointLevel from either side of the edge which
                // might not be the same for coupled edges so sync
                syncTools::syncPointList
                (
                    mesh,
                    pointLevel,
                    maxEqOp<label>(),
                    labelMin
                );

                forAll(map().cellMap(), celli)
                {
                    if (isRefineCell[map().cellMap()[celli]])
                    {
                        dirCellLevel[celli][dir]++;
                    }
                }

                // Do something with the pointLevel. See discussion about the
                // cellLevel. Do we keep min/max ?
                forAll(map().pointMap(), pointi)
                {
                    label oldPointi = map().pointMap()[pointi];
                    if (map().reversePointMap()[oldPointi] != pointi)
                    {
                        // Is added point (splitting an edge)
                        pointLevel[pointi]++;
                    }
                }
            }
        }


        if (nAllRefine == 0)
        {
            Info<< "Stopping refining since no cells selected."
                << nl << endl;
            break;
        }

        meshRefiner_.printMeshInfo
        (
            debug,
            "After directional refinement iteration " + name(iter),
            true
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing directional refinement iteration "
                << iter << " mesh to time " << meshRefiner_.timeName() << endl;
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

    // Adjust cellLevel from dirLevel? As max? Or the min?
    // For now: use max. The idea is that if there is a wall
    // any directional refinement is likely to be aligned with
    // the wall (wall layers) so any snapping/layering would probably
    // want to use this highest refinement level.

    forAll(cellLevel, celli)
    {
        cellLevel[celli] = cmptMax(dirCellLevel[celli]);
    }

    return iter;
}


void Foam::snappyRefineDriver::mergeAndSmoothRatio
(
    const scalarList& allSeedPointDist,
    const label nSmoothExpansion,
    List<Tuple2<scalar, scalar>>&  keyAndValue
)
{
    // Merge duplicate distance from coupled locations to get unique
    // distances to operate on, do on master
    SortableList<scalar> unmergedDist(allSeedPointDist);
    DynamicList<scalar> mergedDist;

    scalar prevDist = GREAT;
    forAll(unmergedDist, i)
    {
        scalar curDist = unmergedDist[i];
        scalar difference = mag(curDist - prevDist);
        if (difference > meshRefiner_.mergeDistance())
        //if (difference > 0.01)
        {
             mergedDist.append(curDist);
             prevDist = curDist;
        }
    }

    // Sort the unique distances
    SortableList<scalar> sortedDist(mergedDist);
    labelList indexSet = sortedDist.indices();

    // Get updated position starting from original (undistorted) mesh
    scalarList seedPointsNewLocation = sortedDist;

    scalar initResidual = 0.0;
    scalar prevIterResidual = GREAT;

    for (label iter = 0; iter < nSmoothExpansion; iter++)
    {

        // Position based edge averaging algorithm operated on
        // all seed plane locations in normalized form.
        //
        //   0   1   2   3   4   5   6  (edge numbers)
        //  ---x---x---x---x---x---x---
        //     0   1   2   3   4   5    (point numbers)
        //
        // Average of edge 1-3 in terms of position
        //  = (point3 - point0)/3
        // Keeping points 0-1 frozen, new position of point 2
        //  = position2 + (average of edge 1-3 as above)
        for(label i = 2; i<mergedDist.size()-1; i++)
        {
            scalar oldX00 = sortedDist[i-2];
            scalar oldX1 = sortedDist[i+1];
            scalar curX0 = seedPointsNewLocation[i-1];
            seedPointsNewLocation[i] = curX0 + (oldX1 - oldX00)/3;
        }

        const scalarField residual(seedPointsNewLocation-sortedDist);
        {
            scalar res(sumMag(residual));

            if (iter == 0)
            {
                initResidual = res;
            }
            res /= initResidual;

            if (mag(prevIterResidual - res) < SMALL)
            {
                if (debug)
                {
                    Pout<< "Converged with iteration " << iter
                        << " initResidual: " << initResidual
                        << " final residual : " << res << endl;
                }
                break;
            }
            else
            {
                prevIterResidual = res;
            }
        }

        // Update the field for next iteration, avoid moving points
        sortedDist = seedPointsNewLocation;

    }

    keyAndValue.setSize(mergedDist.size());

    forAll(mergedDist, i)
    {
        keyAndValue[i].first() = mergedDist[i];
        label index = indexSet[i];
        keyAndValue[i].second() = seedPointsNewLocation[index];
    }
}


Foam::label Foam::snappyRefineDriver::directionalSmooth
(
    const refinementParameters& refineParams
)
{
    addProfiling(split, "snappyHexMesh::refine::smooth");
    Info<< nl
        << "Directional expansion ratio smoothing" << nl
        << "-------------------------------------" << nl
        << endl;

    fvMesh& baseMesh = meshRefiner_.mesh();
    const searchableSurfaces& geometry = meshRefiner_.surfaces().geometry();
    const shellSurfaces& shells = meshRefiner_.shells();

    label iter = 0;

    forAll(shells.nSmoothExpansion(), shellI)
    {
        if
        (
            shells.nSmoothExpansion()[shellI] > 0
         || shells.nSmoothPosition()[shellI] > 0
        )
        {
            label surfi = shells.shells()[shellI];
            const vector& userDirection = shells.smoothDirection()[shellI];


            // Extract inside points
            labelList pointLabels;
            {
                // Get inside points
                List<volumeType> volType;
                geometry[surfi].getVolumeType(baseMesh.points(), volType);

                label nInside = 0;
                forAll(volType, pointi)
                {
                    if (volType[pointi] == volumeType::INSIDE)
                    {
                        nInside++;
                    }
                }
                pointLabels.setSize(nInside);
                nInside = 0;
                forAll(volType, pointi)
                {
                    if (volType[pointi] == volumeType::INSIDE)
                    {
                        pointLabels[nInside++] = pointi;
                    }
                }

                //bitSet isInsidePoint(baseMesh.nPoints());
                //forAll(volType, pointi)
                //{
                //    if (volType[pointi] == volumeType::INSIDE)
                //    {
                //        isInsidePoint.set(pointi);
                //    }
                //}
                //pointLabels = isInsidePoint.used();
            }

            // Mark all directed edges
            bitSet isXEdge(baseMesh.edges().size());
            {
                const edgeList& edges = baseMesh.edges();
                forAll(edges, edgei)
                {
                    const edge& e = edges[edgei];
                    vector eVec(e.vec(baseMesh.points()));
                    eVec /= mag(eVec);
                    if (mag(eVec&userDirection) > 0.9)
                    {
                        isXEdge.set(edgei);
                    }
                }
            }

            // Get the extreme of smoothing region and
            // normalize all points within
            const scalar totalLength =
                geometry[surfi].bounds().span()
              & userDirection;
            const scalar startPosition =
                geometry[surfi].bounds().min()
              & userDirection;

            scalarField normalizedPosition(pointLabels.size(), Zero);
            forAll(pointLabels, i)
            {
                label pointi = pointLabels[i];
                normalizedPosition[i] =
                  (
                    ((baseMesh.points()[pointi]&userDirection) - startPosition)
                  / totalLength
                  );
            }

            // Sort the normalized position
            labelList order(sortedOrder(normalizedPosition));

            DynamicList<scalar> seedPointDist;

            // Select points from finest refinement (one point-per plane)
            scalar prevDist = GREAT;
            forAll(order, i)
            {
                label pointi = order[i];
                scalar curDist = normalizedPosition[pointi];
                if (mag(curDist - prevDist) > meshRefiner_.mergeDistance())
                {
                    seedPointDist.append(curDist);
                    prevDist = curDist;
                }
            }

            // Collect data from all processors
            scalarList allSeedPointDist;
            {
                List<scalarList> gatheredDist(Pstream::nProcs());
                gatheredDist[Pstream::myProcNo()] = seedPointDist;
                Pstream::gatherList(gatheredDist);

                // Combine processor lists into one big list.
                allSeedPointDist =
                    ListListOps::combine<scalarList>
                    (
                        gatheredDist, accessOp<scalarList>()
                    );
            }

            // Pre-set the points not to smooth (after expansion)
            bitSet isFrozenPoint(baseMesh.nPoints(), true);

            {
                scalar minSeed = min(allSeedPointDist);
                scalar maxSeed = max(allSeedPointDist);
                Pstream::broadcasts(UPstream::worldComm, minSeed, maxSeed);

                forAll(normalizedPosition, posI)
                {
                    const scalar pos = normalizedPosition[posI];
                    if
                    (
                        (mag(pos-minSeed) < meshRefiner_.mergeDistance())
                     || (mag(pos-maxSeed) < meshRefiner_.mergeDistance())
                    )
                    {
                        // Boundary point: freeze
                        isFrozenPoint.set(pointLabels[posI]);
                    }
                    else
                    {
                        // Internal to moving region
                        isFrozenPoint.unset(pointLabels[posI]);
                    }
                }
            }

            Info<< "Smoothing " << geometry[surfi].name() << ':' << nl
                << "    Direction                   : " << userDirection << nl
                << "    Number of points            : "
                << returnReduce(pointLabels.size(), sumOp<label>())
                << " (out of " << baseMesh.globalData().nTotalPoints()
                << ")" << nl
                << "    Smooth expansion iterations : "
                << shells.nSmoothExpansion()[shellI] << nl
                << "    Smooth position iterations  : "
                << shells.nSmoothPosition()[shellI] << nl
                << "    Number of planes            : "
                << allSeedPointDist.size()
                << endl;

            // Make lookup from original normalized distance to new value
            List<Tuple2<scalar, scalar>> keyAndValue(allSeedPointDist.size());

            // Filter unique seed distances and iterate for user given steps
            // or until convergence. Then get back map from old to new distance
            if (Pstream::master())
            {
                mergeAndSmoothRatio
                (
                    allSeedPointDist,
                    shells.nSmoothExpansion()[shellI],
                    keyAndValue
                );
            }

            Pstream::broadcast(keyAndValue);

            // Construct an iterpolation table for further queries
            // - although normalized values are used for query,
            //   it might flow out of bounds due to precision, hence clamped
            const interpolationTable<scalar> table
            (
                keyAndValue,
                bounds::repeatableBounding::CLAMP,
                "undefined"
            );

            // Now move the points directly on the baseMesh.
            pointField baseNewPoints(baseMesh.points());
            forAll(pointLabels, i)
            {
                label pointi = pointLabels[i];
                const point& curPoint = baseMesh.points()[pointi];
                scalar curDist = normalizedPosition[i];
                scalar newDist = table(curDist);
                scalar newPosition = startPosition + newDist*totalLength;
                baseNewPoints[pointi] +=
                    userDirection * (newPosition - (curPoint &userDirection));
            }

            // Moving base mesh with expansion ratio smoothing
            vectorField disp(baseNewPoints-baseMesh.points());
            syncTools::syncPointList
            (
                baseMesh,
                disp,
                maxMagSqrEqOp<vector>(),
                point::zero
            );
            baseMesh.movePoints(baseMesh.points()+disp);

            // Reset any moving state
            baseMesh.moving(false);

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(baseMesh.time())++;

                Pout<< "Writing directional expansion ratio smoothed"
                    << " mesh to time " << meshRefiner_.timeName() << endl;

                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    baseMesh.time().path()/meshRefiner_.timeName()
                );
            }

            // Now we have moved the points in user specified region. Smooth
            // them with neighbour points to avoid skewed cells
            // Instead of moving actual mesh, operate on copy
            pointField baseMeshPoints(baseMesh.points());
            scalar initResidual = 0.0;
            scalar prevIterResidual = GREAT;
            for (iter = 0; iter < shells.nSmoothPosition()[shellI]; iter++)
            {
                {
                    const edgeList& edges = baseMesh.edges();
                    const labelListList& pointEdges = baseMesh.pointEdges();

                    pointField unsmoothedPoints(baseMeshPoints);

                    scalarField sumOther(baseMesh.nPoints(), Zero);
                    labelList nSumOther(baseMesh.nPoints(), Zero);
                    labelList nSumXEdges(baseMesh.nPoints(), Zero);
                    forAll(edges, edgei)
                    {
                        const edge& e = edges[edgei];
                        sumOther[e[0]] +=
                            (unsmoothedPoints[e[1]]&userDirection);
                        nSumOther[e[0]]++;
                        sumOther[e[1]] +=
                            (unsmoothedPoints[e[0]]&userDirection);
                        nSumOther[e[1]]++;
                        if (isXEdge[edgei])
                        {
                            nSumXEdges[e[0]]++;
                            nSumXEdges[e[1]]++;
                        }
                    }

                    syncTools::syncPointList
                    (
                        baseMesh,
                        nSumXEdges,
                        plusEqOp<label>(),
                        label(0)
                    );

                    forAll(pointLabels, i)
                    {
                        label pointi = pointLabels[i];

                        if (nSumXEdges[pointi] < 2)
                        {
                            // Hanging node. Remove the (single!) X edge so it
                            // will follow points above or below instead
                            const labelList& pEdges = pointEdges[pointi];
                            forAll(pEdges, pE)
                            {
                                label edgei = pEdges[pE];
                                if (isXEdge[edgei])
                                {
                                    const edge& e = edges[edgei];
                                    label otherPt = e.otherVertex(pointi);
                                    nSumOther[pointi]--;
                                    sumOther[pointi] -=
                                      (
                                          unsmoothedPoints[otherPt]
                                        & userDirection
                                      );
                                }
                            }
                        }
                    }

                    syncTools::syncPointList
                    (
                        baseMesh,
                        sumOther,
                        plusEqOp<scalar>(),
                        scalar(0)
                    );
                    syncTools::syncPointList
                    (
                        baseMesh,
                        nSumOther,
                        plusEqOp<label>(),
                        label(0)
                    );

                    forAll(pointLabels, i)
                    {
                        label pointi = pointLabels[i];

                        if ((nSumOther[pointi] >= 2) && !isFrozenPoint[pointi])
                        {
                            scalar smoothPos =
                                0.5
                               *(
                                    (unsmoothedPoints[pointi]&userDirection)
                                   +sumOther[pointi]/nSumOther[pointi]
                                );

                            vector& v = baseNewPoints[pointi];
                            v += (smoothPos-(v&userDirection))*userDirection;
                        }
                    }

                    const vectorField residual(baseNewPoints - baseMeshPoints);
                    {
                        scalar res(gSum(mag(residual)));

                        if (iter == 0)
                        {
                            initResidual = res;
                        }
                        res /= initResidual;

                        if (mag(prevIterResidual - res) < SMALL)
                        {
                            Info<< "Converged smoothing in iteration " << iter
                                << " initResidual: " << initResidual
                                << " final residual : " << res << endl;
                            break;
                        }
                        else
                        {
                            prevIterResidual = res;
                        }
                    }

                    // Just copy new location instead of moving base mesh
                    baseMeshPoints = baseNewPoints;
                }
            }

            // Move mesh to new location
            vectorField dispSmooth(baseMeshPoints-baseMesh.points());
            syncTools::syncPointList
            (
                baseMesh,
                dispSmooth,
                maxMagSqrEqOp<vector>(),
                point::zero
            );
            baseMesh.movePoints(baseMesh.points()+dispSmooth);

            // Reset any moving state
            baseMesh.moving(false);

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(baseMesh.time())++;

                Pout<< "Writing positional smoothing iteration "
                    << iter << " mesh to time " << meshRefiner_.timeName()
                    << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    baseMesh.time().path()/meshRefiner_.timeName()
                );
            }
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::baffleAndSplitMesh
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    if (dryRun_)
    {
        return;
    }

    addProfiling(split, "snappyHexMesh::refine::splitting");
    Info<< nl
        << "Splitting mesh at surface intersections" << nl
        << "---------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    // Introduce baffles at surface intersections. Note:
    // meshRefinement::surfaceIndex() will
    // be like boundary face from now on so not coupled anymore.
    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,             // detect&remove potential snap problem

        // Snap problem cell detection
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        false,                          // perpendicular edge connected cells
        scalarField(0),                 // per region perpendicular angle
        refineParams.nErodeCellZone(),

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh(),
        !refineParams.useLeakClosure(),
        setFormatter_,
        surfFormatter_
    );


    if (!handleSnapProblems) // merge free standing baffles?
    {
        // By default only merge baffles if on same patch. However the patching
        // of the castellated mesh is not very accurate. For backwards
        // compatibility disable check for new mode only.
        const auto mt = meshRefiner_.meshType();
        const bool mergeSameOnly
        (
            (
                mt == meshRefinement::CASTELLATEDBUFFERLAYER
             || mt == meshRefinement::CASTELLATEDBUFFERLAYER2
            )
          ? false
          : true
        );

        meshRefiner_.mergeFreeStandingBaffles
        (
            mergeSameOnly,      // whether to merge baffles on different patches
            snapParams,
            refineParams.useTopologicalSnapDetection(),
            false,                  // perpendicular edge connected cells
            scalarField(0),         // per region perpendicular angle
            refineParams.planarAngle(),
            motionDict,
            const_cast<Time&>(mesh.time()),
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh()
        );
    }
}


void Foam::snappyRefineDriver::zonify
(
    const refinementParameters& refineParams,
    wordPairHashTable& zonesToFaceZone
)
{
    // Mesh is at its finest. Do zoning
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This puts all faces with intersection across a zoneable surface
    // into that surface's faceZone. All cells inside faceZone get given the
    // same cellZone.

    const labelList namedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(meshRefiner_.surfaces().surfZones());

    if
    (
        namedSurfaces.size()
     || refineParams.zonesInMesh().size()
    )
    {
        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;

        const fvMesh& mesh = meshRefiner_.mesh();

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        meshRefiner_.zonify
        (
            refineParams.allowFreeStandingZoneFaces(),
            refineParams.nErodeCellZone(),
            refineParams.locationsInMesh(),
            refineParams.zonesInMesh(),
            refineParams.locationsOutsideMesh(),
            !refineParams.useLeakClosure(),
            setFormatter_,
            surfFormatter_,
            zonesToFaceZone
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing zoned mesh to time "
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

        // Check that all faces are synced
        meshRefinement::checkCoupledFaceZones(mesh);
    }
}


void Foam::snappyRefineDriver::splitAndMergeBaffles
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    if (dryRun_)
    {
        return;
    }

    Info<< nl
        << "Handling cells with snap problems" << nl
        << "---------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Introduce baffles and split mesh
    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    const scalarField& perpAngle = meshRefiner_.surfaces().perpendicularAngle();

    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,

        // Snap problem cell detection
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        handleSnapProblems,                 // remove perp edge connected cells
        perpAngle,                          // perp angle
        refineParams.nErodeCellZone(),

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.zonesInMesh(),
        refineParams.locationsOutsideMesh(),
        !refineParams.useLeakClosure(),
        setFormatter_,
        surfFormatter_
    );


    // By default only merge baffles if on same patch. However the patching
    // of the castellated mesh is not very accurate. For backwards
    // compatibility disable check for new mode only.
    const auto mt = meshRefiner_.meshType();
    const bool mergeSameOnly
    (
        (
            mt == meshRefinement::CASTELLATEDBUFFERLAYER
         || mt == meshRefinement::CASTELLATEDBUFFERLAYER2
        )
      ? false
      : true
    );

    meshRefiner_.mergeFreeStandingBaffles
    (
        mergeSameOnly,
        snapParams,
        refineParams.useTopologicalSnapDetection(),
        handleSnapProblems,
        perpAngle,
        refineParams.planarAngle(),
        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams.locationsInMesh(),
        refineParams.locationsOutsideMesh()
    );

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();


    // Merge all baffles that are still remaining after duplicating points.
    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh));

    const label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected unsplittable baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        meshRefiner_.mergeBaffles(couples, Map<label>(0));

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        // Remove any now dangling parts
        meshRefiner_.splitMeshRegions
        (
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh(),
            true,   // exit if any path to outside points
            setFormatter_
        );

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        Info<< "Merged free-standing baffles in = "
            << mesh.time().cpuTimeIncrement() << " s." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing handleProblemCells mesh to time "
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


void Foam::snappyRefineDriver::erodeNonManifoldZoneFaces
(
    const refinementParameters& refineParams
)
{
    if (dryRun_)
    {
        return;
    }

    addProfiling(merge, "snappyHexMesh::refine::erode");
    Info<< nl
        << "Erode non-manifold zone faces" << nl
        << "-----------------------------" << nl
        << endl;

    auto& mesh = meshRefiner_.mesh();
    auto& fzs = mesh.faceZones();

    // Detect a faceZone face that:
    // - sits on the edge of the faceZone (i.e. has an open edge)
    // - has a non-manifold edge
    //
    // Limitations:
    // - analyses across all zones - does not support overlapping zones
    // - requires zones to be consistent across processor patches

    bitSet isZoneFace(mesh.nFaces());
    for (const auto& fz : fzs)
    {
        isZoneFace.set(fz.addressing());
    }
    // Unmark non-owner side of processor patches so they don't get double
    // counted.
    const auto& pbm = mesh.boundaryMesh();
    for (label patchi = pbm.nNonProcessor(); patchi < pbm.size(); patchi++)
    {
        if (!refCast<const processorPolyPatch>(pbm[patchi]).owner())
        {
            isZoneFace.unset(pbm[patchi].range());
        }
    }


    // TBD: replace with mesh.globalData().globalMeshFaceAddr()
    const globalIndex globalFaces(mesh.nFaces());

    label nChanged = 0;

    while (true)
    {
        const labelList meshFaces(isZoneFace.sortedToc());

        const indirectPrimitivePatch pp
        (
            IndirectList<face>
            (
                mesh.faces(),
                meshFaces
            ),
            mesh.points()
        );

        const labelListList edgeGlobalFaces
        (
            addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                pp
            )
        );

        const label nOldChanged = nChanged;
        forAll(pp, patchFacei)
        {
            // Detect at least one open edge and one non-manifold edge
            label nOpen = 0;
            label nNonManif = 0;
            for (const label edgei : pp.faceEdges()[patchFacei])
            {
                if (edgeGlobalFaces[edgei].size() == 1)
                {
                    nOpen++;
                }
                else if (edgeGlobalFaces[edgei].size() > 2)
                {
                    nNonManif++;
                }
            }
            if (nNonManif > 0 && nOpen > 0)
            {
                isZoneFace.unset(meshFaces[patchFacei]);
                nChanged++;
            }
        }

        if (returnReduce((nChanged-nOldChanged), sumOp<label>()) == 0)
        {
            break;
        }
    }


    reduce(nChanged, sumOp<label>());

    if (nChanged)
    {
        // Filter out eroded faceZone faces
        fzs.clearAddressing();

        for (auto& fz : fzs)
        {
            DynamicList<label> addressing(fz.size());
            DynamicList<bool> flipMap(fz.size());

            forAll(fz.addressing(), i)
            {
                if (isZoneFace[fz.addressing()[i]])
                {
                    addressing.append(fz.addressing()[i]);
                    flipMap.append(fz.flipMap()[i]);
                }
            }

            //Pout<< "** on faceZone:" << fz.name()
            //  << " from:" << fz.size()
            //    << " to:" << addressing.size() << endl;

            fz.resetAddressing(addressing, flipMap);
        }
    }

    Info<< "Eroded " << nChanged
        << " free-standing zone faces in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


void Foam::snappyRefineDriver::addFaceZones
(
    meshRefinement& meshRefiner,
    const refinementParameters& refineParams,
    const HashTable<Pair<word>>& faceZoneToPatches
)
{
    if (faceZoneToPatches.size())
    {
        Info<< nl
            << "Adding patches for face zones" << nl
            << "-----------------------------" << nl
            << endl;

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Name"
            << setw(30) << "FaceZone"
            << setw(10) << "FaceType"
            << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "----"
            << setw(30) << "--------"
            << setw(10) << "--------"
            << endl;

        const polyMesh& mesh = meshRefiner.mesh();

        // Add patches for added inter-region faceZones
        forAllConstIters(faceZoneToPatches, iter)
        {
            const word& fzName = iter.key();
            const Pair<word>& patchNames = iter.val();

            // Get any user-defined faceZone data
            surfaceZonesInfo::faceZoneType fzType;
            dictionary patchInfo = refineParams.getZoneInfo(fzName, fzType);

            const word& masterName = fzName;
            //const word slaveName = fzName + "_slave";
            //const word slaveName = czNames.second()+"_to_"+czNames.first();
            const word& slaveName = patchNames.second();

            label mpi = meshRefiner.addMeshedPatch(masterName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << mpi
                << setw(20) << mesh.boundaryMesh()[mpi].type()
                << setw(30) << masterName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;


            label sli = meshRefiner.addMeshedPatch(slaveName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << sli
                << setw(20) << mesh.boundaryMesh()[sli].type()
                << setw(30) << slaveName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;

            meshRefiner.addFaceZone(fzName, masterName, slaveName, fzType);
        }

        Info<< endl;
    }
}


void Foam::snappyRefineDriver::mergePatchFaces
(
    const meshRefinement::FaceMergeType mergeType,
    const refinementParameters& refineParams,
    const dictionary& motionDict
)
{
    if (dryRun_)
    {
        return;
    }

    addProfiling(merge, "snappyHexMesh::refine::merge");
    Info<< nl
        << "Merge refined boundary faces" << nl
        << "----------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if
    (
        mergeType == meshRefinement::FaceMergeType::GEOMETRIC
     || mergeType == meshRefinement::FaceMergeType::IGNOREPATCH
    )
    {
        meshRefiner_.mergePatchFacesUndo
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            meshRefiner_.meshedPatches(),
            motionDict,
            labelList(mesh.nFaces(), -1),
            mergeType
        );
    }
    else
    {
        // Still merge refined boundary faces if all four are on same patch
        meshRefiner_.mergePatchFaces
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            4,          // only merge faces split into 4
            meshRefiner_.meshedPatches(),
            meshRefinement::FaceMergeType::GEOMETRIC // no merge across patches
        );
    }

    if (debug)
    {
        meshRefiner_.checkData();
    }

    meshRefiner_.mergeEdgesUndo(Foam::cos(degToRad(45.0)), motionDict);

    if (debug)
    {
        meshRefiner_.checkData();
    }
}


void Foam::snappyRefineDriver::deleteSmallRegions
(
    const refinementParameters& refineParams
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const cellZoneMesh& czm = mesh.cellZones();

    //const labelList zoneIDs
    //(
    //    meshRefiner_.getZones
    //    (
    //        List<surfaceZonesInfo::faceZoneType> fzTypes
    //        ({
    //            surfaceZonesInfo::BAFFLE,
    //            surfaceZonesInfo::BOUNDARY,
    //        });
    //    )
    //);
    const labelList zoneIDs(identity(mesh.faceZones().size()));

    // Get faceZone and patch(es) per face (or -1 if face not on faceZone)
    labelList faceZoneID;
    labelList ownPatch;
    labelList neiPatch;
    labelList nB;       // local count of faces per faceZone
    meshRefiner_.getZoneFaces(zoneIDs, faceZoneID, ownPatch, neiPatch, nB);


    // Mark all faces on outside of zones. Note: assumes that faceZones
    // are consistent with the outside of cellZones ...

    boolList isBlockedFace(mesh.nFaces(), false);
    meshRefiner_.selectSeparatedCoupledFaces(isBlockedFace);

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1)
        {
            isBlockedFace[facei] = true;
        }
    }

    // Map from cell to zone. Not that background cells are not -1 but
    // cellZones.size()
    labelList cellToZone(mesh.nCells(), czm.size());
    for (const auto& cz : czm)
    {
        UIndirectList<label>(cellToZone, cz) = cz.index();
    }

    // Walk to split into regions
    const regionSplit cellRegion(mesh, isBlockedFace);

    // Count number of cells per zone and per region
    labelList nCellsPerRegion(cellRegion.nRegions(), 0);
    labelList regionToZone(cellRegion.nRegions(), -2);
    labelList nCellsPerZone(czm.size()+1, 0);
    forAll(cellRegion, celli)
    {
        const label regioni = cellRegion[celli];
        const label zonei = cellToZone[celli];

        // Zone for this region
        regionToZone[regioni] = zonei;

        nCellsPerRegion[regioni]++;
        nCellsPerZone[zonei]++;
    }
    Pstream::listCombineReduce(nCellsPerRegion, plusEqOp<label>());
    Pstream::listCombineReduce(regionToZone, maxEqOp<label>());
    Pstream::listCombineReduce(nCellsPerZone, plusEqOp<label>());


    // Mark small regions. Note that all processors have the same information
    // so should take the same decision. Alternative is to do master only
    // and scatter the decision.
    forAll(nCellsPerRegion, regioni)
    {
        const label zonei = regionToZone[regioni];
        label& nRegionCells = nCellsPerRegion[regioni];

        if
        (
            nRegionCells < refineParams.minCellFraction()*nCellsPerZone[zonei]
         || nRegionCells < refineParams.nMinCells()
        )
        {
            Info<< "Deleting region " << regioni
                << " (size " << nRegionCells
                << ") of zone size " << nCellsPerZone[zonei]
                << endl;

            // Mark region to be deleted. 0 size (= global) should never
            // occur.
            nRegionCells = 0;
        }
    }

    DynamicList<label> cellsToRemove(mesh.nCells()/128);
    forAll(cellRegion, celli)
    {
        if (nCellsPerRegion[cellRegion[celli]] == 0)
        {
            cellsToRemove.append(celli);
        }
    }
    const label nTotCellsToRemove = returnReduce
    (
        cellsToRemove.size(),
        sumOp<label>()
    );
    if (nTotCellsToRemove > 0)
    {
        Info<< "Deleting " << nTotCellsToRemove
            << " cells in small regions" << endl;

        removeCells cellRemover(mesh);

        cellsToRemove.shrink();
        const labelList exposedFaces
        (
            cellRemover.getExposedFaces(cellsToRemove)
        );
        const labelList exposedPatch
        (
            UIndirectList<label>(ownPatch, exposedFaces)
        );
        (void)meshRefiner_.doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatch,
            cellRemover
        );
    }
}


void Foam::snappyRefineDriver::doRefine
(
    const dictionary& refineDict,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool prepareForSnapping,
    const meshRefinement::FaceMergeType mergeType,
    const dictionary& motionDict
)
{
    addProfiling(refine, "snappyHexMesh::refine");
    Info<< nl
        << "Refinement phase" << nl
        << "----------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();


    // Check that all the keep points are inside the mesh.
    refineParams.findCells(true, mesh, refineParams.locationsInMesh());

    // Check that all the keep points are inside the mesh.
    if (dryRun_)
    {
        refineParams.findCells(true, mesh, refineParams.locationsOutsideMesh());
    }

    // Estimate cell sizes
    if (dryRun_)
    {
        snappyVoxelMeshDriver voxelDriver
        (
            meshRefiner_,
            globalToMasterPatch_,
            globalToSlavePatch_
        );
        voxelDriver.doRefine(refineParams);
    }


    // Refine around feature edges
    featureEdgeRefine
    (
        refineParams,
        100,    // maxIter
        0       // min cells to refine
    );


    // Initial automatic gap-level refinement: gaps smaller than the cell size
    if
    (
        max(meshRefiner_.surfaces().maxGapLevel()) > 0
     || max(meshRefiner_.shells().maxGapLevel()) > 0
     || max(meshRefiner_.surfaces().maxCurvatureLevel()) > 0
    )
    {
        // In case we use automatic gap level refinement do some pre-refinement
        // (fast) since it is so slow.

        // Refine based on surface
        surfaceOnlyRefine
        (
            refineParams,
            20,         // maxIter
            labelMax    // no leak detection yet since all in same cell
        );

        // Refine cells that contain a gap
        smallFeatureRefine
        (
            refineParams,
            100     // maxIter
        );
    }


    // Refine based on surface
    surfaceOnlyRefine
    (
        refineParams,
        100,    // maxIter
        0       // Iteration to start leak detection and closure
    );

    // Pass1 of automatic gap-level refinement: surface-intersected cells
    // in narrow gaps. Done early so we can remove the inside
    gapOnlyRefine
    (
        refineParams,
        100     // maxIter
    );

    // Remove cells inbetween two surfaces
    surfaceProximityBlock
    (
        refineParams,
        1  //100     // maxIter
    );

    // Remove cells (a certain distance) beyond surface intersections
    removeInsideCells
    (
        refineParams,
        1       // nBufferLayers
    );

    // Pass2 of automatic gap-level refinement: all cells in gaps
    bigGapOnlyRefine
    (
        refineParams,
        true,   // spreadGapSize
        100     // maxIter
    );

    // Internal mesh refinement
    shellRefine
    (
        refineParams,
        100    // maxIter
    );

    // Remove any extra cells from limitRegion with level -1, without
    // adding any buffer layer this time
    removeInsideCells
    (
        refineParams,
        0       // nBufferLayers
    );

    // Refine any hexes with 5 or 6 faces refined to make smooth edges
    danglingCellRefine
    (
        refineParams,
        21,     // 1 coarse face + 5 refined faces
        100     // maxIter
    );
    danglingCellRefine
    (
        refineParams,
        24,     // 0 coarse faces + 6 refined faces
        100     // maxIter
    );

    // Refine any cells with differing refinement level on either side
    refinementInterfaceRefine
    (
        refineParams,
        10      // maxIter
    );

    // Directional shell refinement
    directionalShellRefine
    (
        refineParams,
        100    // maxIter
    );

    // Block gaps (always, ignore surface leakLevel)
    if (refineParams.locationsOutsideMesh().size())
    {
        // For now: only check leaks on meshed surfaces. The problem is that
        // blockLeakFaces always generates baffles any not just faceZones ...
        const labelList unnamedSurfaces
        (
            surfaceZonesInfo::getUnnamedSurfaces
            (
                meshRefiner_.surfaces().surfZones()
            )
        );
        if (unnamedSurfaces.size() && !dryRun_)
        {
            meshRefiner_.blockLeakFaces
            (
                globalToMasterPatch_,
                globalToSlavePatch_,
                refineParams.locationsInMesh(),
                refineParams.zonesInMesh(),
                refineParams.locationsOutsideMesh(),
                unnamedSurfaces
            );
        }
    }

    if
    (
        max(meshRefiner_.shells().nSmoothExpansion()) > 0
     || max(meshRefiner_.shells().nSmoothPosition()) > 0
    )
    {
        directionalSmooth(refineParams);
    }


    // Introduce baffles at surface intersections. Remove sections unreachable
    // from keepPoint.
    baffleAndSplitMesh
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );


    //- Commented out for now since causes zoning errors (sigsegv) on
    //  case with faceZones. TBD.
    //// Refine any cells are point/edge connected to the boundary and have a
    //// lower refinement level than the neighbouring cells
    //boundaryRefinementInterfaceRefine
    //(
    //    refineParams,
    //    10      // maxIter
    //);


    // Mesh is at its finest. Do optional zoning (cellZones and faceZones)
    wordPairHashTable zonesToFaceZone;
    zonify(refineParams, zonesToFaceZone);

    // Create pairs of patches for faceZones
    {
        HashTable<Pair<word>> faceZoneToPatches(zonesToFaceZone.size());

        //    Note: zonesToFaceZone contains the same data on different
        //          processors but in different order. We could sort the
        //          contents but instead just loop in sortedToc order.
        List<Pair<word>> czs(zonesToFaceZone.sortedToc());

        forAll(czs, i)
        {
            const Pair<word>& czNames = czs[i];
            const word& fzName = zonesToFaceZone[czNames];

            const word& masterName = fzName;
            const word slaveName = czNames.second() + "_to_" + czNames.first();
            Pair<word> patches(masterName, slaveName);
            faceZoneToPatches.insert(fzName, patches);
        }
        addFaceZones(meshRefiner_, refineParams, faceZoneToPatches);
    }

    // Pull baffles apart
    splitAndMergeBaffles
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );


    // Do something about cells with refined faces on the boundary
    if (prepareForSnapping)
    {
        const auto mt = meshRefiner_.meshType();
        if
        (
            meshRefiner_.mesh().faceZones().size()
         && (
                mt == meshRefinement::CASTELLATEDBUFFERLAYER
             || mt == meshRefinement::CASTELLATEDBUFFERLAYER2
            )
        )
        {
            erodeNonManifoldZoneFaces(refineParams);
        }

        mergePatchFaces(mergeType, refineParams, motionDict);
    }


    if (refineParams.minCellFraction() > 0 || refineParams.nMinCells() > 0)
    {
        // Some small disconnected bits of mesh might remain since at
        // this point faceZones have not been converted into e.g. baffles.
        // We don't know whether e.g. the baffles are reset to be cyclicAMI
        // thus reconnecting. For now check if there are any particularly
        // small regions.
        deleteSmallRegions(refineParams);
    }


    if (!dryRun_ && Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;


        /*
        const bool hasBufferLayer
        (
            (meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER)
         || (meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER2)
        );
        */
        const bool hasBufferLayer = false;

        labelList singleProcPoints;

        /*
        if (hasBufferLayer)
        {
            //- Needed since buffer layer addition did not handle
            //- inter-processor extrusion. Fixed now. This code can be removed.

            // Pick up all points on surfaces with specified buffer layers
            const labelList& surfIndex = meshRefiner_.surfaceIndex();
            const auto& addLayers = meshRefiner_.surfaces().addBufferLayers();

            bitSet isSelected(mesh.nPoints());
            forAll(surfIndex, facei)
            {
                const label surfi = surfIndex[facei];
                if (surfIndex[facei] != -1)
                {
                    const label globalRegioni =
                        meshRefiner_.surfaces().globalRegion(surfi, 0);
                    if (addLayers[globalRegioni])
                    {
                        isSelected.set(mesh.faces()[facei]);
                    }
                }
            }
            syncTools::syncPointList
            (
                mesh,
                isSelected,
                orEqOp<unsigned int>(),
                0
            );
            singleProcPoints = isSelected.sortedToc();
        }
        */

        // Do final balancing. Keep zoned faces on one processor since the
        // snap phase will convert them to baffles and this only works for
        // internal faces.
        meshRefiner_.balance
        (
            true,                           // keepZoneFaces
            hasBufferLayer,                 // keepBaffles
            singleProcPoints,               // keepZonePoints
            scalarField(mesh.nCells(), 1),  // cellWeights
            decomposer_,
            distributor_
        );
    }
}


// ************************************************************************* //
