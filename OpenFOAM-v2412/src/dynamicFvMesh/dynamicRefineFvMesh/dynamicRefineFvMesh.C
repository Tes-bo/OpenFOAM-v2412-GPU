/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "dynamicRefineFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "sigFpe.H"
#include "cellSet.H"
#include "HashOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineFvMesh, IOobject);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineFvMesh, doInit);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicRefineFvMesh::calculateProtectedCells
(
    bitSet& unrefineableCell
) const
{
    if (protectedCell_.empty())
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    unrefineableCell = protectedCell_;

    // Get neighbouring cell level
    labelList neiLevel(nBoundaryFaces());

    for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
    {
        neiLevel[facei-nInternalFaces()] = cellLevel[faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);


    bitSet seedFace;

    while (true)
    {
        // Pick up faces on border of protected cells
        seedFace.reset();
        seedFace.resize(nFaces());

        for (label facei = 0; facei < nInternalFaces(); ++facei)
        {
            const label own = faceOwner()[facei];
            const label nei = faceNeighbour()[facei];

            if
            (
                // Protected owner
                (
                    unrefineableCell.test(own)
                 && (cellLevel[nei] > cellLevel[own])
                )
             ||
                // Protected neighbour
                (
                    unrefineableCell.test(nei)
                 && (cellLevel[own] > cellLevel[nei])
                )
            )
            {
                seedFace.set(facei);
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            const label own = faceOwner()[facei];

            if
            (
                // Protected owner
                (
                    unrefineableCell.test(own)
                 && (neiLevel[facei-nInternalFaces()] > cellLevel[own])
                )
            )
            {
                seedFace.set(facei);
            }
        }

        syncTools::syncFaceList(*this, seedFace, orEqOp<unsigned int>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < nInternalFaces(); ++facei)
        {
            if (seedFace.test(facei))
            {
                if (unrefineableCell.set(faceOwner()[facei]))
                {
                    hasExtended = true;
                }
                if (unrefineableCell.set(faceNeighbour()[facei]))
                {
                    hasExtended = true;
                }
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
        {
            if (seedFace.test(facei))
            {
                const label own = faceOwner()[facei];

                if (unrefineableCell.set(own))
                {
                    hasExtended = true;
                }
            }
        }

        if (!returnReduceOr(hasExtended))
        {
            break;
        }
    }
}


void Foam::dynamicRefineFvMesh::readDict()
{
    const dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).optionalSubDict(typeName + "Coeffs")
    );

    auto fluxVelocities = refineDict.get<List<Pair<word>>>("correctFluxes");

    // Rework into hashtable.
    correctFluxes_.reserve(fluxVelocities.size());
    for (const auto& pr : fluxVelocities)
    {
        correctFluxes_.insert(pr.first(), pr.second());
    }

    refineDict.readEntry("dumpLevel", dumpLevel_);
}


void Foam::dynamicRefineFvMesh::mapFields(const mapPolyMesh& mpm)
{
    //dynamicFvMesh::mapFields(mpm);
    dynamicMotionSolverListFvMesh::mapFields(mpm);


    // Correct old-time volumes for refined/unrefined cells. We know at this
    // point that the points have not moved and the cells have only been split
    // or merged. We hope that dynamicMotionSolverListFvMesh::mapFields
    // does not use old-time volumes ...
    {
        const labelList& cellMap = mpm.cellMap();
        const labelList& reverseCellMap = mpm.reverseCellMap();

        // Each split cell becomes original + 7 additional
        labelList nSubCells(mpm.nOldCells(), 0);

        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && reverseCellMap[oldCelli] >= 0)
            {
                // Found master cell.
                nSubCells[oldCelli]++;
            }
        }

        // Start off from mapped old volumes
        scalarField correctedV0(this->V0());

        // Correct any split cells
        const auto& V = this->V();
        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && nSubCells[oldCelli] == 8)
            {
                // Found split cell. Use current volume instead of mapped
                // old one
                correctedV0[celli] = V[celli];
            }
        }

        const auto& cellsFromCells = mpm.cellsFromCellsMap();
        for (const auto& s : cellsFromCells)
        {
            // Reset unsplit cell
            const label celli = s.index();
            correctedV0[celli] = V[celli];
            //? Or sum up old volumes?
            //const labelList& oldCells = s.masterObjects();
            //for (const label oldCelli : oldCells)
            //{
            //}
        }

        setV0().field() = correctedV0;
    }


    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = mpm.faceMap();
        const labelList& reverseFaceMap = mpm.reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        // Estimate number of faces created

        bitSet masterFaces(nFaces());

        forAll(faceMap, facei)
        {
            const label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                const label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei << endl
                        << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.set(masterFacei);
                }
            }
        }

        if (debug)
        {
            Pout<< "Found " << masterFaces.count() << " split faces " << endl;
        }

        UPtrList<surfaceScalarField> fluxes
        (
            this->objectRegistry::sorted<surfaceScalarField>()
        );

        // Remove surfaceInterpolation to allow re-calculation on demand
        // This could be done in fvMesh::updateMesh but some dynamicFvMesh
        // might need the old interpolation fields (weights, etc).
        surfaceInterpolation::clearOut();

        for (surfaceScalarField& phi : fluxes)
        {
            const word& UName = correctFluxes_.lookup(phi.name(), word::null);

            if (UName.empty())
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << phi.name()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << phi.name() << " none) to suppress this warning."
                    << endl;
                continue;
            }
            if (UName == "none")
            {
                continue;
            }
            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << phi.name()
                    << " to NaN" << endl;

                sigFpe::fillNan(phi.primitiveFieldRef());
                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << phi.name()
                    << " using interpolated flux " << UName
                    << endl;
            }

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < nInternalFaces(); ++facei)
            {
                const label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    phi[facei] = phiU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    phi[facei] = phiU[facei];
                }
            }

            // Recalculate new boundary faces.
            auto& phiBf = phi.boundaryFieldRef();

            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    const label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    ++facei;
                }
            }

            // Update master faces
            for (const label facei : masterFaces)
            {
                if (isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    const label patchi = boundaryMesh().whichPatch(facei);
                    const label i = facei - boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }

    // Correct the flux for injected faces - these are the faces which have
    // no correspondence to the old mesh (i.e. added without a masterFace, edge
    // or point). An example is the internal faces from hexRef8.
    {
        const labelList& faceMap = mpm.faceMap();

        mapNewInternalFaces<scalar>(this->Sf(), this->magSf(), faceMap);
        mapNewInternalFaces<vector>(this->Sf(), this->magSf(), faceMap);

        // No oriented fields of more complex type
        mapNewInternalFaces<sphericalTensor>(faceMap);
        mapNewInternalFaces<symmTensor>(faceMap);
        mapNewInternalFaces<tensor>(faceMap);
    }
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);

    // Clear moving flag. This is currently required since geometry calculation
    // might get triggered when doing processor patches.
    // (TBD: should be in changeMesh if no inflation?)
    moving(false);
    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < nInternalFaces(); ++facei)
        {
            const label oldFacei = map().faceMap()[facei];

            if (oldFacei >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    //    // Remove the stored tet base points
    //    tetBasePtIsPtr_.clear();
    //    // Remove the cell tree
    //    cellTreePtr_.clear();

    // Update fields
    updateMesh(*map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */



    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(*map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        bitSet newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            const label oldCelli = map().cellMap()[celli];
            if (protectedCell_.test(oldCelli))
            {
                newProtectedCell.set(celli);
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList());

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineFvMesh::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setUnrefinement(splitPoints, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        for (const label pointi : splitPoints)
        {
            const labelList& pEdges = pointEdges()[pointi];

            for (const label edgei : pEdges)
            {
                const label otherPointi = edges()[edgei].otherVertex(pointi);

                const labelList& pFaces = pointFaces()[otherPointi];

                for (const label facei : pFaces)
                {
                    faceToSplitPoint.insert(facei, otherPointi);
                }
            }
        }
    }


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);

    // Clear moving flag. This is currently required since geometry calculation
    // might get triggered when doing processor patches.
    // (TBD: should be in changeMesh if no inflation?)
    moving(false);
    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(*map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        UPtrList<surfaceScalarField> fluxes
        (
            this->objectRegistry::sorted<surfaceScalarField>()
        );

        for (surfaceScalarField& phi : fluxes)
        {
            const word& UName = correctFluxes_.lookup(phi.name(), word::null);

            if (UName.empty())
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << phi.name()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << phi.name() << " none) to suppress this warning."
                    << endl;
                continue;
            }
            if (UName == "none")
            {
                continue;
            }

            DebugInfo
                << "Mapping flux " << phi.name()
                << " using interpolated flux " << UName
                << endl;

            auto& phiBf = phi.boundaryFieldRef();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIters(faceToSplitPoint, iter)
            {
                const label oldFacei = iter.key();
                const label oldPointi = iter.val();

                if (reversePointMap[oldPointi] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    const label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (isInternalFace(facei))
                        {
                            phi[facei] = phiU[facei];
                        }
                        else
                        {
                            label patchi = boundaryMesh().whichPatch(facei);
                            label i = facei - boundaryMesh()[patchi].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchi];
                            fvsPatchScalarField& patchPhi = phiBf[patchi];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(*map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        bitSet newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            const label oldCelli = map().cellMap()[celli];
            if (protectedCell_.test(oldCelli))
            {
                newProtectedCell.set(celli);
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList());

    return map;
}


Foam::scalarField
Foam::dynamicRefineFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        for (const label celli : pCells)
        {
            vFld[celli] = max(vFld[celli], pFld[pointi]);
        }
    }
    return vFld;
}


Foam::scalarField
Foam::dynamicRefineFvMesh::maxCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        for (const label celli : pCells)
        {
            pFld[pointi] = max(pFld[pointi], vFld[celli]);
        }
    }
    return pFld;
}


Foam::scalarField
Foam::dynamicRefineFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar sum = 0.0;
        for (const label celli : pCells)
        {
            sum += vFld[celli];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


Foam::scalarField Foam::dynamicRefineFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), scalar(-1));

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::dynamicRefineFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    bitSet& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCell.set(celli);
        }
    }
}


Foam::labelList Foam::dynamicRefineFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const bitSet& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    bitSet unrefineableCell;
    calculateProtectedCells(unrefineableCell);

    // Count current selection
    label nLocalCandidates = candidateCell.count();
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    if (nCandidates < nTotToRefine)
    {
        for (const label celli : candidateCell)
        {
            if
            (
                (!unrefineableCell.test(celli))
             && cellLevel[celli] < maxRefinement
            )
            {
                candidates.append(celli);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; ++level)
        {
            for (const label celli : candidateCell)
            {
                if
                (
                    (!unrefineableCell.test(celli))
                 && cellLevel[celli] == level
                )
                {
                    candidates.append(celli);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::dynamicRefineFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const bitSet& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());


    const labelListList& pointCells = this->pointCells();

    // If we have any protected cells make sure they also are not being
    // unrefined

    bitSet protectedPoint(nPoints());

    if (protectedCell_.size())
    {
        // Get all points on a protected cell
        forAll(pointCells, pointi)
        {
            for (const label celli : pointCells[pointi])
            {
                if (protectedCell_.test(celli))
                {
                    protectedPoint.set(pointi);
                    break;
                }
            }
        }

        syncTools::syncPointList
        (
            *this,
            protectedPoint,
            orEqOp<unsigned int>(),
            0u
        );

        DebugInfo<< "From "
            << returnReduce(protectedCell_.count(), sumOp<label>())
            << " protected cells found "
            << returnReduce(protectedPoint.count(), sumOp<label>())
            << " protected points." << endl;
    }


    DynamicList<label> newSplitPoints(splitPoints.size());

    for (const label pointi : splitPoints)
    {
        if (!protectedPoint[pointi] && pFld[pointi] < unrefineLevel)
        {
            // Check that all cells are not marked
            bool hasMarked = false;

            for (const label celli : pointCells[pointi])
            {
                if (markedCell.test(celli))
                {
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked)
            {
                newSplitPoints.append(pointi);
            }
        }
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void Foam::dynamicRefineFvMesh::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace
    for (label facei = 0; facei < nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(faceOwner()[facei]);
            markedCell.set(faceNeighbour()[facei]);
        }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(faceOwner()[facei]);
        }
    }
}


void Foam::dynamicRefineFvMesh::checkEightAnchorPoints
(
    bitSet& protectedCell
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(nCells(), Zero);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = pointCells(pointi);

        for (const label celli : pCells)
        {
            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[celli] == 8)
                {
                    protectedCell.set(celli);
                }

                if (!protectedCell.test(celli))
                {
                    ++nAnchorPoints[celli];
                }
            }
        }
    }


    forAll(protectedCell, celli)
    {
        if (nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineFvMesh::dynamicRefineFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    //dynamicFvMesh(io, doInit),
    dynamicMotionSolverListFvMesh(io, doInit),
    meshCutter_(*this)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicRefineFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        //dynamicFvMesh::init(doInit);
        // Note: allow zero motion solvers
        dynamicMotionSolverListFvMesh::init(doInit, false);
    }

    protectedCell_.setSize(nCells());
    nRefinementIterations_ = 0;
    dumpLevel_ = false;

    // Read static part of dictionary
    readDict();


    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(nCells(), Zero);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        for (const label celli : pCells)
        {
            if (!protectedCell_.test(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    ++nAnchors[celli];

                    if (nAnchors[celli] > 8)
                    {
                        protectedCell_.set(celli);
                    }
                }
            }
        }
    }


    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(nFaces());

        for (label facei = 0; facei < nInternalFaces(); ++facei)
        {
            neiLevel[facei] = cellLevel[faceNeighbour()[facei]];
        }
        for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
        {
            neiLevel[facei] = cellLevel[faceOwner()[facei]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        bitSet protectedFace(nFaces());

        forAll(faceOwner(), facei)
        {
            const label faceLevel = max
            (
                cellLevel[faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = faces()[facei];

            label nAnchors = 0;

            for (const label pointi : f)
            {
                if (pointLevel[pointi] <= faceLevel)
                {
                    ++nAnchors;

                    if (nAnchors > 4)
                    {
                        protectedFace.set(facei);
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(*this, protectedFace, orEqOp<unsigned int>());

        for (label facei = 0; facei < nInternalFaces(); ++facei)
        {
            if (protectedFace.test(facei))
            {
                protectedCell_.set(faceOwner()[facei]);
                protectedCell_.set(faceNeighbour()[facei]);
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
        {
            if (protectedFace.test(facei))
            {
                protectedCell_.set(faceOwner()[facei]);
            }
        }

        // Also protect any cells that are less than hex
        forAll(cells(), celli)
        {
            const cell& cFaces = cells()[celli];

            if (cFaces.size() < 6)
            {
                protectedCell_.set(celli);
            }
            else
            {
                for (const label cfacei : cFaces)
                {
                    if (faces()[cfacei].size() < 4)
                    {
                        protectedCell_.set(celli);
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCell_);
    }

    if (!returnReduceOr(protectedCell_.any()))
    {
        protectedCell_.clear();
    }
    else
    {
        cellSet protectedCells
        (
            *this,
            "protectedCells",
            HashSetOps::used(protectedCell_)
        );

        Info<< "Detected "
            << returnReduce(protectedCells.size(), sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        protectedCells.write();
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineFvMesh::updateTopology()
{
    // Re-read dictionary. Chosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).optionalSubDict(typeName + "Coeffs")
    );

    const label refineInterval = refineDict.get<label>("refineInterval");

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }


    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        const label maxCells = refineDict.get<label>("maxCells");

        if (maxCells <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const label maxRefinement = refineDict.get<label>("maxRefinement");

        if (maxRefinement <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const word fieldName(refineDict.get<word>("field"));

        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        const scalar lowerRefineLevel =
            refineDict.get<scalar>("lowerRefineLevel");
        const scalar upperRefineLevel =
            refineDict.get<scalar>("upperRefineLevel");
        const scalar unrefineLevel = refineDict.getOrDefault<scalar>
        (
            "unrefineLevel",
            GREAT
        );
        const label nBufferLayers = refineDict.get<label>("nBufferLayers");

        // Cells marked for refinement or otherwise protected from unrefinement.
        bitSet refineCell(nCells());

        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            vFld,
            refineCell
        );

        if (globalData().nTotalCells() < maxCells)
        {
            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
                    maxRefinement,
                    refineCell
                )
            );

            if (returnReduceOr(cellsToRefine.size()))
            {
                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    bitSet newRefineCell(cellMap.size());

                    forAll(cellMap, celli)
                    {
                        const label oldCelli = cellMap[celli];

                        if
                        (
                            (oldCelli < 0)
                         || (reverseCellMap[oldCelli] != celli)
                         || (refineCell.test(oldCelli))
                        )
                        {
                            newRefineCell.set(celli);
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                // Extend with a buffer layer to prevent neighbouring points
                // being unrefined.
                for (label i = 0; i < nBufferLayers; ++i)
                {
                    extendMarkedCells(refineCell);
                }

                hasChanged = true;
            }
        }


        {
            // Select unrefineable points that are not marked in refineCell
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(vFld)
                )
            );

            if (returnReduceOr(pointsToUnrefine.size()))
            {
                // Refine/update mesh
                unrefine(pointsToUnrefine);

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occasionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
    }

    topoChanging(hasChanged);
    if (hasChanged)
    {
        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }

    return hasChanged;
}


bool Foam::dynamicRefineFvMesh::update()
{
    bool hasChanged = updateTopology();
    // Do any mesh motion (resets mesh.moving() if it does any mesh motion)
    hasChanged = dynamicMotionSolverListFvMesh::update() && hasChanged;

    return hasChanged;
}


bool Foam::dynamicRefineFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
    (
        //dynamicFvMesh::writeObject(streamOpt, writeOnProc)
        dynamicMotionSolverListFvMesh::writeObject(streamOpt, writeOnProc)
     && meshCutter_.write(writeOnProc)
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimensionedScalar(dimless, Zero)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, celli)
        {
            scalarCellLevel[celli] = cellLevel[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
