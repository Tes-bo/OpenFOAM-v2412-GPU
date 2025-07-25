/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "polyTopoChange.H"
#include "polyMesh.H"
#include "polyAddPoint.H"
#include "polyModifyPoint.H"
#include "polyRemovePoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyAddCell.H"
#include "polyModifyCell.H"
#include "polyRemoveCell.H"
#include "CircularBuffer.H"
#include "CompactListList.H"
#include "objectMap.H"
#include "processorPolyPatch.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTopoChange, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Renumber with special handling for merged items (marked with <-1)
void Foam::polyTopoChange::renumberReverseMap
(
    const labelUList& oldToNew,
    DynamicList<label>& elems
)
{
    forAll(elems, elemI)
    {
        const label val = elems[elemI];

        if (val >= 0)
        {
            elems[elemI] = oldToNew[val];
        }
        else if (val < -1)
        {
            const label mergedVal = -val-2;
            elems[elemI] = -oldToNew[mergedVal]-2;
        }
    }
}


void Foam::polyTopoChange::renumber
(
    const labelUList& oldToNew,
    labelHashSet& labels
)
{
    labelHashSet newSet(labels.capacity());

    for (const label val : labels)
    {
        const label newVal = oldToNew[val];

        if (newVal >= 0)
        {
            newSet.insert(newVal);
        }
    }

    labels.transfer(newSet);
}


// Renumber and remove -1 elements.
void Foam::polyTopoChange::renumberCompact
(
    const labelUList& oldToNew,
    labelList& elems
)
{
    label nElem = 0;

    for (const label val : elems)
    {
        const label newVal = oldToNew[val];

        if (newVal != -1)
        {
            elems[nElem++] = newVal;
        }
    }
    elems.setSize(nElem);
}


void Foam::polyTopoChange::countMap
(
    const labelUList& map,
    const labelUList& reverseMap,
    label& nAdd,
    label& nInflate,
    label& nMerge,
    label& nRemove
)
{
    nAdd = 0;
    nInflate = 0;
    nMerge = 0;
    nRemove = 0;

    forAll(map, newCelli)
    {
        const label oldCelli = map[newCelli];

        if (oldCelli >= 0)
        {
            if
            (
                oldCelli < reverseMap.size()
             && reverseMap[oldCelli] == newCelli
            )
            {
                // unchanged
            }
            else
            {
                // Added (from another cell v.s. inflated from face/point)
                nAdd++;
            }
        }
        else if (oldCelli == -1)
        {
            // Created from nothing
            nInflate++;
        }
        else
        {
            FatalErrorInFunction
                << " new:" << newCelli << abort(FatalError);
        }
    }

    forAll(reverseMap, oldCelli)
    {
        const label newCelli = reverseMap[oldCelli];

        if (newCelli >= 0)
        {
            // unchanged
        }
        else if (newCelli == -1)
        {
            // removed
            nRemove++;
        }
        else
        {
            // merged into -newCelli-2
            nMerge++;
        }
    }
}


void Foam::polyTopoChange::writeMeshStats(const polyMesh& mesh, Ostream& os)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    forAll(patches, patchi)
    {
        patchSizes[patchi] = patches[patchi].size();
        patchStarts[patchi] = patches[patchi].start();
    }

    const auto& czs = mesh.cellZones();
    labelList cellZoneSizes(czs.size(), 0);
    for (const auto& cz : czs)
    {
        cellZoneSizes[cz.index()] = cz.size();
    }
    const auto& fzs = mesh.faceZones();
    labelList faceZoneSizes(fzs.size(), 0);
    for (const auto& fz : fzs)
    {
        faceZoneSizes[fz.index()] = fz.size();
    }
    const auto& pzs = mesh.pointZones();
    labelList pointZoneSizes(pzs.size(), 0);
    for (const auto& pz : pzs)
    {
        pointZoneSizes[pz.index()] = pz.size();
    }

    os  << "    Points      : " << mesh.nPoints() << nl
        << "    Faces       : " << mesh.nFaces() << nl
        << "    Cells       : " << mesh.nCells() << nl
        << "    PatchSizes  : " << flatOutput(patchSizes) << nl
        << "    PatchStarts : " << flatOutput(patchStarts) << nl
        << "    cZoneSizes  : " << flatOutput(cellZoneSizes) << nl
        << "    fZoneSizes  : " << flatOutput(faceZoneSizes) << nl
        << "    pZoneSizes  : " << flatOutput(pointZoneSizes) << nl
        << endl;
}


void Foam::polyTopoChange::getMergeSets
(
    const labelUList& reverseCellMap,
    const labelUList& cellMap,
    List<objectMap>& cellsFromCells
)
{
    // Per new cell the number of old cells that have been merged into it
    labelList nMerged(cellMap.size(), 1);

    forAll(reverseCellMap, oldCelli)
    {
        const label newCelli = reverseCellMap[oldCelli];

        if (newCelli < -1)
        {
            label mergeCelli = -newCelli-2;

            nMerged[mergeCelli]++;
        }
    }

    // From merged cell to set index
    labelList cellToMergeSet(cellMap.size(), -1);

    label nSets = 0;

    forAll(nMerged, celli)
    {
        if (nMerged[celli] > 1)
        {
            cellToMergeSet[celli] = nSets++;
        }
    }

    // Collect cell labels.
    // Each objectMap will have
    // - index : new mesh cell label
    // - masterObjects : list of old cells that have been merged. Element 0
    //                   will be the original destination cell label.

    cellsFromCells.setSize(nSets);

    forAll(reverseCellMap, oldCelli)
    {
        const label newCelli = reverseCellMap[oldCelli];

        if (newCelli < -1)
        {
            const label mergeCelli = -newCelli-2;

            // oldCelli was merged into mergeCelli

            const label setI = cellToMergeSet[mergeCelli];

            objectMap& mergeSet = cellsFromCells[setI];

            if (mergeSet.masterObjects().empty())
            {
                // First occurrence of master cell mergeCelli

                mergeSet.index() = mergeCelli;
                mergeSet.masterObjects().setSize(nMerged[mergeCelli]);

                // old master label
                mergeSet.masterObjects()[0] = cellMap[mergeCelli];

                // old slave label
                mergeSet.masterObjects()[1] = oldCelli;

                nMerged[mergeCelli] = 2;
            }
            else
            {
                mergeSet.masterObjects()[nMerged[mergeCelli]++] = oldCelli;
            }
        }
    }
}


bool Foam::polyTopoChange::hasValidPoints(const face& f) const
{
    for (const label fp : f)
    {
        if (fp < 0 || fp >= points_.size())
        {
            return false;
        }
    }
    return true;
}


Foam::pointField Foam::polyTopoChange::facePoints(const face& f) const
{
    pointField points(f.size());
    forAll(f, fp)
    {
        if (f[fp] < 0 && f[fp] >= points_.size())
        {
            FatalErrorInFunction
                << "Problem." << abort(FatalError);
        }
        points[fp] = points_[f[fp]];
    }
    return points;
}


void Foam::polyTopoChange::checkFace
(
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const label patchi,
    const label zoneI
) const
{
    if (nei == -1)
    {
        if (own == -1 && zoneI != -1)
        {
            // retired face
        }
        else if (patchi == -1 || patchi >= nPatches_)
        {
            FatalErrorInFunction
                << "Face has no neighbour (so external) but does not have"
                << " a valid patch" << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }
    else
    {
        if (patchi != -1)
        {
            FatalErrorInFunction
                << "Cannot both have valid patchi and neighbour" << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }

        if (nei <= own)
        {
            FatalErrorInFunction
                << "Owner cell label should be less than neighbour cell label"
                << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }

    if (f.size() < 3 || f.found(-1))
    {
        FatalErrorInFunction
            << "Illegal vertices in face"
            << nl
            << "f:" << f
            << " facei(-1 if added face):" << facei
            << " own:" << own << " nei:" << nei
            << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
    }
    if (facei >= 0 && facei < faces_.size() && faceRemoved(facei))
    {
        FatalErrorInFunction
            << "Face already marked for removal"
            << nl
            << "f:" << f
            << " facei(-1 if added face):" << facei
            << " own:" << own << " nei:" << nei
            << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
    }
    forAll(f, fp)
    {
        if (f[fp] < points_.size() && pointRemoved(f[fp]))
        {
            FatalErrorInFunction
                << "Face uses removed vertices"
                << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }
}


void Foam::polyTopoChange::makeCells
(
    const label nActiveFaces,
    labelList& cellFaces,
    labelList& cellFaceOffsets
) const
{
    cellFaces.setSize(2*nActiveFaces);
    cellFaceOffsets.setSize(cellMap_.size() + 1);

    // Faces per cell
    labelList nNbrs(cellMap_.size(), Zero);

    // 1. Count faces per cell

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceOwner_[facei] < 0)
        {
            pointField newPoints;
            if (facei < faces_.size())
            {
                const face& f = faces_[facei];
                newPoints.setSize(f.size(), vector::max);
                forAll(f, fp)
                {
                    if (f[fp] < points_.size())
                    {
                        newPoints[fp] = points_[f[fp]];
                    }
                }
            }


            FatalErrorInFunction
                << "Face " << facei << " is active but its owner has"
                << " been deleted. This is usually due to deleting cells"
                << " without modifying exposed faces to be boundary faces."
                << exit(FatalError);
        }
        nNbrs[faceOwner_[facei]]++;
    }
    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceNeighbour_[facei] >= 0)
        {
            nNbrs[faceNeighbour_[facei]]++;
        }
    }

    // 2. Calculate offsets

    cellFaceOffsets[0] = 0;
    forAll(nNbrs, celli)
    {
        cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + nNbrs[celli];
    }

    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label celli = faceOwner_[facei];

        cellFaces[cellFaceOffsets[celli] + nNbrs[celli]++] = facei;
    }

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label celli = faceNeighbour_[facei];

        if (celli >= 0)
        {
            cellFaces[cellFaceOffsets[celli] + nNbrs[celli]++] = facei;
        }
    }

    // Last offset points to beyond end of cellFaces.
    cellFaces.setSize(cellFaceOffsets[cellMap_.size()]);
}


// Create cell-cell addressing. Called after compaction (but before ordering)
// of faces
void Foam::polyTopoChange::makeCellCells
(
    const label nActiveFaces,
    CompactListList<label>& cellCells
) const
{
    // Neighbours per cell
    labelList nNbrs(cellMap_.size(), Zero);

    // 1. Count neighbours (through internal faces) per cell

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceNeighbour_[facei] >= 0)
        {
            nNbrs[faceOwner_[facei]]++;
            nNbrs[faceNeighbour_[facei]]++;
        }
    }

    // 2. Construct csr
    cellCells.setSize(nNbrs);


    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label nei = faceNeighbour_[facei];

        if (nei >= 0)
        {
            label own = faceOwner_[facei];
            cellCells(own, nNbrs[own]++) = nei;
            cellCells(nei, nNbrs[nei]++) = own;
        }
    }
}


// Cell ordering (based on bandCompression).
// Handles removed cells. Returns number of remaining cells.
Foam::label Foam::polyTopoChange::getCellOrder
(
    const CompactListList<label>& cellCellAddressing,
    labelList& oldToNew
) const
{
    const label nOldCells(cellCellAddressing.size());

    // Which cells are visited/unvisited
    bitSet unvisited(nOldCells, true);

    // Exclude removed cells
    for (label celli = 0; celli < nOldCells; ++celli)
    {
        if (cellRemoved(celli))
        {
            unvisited.unset(celli);
        }
    }

    // The new output order
    labelList newOrder(unvisited.count());


    // Various work arrays
    // ~~~~~~~~~~~~~~~~~~~

    //- Neighbour cells
    DynamicList<label> nbrCells;

    //- Neighbour ordering
    DynamicList<label> nbrOrder;

    //- Corresponding weights for neighbour cells
    DynamicList<label> weights;

    // FIFO buffer for renumbering.
    CircularBuffer<label> queuedCells(1024);


    label cellInOrder = 0;

    while (true)
    {
        // For a disconnected region find the lowest connected cell.
        label currCelli = -1;
        label minCount = labelMax;

        for (const label celli : unvisited)
        {
            const label nbrCount = cellCellAddressing[celli].size();

            if (minCount > nbrCount)
            {
                minCount = nbrCount;
                currCelli = celli;
            }
        }

        if (currCelli == -1)
        {
            break;
        }


        // Starting from currCelli - walk breadth-first

        queuedCells.push_back(currCelli);

        // Loop through queuedCells list. Add the first cell into the
        // cell order if it has not already been visited and ask for its
        // neighbours. If the neighbour in question has not been visited,
        // add it to the end of the queuedCells list

        while (!queuedCells.empty())
        {
            // Process as FIFO
            currCelli = queuedCells.front();
            queuedCells.pop_front();

            if (unvisited.test(currCelli))
            {
                // First visit...
                unvisited.unset(currCelli);

                // Add into cellOrder
                newOrder[cellInOrder] = currCelli;
                ++cellInOrder;

                // Find if the neighbours have been visited
                const auto& neighbours = cellCellAddressing[currCelli];

                // Add in increasing order of connectivity

                // 1. Count neighbours of unvisited neighbours
                nbrCells.clear();
                weights.clear();

                for (const label nbr : neighbours)
                {
                    const label nbrCount = cellCellAddressing[nbr].size();

                    if (unvisited.test(nbr))
                    {
                        // Not visited (or removed), add to the list
                        nbrCells.push_back(nbr);
                        weights.push_back(nbrCount);
                    }
                }

                // Resize DynamicList prior to sortedOrder
                nbrOrder.resize_nocopy(weights.size());

                // 2. Ascending order
                Foam::sortedOrder(weights, nbrOrder);

                // 3. Add to FIFO in sorted order
                for (const label nbrIdx : nbrOrder)
                {
                    queuedCells.push_back(nbrCells[nbrIdx]);
                }
            }
        }
    }

    // Now we have new-to-old in newOrder.
    newOrder.resize(cellInOrder);  // Extra safety, but should be a no-op

    // Invert to get old-to-new. Make sure removed (i.e. unmapped) cells are -1.
    oldToNew = invert(nOldCells, newOrder);

    return cellInOrder;
}


// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces and in patch order.
void Foam::polyTopoChange::getFaceOrder
(
    const label nActiveFaces,
    const labelUList& cellFaces,
    const labelUList& cellFaceOffsets,

    labelList& oldToNew,
    labelList& patchSizes,
    labelList& patchStarts
) const
{
    oldToNew.setSize(faceOwner_.size());
    oldToNew = -1;

    // First unassigned face
    label newFacei = 0;

    labelList nbr;
    labelList order;

    forAll(cellMap_, celli)
    {
        label startOfCell = cellFaceOffsets[celli];
        label nFaces = cellFaceOffsets[celli+1] - startOfCell;

        // Neighbouring cells
        nbr.setSize(nFaces);

        for (label i = 0; i < nFaces; i++)
        {
            label facei = cellFaces[startOfCell + i];

            label nbrCelli = faceNeighbour_[facei];

            if (facei >= nActiveFaces)
            {
                // Retired face.
                nbr[i] = -1;
            }
            else if (nbrCelli != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCelli == celli)
                {
                    nbrCelli = faceOwner_[facei];
                }

                if (celli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        sortedOrder(nbr, order);

        for (const label index : order)
        {
            if (nbr[index] != -1)
            {
                oldToNew[cellFaces[startOfCell + index]] = newFacei++;
            }
        }
    }


    // Pick up all patch faces in patch face order.
    patchStarts.setSize(nPatches_);
    patchStarts = 0;
    patchSizes.setSize(nPatches_);
    patchSizes = 0;

    if (nPatches_ > 0)
    {
        patchStarts[0] = newFacei;

        for (label facei = 0; facei < nActiveFaces; facei++)
        {
            if (region_[facei] >= 0)
            {
                patchSizes[region_[facei]]++;
            }
        }

        label facei = patchStarts[0];

        forAll(patchStarts, patchi)
        {
            patchStarts[patchi] = facei;
            facei += patchSizes[patchi];
        }
    }

    //if (debug)
    //{
    //    Pout<< "patchSizes:" << patchSizes << nl
    //        << "patchStarts:" << patchStarts << endl;
    //}

    labelList workPatchStarts(patchStarts);

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (region_[facei] >= 0)
        {
            oldToNew[facei] = workPatchStarts[region_[facei]]++;
        }
    }

    // Retired faces.
    for (label facei = nActiveFaces; facei < oldToNew.size(); facei++)
    {
        oldToNew[facei] = facei;
    }

    // Check done all faces.
    forAll(oldToNew, facei)
    {
        if (oldToNew[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position"
                << " for face " << facei
                << " owner " << faceOwner_[facei]
                << " neighbour " << faceNeighbour_[facei]
                << " region " << region_[facei] << endl
                << "This is usually caused by not specifying a patch for"
                << " a boundary face." << nl
                << "Switch on the polyTopoChange::debug flag to catch"
                << " this error earlier." << nl;
            if (hasValidPoints(faces_[facei]))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") " << facePoints(faces_[facei]);
            }
            FatalError << abort(FatalError);
        }
    }
}


// Reorder and compact faces according to map.
void Foam::polyTopoChange::reorderCompactFaces
(
    const label newSize,
    const labelUList& oldToNew
)
{
    reorder(oldToNew, faces_);
    faces_.setCapacity(newSize);

    reorder(oldToNew, region_);
    region_.setCapacity(newSize);

    reorder(oldToNew, faceOwner_);
    faceOwner_.setCapacity(newSize);

    reorder(oldToNew, faceNeighbour_);
    faceNeighbour_.setCapacity(newSize);

    // Update faceMaps.
    reorder(oldToNew, faceMap_);
    faceMap_.setCapacity(newSize);

    renumberReverseMap(oldToNew, reverseFaceMap_);

    renumberKey(oldToNew, faceFromPoint_);
    renumberKey(oldToNew, faceFromEdge_);
    inplaceReorder(oldToNew, flipFaceFlux_);
    flipFaceFlux_.setCapacity(newSize);
    renumberKey(oldToNew, faceZone_);
    inplaceReorder(oldToNew, faceZoneFlip_);
    faceZoneFlip_.setCapacity(newSize);
    if (faceAdditionalZones_.size())
    {
        // Extend to number of faces so oldToNew can be used
        faceAdditionalZones_.setSize(newSize);
        inplaceReorder(oldToNew, faceAdditionalZones_);
        //faceAdditionalZones_.setCapacity(newSize);
    }
}


void Foam::polyTopoChange::shrink()
{
    points_.shrink();
    pointMap_.shrink();
    reversePointMap_.shrink();
    pointAdditionalZones_.shrink();

    faces_.shrink();
    region_.shrink();
    faceOwner_.shrink();
    faceNeighbour_.shrink();
    faceMap_.shrink();
    reverseFaceMap_.shrink();
    faceAdditionalZones_.shrink();

    cellMap_.shrink();
    reverseCellMap_.shrink();
    cellZone_.shrink();
    cellAdditionalZones_.shrink();
}


// Compact all and orders points and faces:
// - points into internal followed by external points
// - internalfaces upper-triangular
// - externalfaces after internal ones.
void Foam::polyTopoChange::compact
(
    const bool orderCells,
    const bool orderPoints,
    label& nInternalPoints,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    shrink();

    // Compact points
    label nActivePoints = 0;
    {
        labelList localPointMap(points_.size(), -1);
        label newPointi = 0;

        if (!orderPoints)
        {
            nInternalPoints = -1;

            forAll(points_, pointi)
            {
                if (!pointRemoved(pointi) && !retiredPoints_.found(pointi))
                {
                    localPointMap[pointi] = newPointi++;
                }
            }
            nActivePoints = newPointi;
        }
        else
        {
            forAll(points_, pointi)
            {
                if (!pointRemoved(pointi) && !retiredPoints_.found(pointi))
                {
                    nActivePoints++;
                }
            }

            // Mark boundary points
            forAll(faceOwner_, facei)
            {
                if
                (
                   !faceRemoved(facei)
                 && faceOwner_[facei] >= 0
                 && faceNeighbour_[facei] < 0
                )
                {
                    // Valid boundary face
                    const face& f = faces_[facei];

                    forAll(f, fp)
                    {
                        label pointi = f[fp];

                        if (localPointMap[pointi] == -1)
                        {
                            if
                            (
                                pointRemoved(pointi)
                             || retiredPoints_.found(pointi)
                            )
                            {
                                FatalErrorInFunction
                                    << "Removed or retired point " << pointi
                                    << " in face " << f
                                    << " at position " << facei << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointi] = newPointi++;
                        }
                    }
                }
            }

            label nBoundaryPoints = newPointi;
            nInternalPoints = nActivePoints - nBoundaryPoints;

            // Move the boundary addressing up
            forAll(localPointMap, pointi)
            {
                if (localPointMap[pointi] != -1)
                {
                    localPointMap[pointi] += nInternalPoints;
                }
            }

            newPointi = 0;

            // Mark internal points
            forAll(faceOwner_, facei)
            {
                if
                (
                   !faceRemoved(facei)
                 && faceOwner_[facei] >= 0
                 && faceNeighbour_[facei] >= 0
                )
                {
                    // Valid internal face
                    const face& f = faces_[facei];

                    forAll(f, fp)
                    {
                        label pointi = f[fp];

                        if (localPointMap[pointi] == -1)
                        {
                            if
                            (
                                pointRemoved(pointi)
                             || retiredPoints_.found(pointi)
                            )
                            {
                                FatalErrorInFunction
                                    << "Removed or retired point " << pointi
                                    << " in face " << f
                                    << " at position " << facei << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointi] = newPointi++;
                        }
                    }
                }
            }

            if (newPointi != nInternalPoints)
            {
                FatalErrorInFunction
                    << "Problem." << abort(FatalError);
            }
            newPointi = nActivePoints;
        }

        for (const label pointi : retiredPoints_)
        {
            localPointMap[pointi] = newPointi++;
        }


        if (debug)
        {
            Pout<< "Points : active:" << nActivePoints
                << "  removed:" << points_.size()-newPointi << endl;
        }

        reorder(localPointMap, points_);
        points_.setCapacity(newPointi);

        // Update pointMaps
        reorder(localPointMap, pointMap_);
        pointMap_.setCapacity(newPointi);
        renumberReverseMap(localPointMap, reversePointMap_);

        renumberKey(localPointMap, pointZone_);
        renumber(localPointMap, retiredPoints_);
        if (pointAdditionalZones_.size())
        {
            reorder(localPointMap, pointAdditionalZones_);
        }

        // Use map to relabel face vertices
        forAll(faces_, facei)
        {
            face& f = faces_[facei];

            //labelList oldF(f);
            renumberCompact(localPointMap, f);

            if (!faceRemoved(facei) && f.size() < 3)
            {
                FatalErrorInFunction
                    << "Created illegal face " << f
                    //<< " from face " << oldF
                    << " at position:" << facei
                    << " when filtering removed points"
                    << abort(FatalError);
            }
        }
    }


    // Compact faces.
    {
        labelList localFaceMap(faces_.size(), -1);
        label newFacei = 0;

        forAll(faces_, facei)
        {
            if (!faceRemoved(facei) && faceOwner_[facei] >= 0)
            {
                localFaceMap[facei] = newFacei++;
            }
        }
        nActiveFaces_ = newFacei;

        forAll(faces_, facei)
        {
            if (!faceRemoved(facei) && faceOwner_[facei] < 0)
            {
                // Retired face
                localFaceMap[facei] = newFacei++;
            }
        }

        if (debug)
        {
            Pout<< "Faces : active:" << nActiveFaces_
                << "  removed:" << faces_.size()-newFacei << endl;
        }

        // Reorder faces.
        reorderCompactFaces(newFacei, localFaceMap);
    }

    // Compact cells.
    {
        labelList localCellMap;
        label newCelli;

        if (orderCells)
        {
            // Construct cellCell addressing
            CompactListList<label> cellCells;
            makeCellCells(nActiveFaces_, cellCells);

            // Cell ordering (based on bandCompression). Handles removed cells.
            newCelli = getCellOrder(cellCells, localCellMap);
        }
        else
        {
            // Compact out removed cells
            localCellMap.setSize(cellMap_.size());
            localCellMap = -1;

            newCelli = 0;
            forAll(cellMap_, celli)
            {
                if (!cellRemoved(celli))
                {
                    localCellMap[celli] = newCelli++;
                }
            }
        }

        if (debug)
        {
            Pout<< "Cells : active:" << newCelli
                << "  removed:" << cellMap_.size()-newCelli << endl;
        }

        // Renumber -if cells reordered or -if cells removed
        if (orderCells || (newCelli != cellMap_.size()))
        {
            reorder(localCellMap, cellMap_);
            cellMap_.setCapacity(newCelli);
            renumberReverseMap(localCellMap, reverseCellMap_);

            reorder(localCellMap, cellZone_);
            cellZone_.setCapacity(newCelli);
            if (cellAdditionalZones_.size())
            {
                reorder(localCellMap, cellAdditionalZones_);
                cellAdditionalZones_.setCapacity(newCelli);
            }

            renumberKey(localCellMap, cellFromPoint_);
            renumberKey(localCellMap, cellFromEdge_);
            renumberKey(localCellMap, cellFromFace_);

            // Renumber owner/neighbour. Take into account if neighbour suddenly
            // gets lower cell than owner.
            forAll(faceOwner_, facei)
            {
                label own = faceOwner_[facei];
                label nei = faceNeighbour_[facei];

                if (own >= 0)
                {
                    // Update owner
                    faceOwner_[facei] = localCellMap[own];

                    if (nei >= 0)
                    {
                        // Update neighbour.
                        faceNeighbour_[facei] = localCellMap[nei];

                        // Check if face needs reversing.
                        if
                        (
                            faceNeighbour_[facei] >= 0
                         && faceNeighbour_[facei] < faceOwner_[facei]
                        )
                        {
                            faces_[facei].flip();
                            std::swap(faceOwner_[facei], faceNeighbour_[facei]);
                            flipFaceFlux_.flip(facei);
                            faceZoneFlip_.flip(facei);
                            if (facei < faceAdditionalZones_.size())
                            {
                                for (auto& zas : faceAdditionalZones_[facei])
                                {
                                    // Flip sign
                                    zas = -zas;
                                }
                            }
                        }
                    }
                }
                else if (nei >= 0)
                {
                    // Update neighbour.
                    faceNeighbour_[facei] = localCellMap[nei];
                }
            }
        }
    }

    // Reorder faces into upper-triangular and patch ordering
    {
        // Create cells (packed storage)
        labelList cellFaces;
        labelList cellFaceOffsets;
        makeCells(nActiveFaces_, cellFaces, cellFaceOffsets);

        // Do upper triangular order and patch sorting
        labelList localFaceMap;
        getFaceOrder
        (
            nActiveFaces_,
            cellFaces,
            cellFaceOffsets,

            localFaceMap,
            patchSizes,
            patchStarts
        );

        // Reorder faces.
        reorderCompactFaces(localFaceMap.size(), localFaceMap);
    }
}


// Find faces to interpolate to create value for new face. Only used if
// face was inflated from edge or point. Internal faces should only be
// created from internal faces, external faces only from external faces
// (and ideally the same patch)
// Is bit problematic if there are no faces to select, i.e. in polyDualMesh
// an internal face can be created from a boundary edge with no internal
// faces connected to it.
Foam::labelList Foam::polyTopoChange::selectFaces
(
    const primitiveMesh& mesh,
    const labelUList& faceLabels,
    const bool internalFacesOnly
)
{
    label nFaces = 0;

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        if (internalFacesOnly == mesh.isInternalFace(facei))
        {
            nFaces++;
        }
    }

    labelList collectedFaces;

    if (nFaces == 0)
    {
        // Did not find any faces of the correct type so just use any old
        // face.
        collectedFaces = faceLabels;
    }
    else
    {
        collectedFaces.setSize(nFaces);

        nFaces = 0;

        forAll(faceLabels, i)
        {
            label facei = faceLabels[i];

            if (internalFacesOnly == mesh.isInternalFace(facei))
            {
                collectedFaces[nFaces++] = facei;
            }
        }
    }

    return collectedFaces;
}


// Calculate pointMap per patch (so from patch point label to old patch point
// label)
void Foam::polyTopoChange::calcPatchPointMap
(
    const UList<Map<label>>& oldPatchMeshPointMaps,
    const labelUList& patchMap,
    const polyBoundaryMesh& boundary,
    labelListList& patchPointMap
) const
{
    patchPointMap.setSize(patchMap.size());

    forAll(patchMap, patchi)
    {
        const label oldPatchi = patchMap[patchi];

        if (oldPatchi != -1)
        {
            const labelList& meshPoints = boundary[patchi].meshPoints();

            const Map<label>& oldMeshPointMap =
                oldPatchMeshPointMaps[oldPatchi];

            labelList& curPatchPointRnb = patchPointMap[patchi];

            curPatchPointRnb.setSize(meshPoints.size());

            forAll(meshPoints, i)
            {
                if (meshPoints[i] < pointMap_.size())
                {
                    // Check if old point was part of same patch
                    curPatchPointRnb[i] = oldMeshPointMap.lookup
                    (
                        pointMap_[meshPoints[i]],
                        -1
                    );
                }
                else
                {
                    curPatchPointRnb[i] = -1;
                }
            }
        }
    }
}


void Foam::polyTopoChange::calcFaceInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges,
    List<objectMap>& facesFromFaces
) const
{
    // Faces inflated from points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromPoints.setSize(faceFromPoint_.size());

    if (faceFromPoint_.size())
    {
        label nFacesFromPoints = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIters(faceFromPoint_, iter)
        {
            const label facei = iter.key();
            const label pointi = iter.val();

            // Get internal or patch faces using point on old mesh
            const bool internal = (region_[facei] == -1);

            facesFromPoints[nFacesFromPoints++] = objectMap
            (
                facei,
                selectFaces
                (
                    mesh,
                    mesh.pointFaces()[pointi],
                    internal
                )
            );
        }
    }


    // Faces inflated from edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromEdges.setSize(faceFromEdge_.size());

    if (faceFromEdge_.size())
    {
        label nFacesFromEdges = 0;

        // Collect all still existing faces connected to this edge.
        forAllConstIters(faceFromEdge_, iter)
        {
            const label facei = iter.key();
            const label edgei = iter.val();

            // Get internal or patch faces using edge on old mesh
            const bool internal = (region_[facei] == -1);

            facesFromEdges[nFacesFromEdges++] = objectMap
            (
                facei,
                selectFaces
                (
                    mesh,
                    mesh.edgeFaces(edgei),
                    internal
                )
            );
        }
    }


    // Faces from face merging
    // ~~~~~~~~~~~~~~~~~~~~~~~

    getMergeSets
    (
        reverseFaceMap_,
        faceMap_,
        facesFromFaces
    );
}


void Foam::polyTopoChange::calcCellInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& cellsFromPoints,
    List<objectMap>& cellsFromEdges,
    List<objectMap>& cellsFromFaces,
    List<objectMap>& cellsFromCells
) const
{
    cellsFromPoints.setSize(cellFromPoint_.size());

    if (cellFromPoint_.size())
    {
        label nCellsFromPoints = 0;

        // Collect all still existing cells connected to this point.
        forAllConstIters(cellFromPoint_, iter)
        {
            const label celli = iter.key();
            const label pointi = iter.val();

            cellsFromPoints[nCellsFromPoints++] = objectMap
            (
                celli,
                mesh.pointCells()[pointi]
            );
        }
    }


    cellsFromEdges.setSize(cellFromEdge_.size());

    if (cellFromEdge_.size())
    {
        label nCellsFromEdges = 0;

        // Collect all still existing cells connected to this edge.
        forAllConstIters(cellFromEdge_, iter)
        {
            const label celli = iter.key();
            const label edgei = iter.val();

            cellsFromEdges[nCellsFromEdges++] = objectMap
            (
                celli,
                mesh.edgeCells()[edgei]
            );
        }
    }


    cellsFromFaces.setSize(cellFromFace_.size());

    if (cellFromFace_.size())
    {
        label nCellsFromFaces = 0;

        labelList twoCells(2);

        // Collect all still existing faces connected to this point.
        forAllConstIters(cellFromFace_, iter)
        {
            const label celli = iter.key();
            const label oldFacei = iter.val();

            if (mesh.isInternalFace(oldFacei))
            {
                twoCells[0] = mesh.faceOwner()[oldFacei];
                twoCells[1] = mesh.faceNeighbour()[oldFacei];
                cellsFromFaces[nCellsFromFaces++] = objectMap
                (
                    celli,
                    twoCells
                );
            }
            else
            {
                cellsFromFaces[nCellsFromFaces++] = objectMap
                (
                    celli,
                    labelList(1, mesh.faceOwner()[oldFacei])
                );
            }
        }
    }


    // Cells from cell merging
    // ~~~~~~~~~~~~~~~~~~~~~~~

    getMergeSets
    (
        reverseCellMap_,
        cellMap_,
        cellsFromCells
    );
}


void Foam::polyTopoChange::resetZones
(
    const polyMesh& mesh,
    polyMesh& newMesh,
    labelListList& pointZoneMap,
    labelListList& faceZoneFaceMap,
    labelListList& cellZoneMap
) const
{
    // pointZones
    // ~~~~~~~~~~

    pointZoneMap.setSize(mesh.pointZones().size());
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Count points per zone

        labelList nPoints(pointZones.size(), Zero);

        forAllConstIters(pointZone_, iter)
        {
            const label pointi = iter.key();
            const label zonei = iter.val();

            if (zonei < 0 || zonei >= pointZones.size())
            {
                FatalErrorInFunction
                    << "Illegal zoneID " << zonei << " for point "
                    << pointi << " coord " << mesh.points()[pointi]
                    << abort(FatalError);
            }
            nPoints[zonei]++;
        }
        forAll(pointAdditionalZones_, pointi)
        {
            for (const label zonei : pointAdditionalZones_[pointi])
            {
                if (zonei < 0 || zonei >= pointZones.size())
                {
                    FatalErrorInFunction
                        << "Illegal zoneID " << zonei << " for point "
                        << pointi << abort(FatalError);
                }
                nPoints[zonei]++;
            }
        }

        // Distribute points per zone

        labelListList addressing(pointZones.size());
        forAll(addressing, zonei)
        {
            addressing[zonei].setSize(nPoints[zonei]);
        }
        nPoints = 0;

        forAllConstIters(pointZone_, iter)
        {
            const label pointi = iter.key();
            const label zonei = iter.val();

            addressing[zonei][nPoints[zonei]++] = pointi;
        }
        forAll(pointAdditionalZones_, pointi)
        {
            for (const label zonei : pointAdditionalZones_[pointi])
            {
                addressing[zonei][nPoints[zonei]++] = pointi;
            }
        }
        // Sort the addressing
        forAll(addressing, zonei)
        {
            stableSort(addressing[zonei]);
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get pointZoneMap.
        forAll(addressing, zonei)
        {
            const pointZone& oldZone = pointZones[zonei];
            const labelList& newZoneAddr = addressing[zonei];

            labelList& curPzRnb = pointZoneMap[zonei];
            curPzRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < pointMap_.size())
                {
                    curPzRnb[i] = oldZone.whichPoint(pointMap_[newZoneAddr[i]]);
                }
                else
                {
                    curPzRnb[i] = -1;
                }
            }
        }

        // Reset the addressing on the zone
        newMesh.pointZones().clearAddressing();
        forAll(newMesh.pointZones(), zonei)
        {
            if (debug)
            {
                Pout<< "pointZone:" << zonei
                    << "  name:" << newMesh.pointZones()[zonei].name()
                    << "  size:" << addressing[zonei].size()
                    << endl;
            }

            newMesh.pointZones()[zonei] = addressing[zonei];
        }
    }


    // faceZones
    // ~~~~~~~~~

    faceZoneFaceMap.setSize(mesh.faceZones().size());
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList nFaces(faceZones.size(), Zero);

        forAllConstIters(faceZone_, iter)
        {
            const label facei = iter.key();
            const label zonei = iter.val();

            if (zonei < 0 || zonei >= faceZones.size())
            {
                FatalErrorInFunction
                    << "Illegal zoneID " << zonei << " for face "
                    << facei
                    << abort(FatalError);
            }
            nFaces[zonei]++;
        }
        forAll(faceAdditionalZones_, facei)
        {
            for (const label zoneAndSign : faceAdditionalZones_[facei])
            {
                const label zonei = mag(zoneAndSign)-1;
                if (zonei < 0 || zonei >= faceZones.size())
                {
                    FatalErrorInFunction
                        << "Illegal zoneID " << zonei << " for face "
                        << facei << abort(FatalError);
                }
                nFaces[zonei]++;
            }
        }

        labelListList addressing(faceZones.size());
        boolListList flipMode(faceZones.size());

        forAll(addressing, zonei)
        {
            addressing[zonei].setSize(nFaces[zonei]);
            flipMode[zonei].setSize(nFaces[zonei]);
        }
        nFaces = 0;

        forAllConstIters(faceZone_, iter)
        {
            const label facei = iter.key();
            const label zonei = iter.val();

            const label index = nFaces[zonei]++;
            addressing[zonei][index] = facei;
            flipMode[zonei][index] = faceZoneFlip_.test(facei);

            if (facei < faceAdditionalZones_.size())
            {
                for (const label zoneAndSign : faceAdditionalZones_[facei])
                {
                    const label zonei = mag(zoneAndSign)-1;
                    const bool flip = (zoneAndSign > 0);

                    const label index = nFaces[zonei]++;
                    addressing[zonei][index] = facei;
                    flipMode[zonei][index] = flip;
                }
            }
        }

        // Sort the addressing
        forAll(addressing, zonei)
        {
            labelList newToOld(sortedOrder(addressing[zonei]));
            {
                labelList newAddressing(addressing[zonei].size());
                forAll(newAddressing, i)
                {
                    newAddressing[i] = addressing[zonei][newToOld[i]];
                }
                addressing[zonei].transfer(newAddressing);
            }
            {
                boolList newFlipMode(flipMode[zonei].size());
                forAll(newFlipMode, i)
                {
                    newFlipMode[i] = flipMode[zonei][newToOld[i]];
                }
                flipMode[zonei].transfer(newFlipMode);
            }
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get faceZoneFaceMap.
        forAll(addressing, zonei)
        {
            const faceZone& oldZone = faceZones[zonei];
            const labelList& newZoneAddr = addressing[zonei];

            labelList& curFzFaceRnb = faceZoneFaceMap[zonei];

            curFzFaceRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < faceMap_.size())
                {
                    curFzFaceRnb[i] =
                        oldZone.whichFace(faceMap_[newZoneAddr[i]]);
                }
                else
                {
                    curFzFaceRnb[i] = -1;
                }
            }
        }


        // Reset the addressing on the zone
        newMesh.faceZones().clearAddressing();
        forAll(newMesh.faceZones(), zonei)
        {
            if (debug)
            {
                Pout<< "faceZone:" << zonei
                    << "  name:" << newMesh.faceZones()[zonei].name()
                    << "  size:" << addressing[zonei].size()
                    << endl;
            }

            newMesh.faceZones()[zonei].resetAddressing
            (
                addressing[zonei],
                flipMode[zonei]
            );
        }
    }


    // cellZones
    // ~~~~~~~~~

    cellZoneMap.setSize(mesh.cellZones().size());
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList nCells(cellZones.size(), Zero);

        forAll(cellZone_, celli)
        {
            const label zonei = cellZone_[celli];

            if (zonei >= cellZones.size())
            {
                FatalErrorInFunction
                    << "Illegal zoneID " << zonei << " for cell "
                    << celli << abort(FatalError);
            }

            if (zonei >= 0)
            {
                nCells[zonei]++;
            }
        }
        forAll(cellAdditionalZones_, celli)
        {
            for (const label zonei : cellAdditionalZones_[celli])
            {
                if (zonei < 0 || zonei >= cellZones.size())
                {
                    FatalErrorInFunction
                        << "Illegal zoneID " << zonei << " for cell "
                        << celli << abort(FatalError);
                }
                nCells[zonei]++;
            }
        }


        labelListList addressing(cellZones.size());
        forAll(addressing, zonei)
        {
            addressing[zonei].setSize(nCells[zonei]);
        }
        nCells = 0;

        forAll(cellZone_, celli)
        {
            const label zonei = cellZone_[celli];

            if (zonei >= 0)
            {
                addressing[zonei][nCells[zonei]++] = celli;
            }
        }
        forAll(cellAdditionalZones_, celli)
        {
            for (const label zonei : cellAdditionalZones_[celli])
            {
                addressing[zonei][nCells[zonei]++] = celli;
            }
        }
        // Sort the addressing
        forAll(addressing, zonei)
        {
            stableSort(addressing[zonei]);
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get cellZoneMap.
        forAll(addressing, zonei)
        {
            const cellZone& oldZone = cellZones[zonei];
            const labelList& newZoneAddr = addressing[zonei];

            labelList& curCellRnb = cellZoneMap[zonei];

            curCellRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < cellMap_.size())
                {
                    curCellRnb[i] =
                        oldZone.whichCell(cellMap_[newZoneAddr[i]]);
                }
                else
                {
                    curCellRnb[i] = -1;
                }
            }
        }

        // Reset the addressing on the zone
        newMesh.cellZones().clearAddressing();
        forAll(newMesh.cellZones(), zonei)
        {
            if (debug)
            {
                Pout<< "cellZone:" << zonei
                    << "  name:" << newMesh.cellZones()[zonei].name()
                    << "  size:" << addressing[zonei].size()
                    << endl;
            }

            newMesh.cellZones()[zonei] = addressing[zonei];
        }
    }
}


void Foam::polyTopoChange::calcFaceZonePointMap
(
    const polyMesh& mesh,
    const UList<Map<label>>& oldFaceZoneMeshPointMaps,
    labelListList& faceZonePointMap
) const
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    faceZonePointMap.setSize(faceZones.size());

    forAll(faceZones, zonei)
    {
        const faceZone& newZone = faceZones[zonei];

        const labelList& newZoneMeshPoints = newZone.patch().meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldFaceZoneMeshPointMaps[zonei];

        labelList& curFzPointRnb = faceZonePointMap[zonei];

        curFzPointRnb.setSize(newZoneMeshPoints.size());

        forAll(newZoneMeshPoints, pointi)
        {
            if (newZoneMeshPoints[pointi] < pointMap_.size())
            {
                curFzPointRnb[pointi] =
                    oldZoneMeshPointMap.lookup
                    (
                        pointMap_[newZoneMeshPoints[pointi]],
                        -1
                    );
            }
            else
            {
                curFzPointRnb[pointi] = -1;
            }
        }
    }
}


void Foam::polyTopoChange::reorderCoupledFaces
(
    const bool syncParallel,
    const polyBoundaryMesh& oldBoundary,
    const labelUList& patchMap,     // new to old patches
    const labelUList& patchStarts,
    const labelUList& patchSizes,
    const pointField& points
)
{
    // Mapping for faces (old to new). Extends over all mesh faces for
    // convenience (could be just the external faces)
    labelList oldToNew(identity(faces_.size()));

    // Rotation on new faces.
    labelList rotation(faces_.size(), Zero);

    // Allocate unique tag for all comms
    const int oldTag = UPstream::incrMsgType();

    PstreamBuffers pBufs;

    // Send ordering
    forAll(patchMap, patchi)
    {
        const label oldPatchi = patchMap[patchi];

        if
        (
            oldPatchi != -1
         && (syncParallel || !isA<processorPolyPatch>(oldBoundary[oldPatchi]))
        )
        {
            oldBoundary[oldPatchi].initOrder
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchi],
                        patchStarts[patchi]
                    ),
                    points
                )
            );
        }
    }

    if (syncParallel)
    {
        pBufs.finishedSends();
    }

    // Receive and calculate ordering

    bool anyChanged = false;

    forAll(patchMap, patchi)
    {
        const label oldPatchi = patchMap[patchi];

        if
        (
            oldPatchi != -1
         && (syncParallel || !isA<processorPolyPatch>(oldBoundary[oldPatchi]))
        )
        {
            labelList patchFaceMap(patchSizes[patchi], -1);
            labelList patchFaceRotation(patchSizes[patchi], Zero);

            const bool changed = oldBoundary[oldPatchi].order
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchi],
                        patchStarts[patchi]
                    ),
                    points
                ),
                patchFaceMap,
                patchFaceRotation
            );

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                const label start = patchStarts[patchi];

                forAll(patchFaceMap, patchFacei)
                {
                    oldToNew[patchFacei + start] =
                        start + patchFaceMap[patchFacei];
                }

                forAll(patchFaceRotation, patchFacei)
                {
                    rotation[patchFacei + start] =
                        patchFaceRotation[patchFacei];
                }

                anyChanged = true;
            }
        }
    }

    if (syncParallel ? returnReduceOr(anyChanged) : anyChanged)
    {
        // Reorder faces according to oldToNew.
        reorderCompactFaces(oldToNew.size(), oldToNew);

        // Rotate faces (rotation is already in new face indices).
        forAll(rotation, facei)
        {
            if (rotation[facei] != 0)
            {
                inplaceRotateList<List, label>(faces_[facei], rotation[facei]);
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


void Foam::polyTopoChange::compactAndReorder
(
    const polyMesh& mesh,
    const labelUList& patchMap,     // from new to old patch
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints,

    label& nInternalPoints,
    pointField& newPoints,
    labelList& patchSizes,
    labelList& patchStarts,
    List<objectMap>& pointsFromPoints,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges,
    List<objectMap>& facesFromFaces,
    List<objectMap>& cellsFromPoints,
    List<objectMap>& cellsFromEdges,
    List<objectMap>& cellsFromFaces,
    List<objectMap>& cellsFromCells,
    List<Map<label>>& oldPatchMeshPointMaps,
    labelList& oldPatchNMeshPoints,
    labelList& oldPatchStarts,
    List<Map<label>>& oldFaceZoneMeshPointMaps
)
{
    if (patchMap.size() != nPatches_)
    {
        FatalErrorInFunction
            << "polyTopoChange was constructed with a mesh with "
            << nPatches_ << " patches." << endl
            << "The mesh now provided has a different number of patches "
            << patchMap.size()
            << " which is illegal" << endl
            << abort(FatalError);
    }

    // Remove any holes from points/faces/cells and sort faces.
    // Sets nActiveFaces_.
    compact(orderCells, orderPoints, nInternalPoints, patchSizes, patchStarts);

    // Transfer points to pointField. points_ are now cleared!
    // Only done since e.g. reorderCoupledFaces requires pointField.
    newPoints.transfer(points_);

    // Reorder any coupled faces
    reorderCoupledFaces
    (
        syncParallel,
        mesh.boundaryMesh(),
        patchMap,
        patchStarts,
        patchSizes,
        newPoints
    );


    // Calculate inflation/merging maps
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are for the new face(/point/cell) the old faces whose value
    // needs to be
    // averaged/summed to get the new value. These old faces come either from
    // merged old faces (face remove into other face),
    // the old edgeFaces (inflate from edge) or the old pointFaces (inflate
    // from point). As an additional complexity will use only internal faces
    // to create new value for internal face and vice versa only patch
    // faces to to create patch face value.

    // For point only point merging
    getMergeSets
    (
        reversePointMap_,
        pointMap_,
        pointsFromPoints
    );

    calcFaceInflationMaps
    (
        mesh,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces
    );

    calcCellInflationMaps
    (
        mesh,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells
    );

    // Clear inflation info
    {
        faceFromPoint_.clearStorage();
        faceFromEdge_.clearStorage();

        cellFromPoint_.clearStorage();
        cellFromEdge_.clearStorage();
        cellFromFace_.clearStorage();
    }


    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    // Grab patch mesh point maps
    oldPatchMeshPointMaps.setSize(boundary.size());
    oldPatchNMeshPoints.setSize(boundary.size());
    oldPatchStarts.setSize(boundary.size());

    forAll(boundary, patchi)
    {
        // Copy old face zone mesh point maps
        oldPatchMeshPointMaps[patchi] = boundary[patchi].meshPointMap();
        oldPatchNMeshPoints[patchi] = boundary[patchi].meshPoints().size();
        oldPatchStarts[patchi] = boundary[patchi].start();
    }

    // Grab old face zone mesh point maps.
    // These need to be saved before resetting the mesh and are used
    // later on to calculate the faceZone pointMaps.
    oldFaceZoneMeshPointMaps.setSize(mesh.faceZones().size());

    forAll(mesh.faceZones(), zonei)
    {
        const faceZone& oldZone = mesh.faceZones()[zonei];

        oldFaceZoneMeshPointMaps[zonei] = oldZone.patch().meshPointMap();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyTopoChange::polyTopoChange(const label nPatches, const bool strict)
:
    strict_(strict),
    nPatches_(nPatches),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    pointZone_(0),
    retiredPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    faceFromPoint_(0),
    faceFromEdge_(0),
    flipFaceFlux_(0),
    faceZone_(0),
    faceZoneFlip_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0),
    cellFromPoint_(0),
    cellFromEdge_(0),
    cellFromFace_(0),
    cellZone_(0)
{}


// Construct from components
Foam::polyTopoChange::polyTopoChange
(
    const polyMesh& mesh,
    const bool strict
)
:
    strict_(strict),
    nPatches_(0),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    pointZone_(0),
    retiredPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    faceFromPoint_(0),
    faceFromEdge_(0),
    flipFaceFlux_(0),
    faceZone_(0),
    faceZoneFlip_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0),
    cellFromPoint_(0),
    cellFromEdge_(0),
    cellFromFace_(0),
    cellZone_(0)
{
    addMesh
    (
        mesh,
        identity(mesh.boundaryMesh().size()),
        identity(mesh.pointZones().size()),
        identity(mesh.faceZones().size()),
        identity(mesh.cellZones().size())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyTopoChange::clear()
{
    points_.clearStorage();
    pointMap_.clearStorage();
    reversePointMap_.clearStorage();
    pointZone_.clearStorage();
    pointAdditionalZones_.clearStorage();
    retiredPoints_.clearStorage();

    faces_.clearStorage();
    region_.clearStorage();
    faceOwner_.clearStorage();
    faceNeighbour_.clearStorage();
    faceMap_.clearStorage();
    reverseFaceMap_.clearStorage();
    faceFromPoint_.clearStorage();
    faceFromEdge_.clearStorage();
    flipFaceFlux_.clearStorage();
    faceZone_.clearStorage();
    faceZoneFlip_.clearStorage();
    faceAdditionalZones_.clearStorage();
    nActiveFaces_ = 0;

    cellMap_.clearStorage();
    reverseCellMap_.clearStorage();
    cellZone_.clearStorage();
    cellAdditionalZones_.clearStorage();
    cellFromPoint_.clearStorage();
    cellFromEdge_.clearStorage();
    cellFromFace_.clearStorage();
}


void Foam::polyTopoChange::addMesh
(
    const polyMesh& mesh,
    const labelUList& patchMap,
    const labelUList& pointZoneMap,
    const labelUList& faceZoneMap,
    const labelUList& cellZoneMap
)
{
    label maxRegion = nPatches_ - 1;
    for (const label regioni : patchMap)
    {
        maxRegion = max(maxRegion, regioni);
    }
    nPatches_ = maxRegion + 1;


    // Add points
    {
        const pointField& points = mesh.points();
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Extend
        points_.setCapacity(points_.size() + points.size());
        pointMap_.setCapacity(pointMap_.size() + points.size());
        reversePointMap_.setCapacity(reversePointMap_.size() + points.size());
        pointZone_.resize(pointZone_.size() + points.size()/128);

        // Add points without pointZone
        labelList newPointID(points.size());
        forAll(newPointID, pointi)
        {
            // Add point from point
            newPointID[pointi] = addPoint
            (
                points[pointi],
                pointi,
                 -1,
                true
            );
        }

        // Modify for pointZones
        DynamicList<label> zones;
        forAll(pointZones, zonei)
        {
            for (const label pointi : pointZones[zonei])
            {
                // Get previous zones
                this->pointZones(newPointID[pointi], zones);
                for (auto& zonei : zones)
                {
                    zonei = pointZoneMap[zonei];
                }
                // Add new zone
                zones.append(pointZoneMap[zonei]);
                modifyPoint
                (
                    newPointID[pointi],
                    points_[newPointID[pointi]],
                    zones,
                    true
                );
            }
        }
    }

    // Add cells
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        // Resize

        // Note: polyMesh does not allow retired cells anymore. So allCells
        // always equals nCells
        label nAllCells = mesh.nCells();

        cellMap_.setCapacity(cellMap_.size() + nAllCells);
        reverseCellMap_.setCapacity(reverseCellMap_.size() + nAllCells);
        cellFromPoint_.resize(cellFromPoint_.size() + nAllCells/128);
        cellFromEdge_.resize(cellFromEdge_.size() + nAllCells/128);
        cellFromFace_.resize(cellFromFace_.size() + nAllCells/128);
        cellZone_.setCapacity(cellZone_.size() + nAllCells);


        // Add cells without cellZone
        labelList newCellID(nAllCells);
        for (label celli = 0; celli < nAllCells; celli++)
        {
            // Add cell from cell
            newCellID[celli] = addCell(-1, -1, -1, celli, -1);
        }

        // Modify for cellZones
        DynamicList<label> zones;
        forAll(cellZones, zonei)
        {
            for (const label celli : cellZones[zonei])
            {
                // Get previous zones
                this->cellZones(newCellID[celli], zones);
                for (auto& zonei : zones)
                {
                    zonei = cellZoneMap[zonei];
                }
                // Add new zone
                zones.append(cellZoneMap[zonei]);
                modifyCell(newCellID[celli], zones);
            }
        }
    }

    // Add faces
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const faceList& faces = mesh.faces();
        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const faceZoneMesh& faceZones = mesh.faceZones();

        // Resize
        label nAllFaces = mesh.faces().size();

        faces_.setCapacity(faces_.size() + nAllFaces);
        region_.setCapacity(region_.size() + nAllFaces);
        faceOwner_.setCapacity(faceOwner_.size() + nAllFaces);
        faceNeighbour_.setCapacity(faceNeighbour_.size() + nAllFaces);
        faceMap_.setCapacity(faceMap_.size() + nAllFaces);
        reverseFaceMap_.setCapacity(reverseFaceMap_.size() + nAllFaces);
        faceFromPoint_.resize(faceFromPoint_.size() + nAllFaces/128);
        faceFromEdge_.resize(faceFromEdge_.size() + nAllFaces/128);
        flipFaceFlux_.setCapacity(faces_.size() + nAllFaces);
        faceZone_.resize(faceZone_.size() + nAllFaces/128);
        faceZoneFlip_.setCapacity(faceMap_.capacity());


        // Add faces (in mesh order) without faceZone
        labelList newFaceID(nAllFaces);

        // 1. Internal faces
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            newFaceID[facei] = addFace
            (
                faces[facei],
                faceOwner[facei],
                faceNeighbour[facei],
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                facei,                      // masterFaceID
                false,                      // flipFaceFlux
                -1,                         // patchID
                -1,                         // zoneID
                false                       // zoneFlip
            );
        }

        // 2. Patch faces
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.start() != faces_.size())
            {
                FatalErrorInFunction
                    << "Problem : "
                    << "Patch " << pp.name() << " starts at " << pp.start()
                    << endl
                    << "Current face counter at " << faces_.size() << endl
                    << "Are patches in incremental order?"
                    << abort(FatalError);
            }
            forAll(pp, patchFacei)
            {
                const label facei = pp.start() + patchFacei;

                newFaceID[facei] = addFace
                (
                    faces[facei],
                    faceOwner[facei],
                    -1,                         // neighbour
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    facei,                      // masterFaceID
                    false,                      // flipFaceFlux
                    patchMap[patchi],           // patchID
                    -1,                         // zoneID
                    false                       // zoneFlip
                );
            }
        }


        // Modify for faceZones

        DynamicList<label> zones;
        DynamicList<bool> flips;
        forAll(faceZones, zonei)
        {
            const auto& fz = faceZones[zonei];

            forAll(fz, i)
            {
                const label facei = fz[i];
                const bool flip = fz.flipMap()[i];

                // Modify face index for combined zones
                const label newFacei = newFaceID[facei];

                // Get previous zones
                this->faceZones(newFacei, zones, flips);
                for (auto& zonei : zones)
                {
                    zonei = faceZoneMap[zonei];
                }
                // Add new zone
                zones.append(faceZoneMap[zonei]);
                flips.append(flip);

                modifyFace
                (
                    faces_[newFacei],
                    newFacei,
                    faceOwner_[newFacei],
                    faceNeighbour_[newFacei],
                    flipFaceFlux_(newFacei),
                    region_[newFacei],
                    zones,
                    flips
                );
            }
        }
    }
}


void Foam::polyTopoChange::setCapacity
(
    const label nPoints,
    const label nFaces,
    const label nCells
)
{
    points_.setCapacity(nPoints);
    pointMap_.setCapacity(nPoints);
    reversePointMap_.setCapacity(nPoints);
    pointZone_.resize(pointZone_.size() + nPoints/128);

    faces_.setCapacity(nFaces);
    region_.setCapacity(nFaces);
    faceOwner_.setCapacity(nFaces);
    faceNeighbour_.setCapacity(nFaces);
    faceMap_.setCapacity(nFaces);
    reverseFaceMap_.setCapacity(nFaces);
    faceFromPoint_.resize(faceFromPoint_.size() + nFaces/128);
    faceFromEdge_.resize(faceFromEdge_.size() + nFaces/128);
    flipFaceFlux_.setCapacity(nFaces);
    faceZone_.resize(faceZone_.size() + nFaces/128);
    faceZoneFlip_.setCapacity(nFaces);

    cellMap_.setCapacity(nCells);
    reverseCellMap_.setCapacity(nCells);
    cellFromPoint_.resize(cellFromPoint_.size() + nCells/128);
    cellFromEdge_.resize(cellFromEdge_.size() + nCells/128);
    cellFromFace_.resize(cellFromFace_.size() + nCells/128);
    cellZone_.setCapacity(nCells);
}


Foam::label Foam::polyTopoChange::setAction(const topoAction& action)
{
    if (isType<polyAddPoint>(action))
    {
        const polyAddPoint& pap = refCast<const polyAddPoint>(action);

        return addPoint
        (
            pap.newPoint(),
            pap.masterPointID(),
            pap.zoneID(),
            pap.inCell()
        );
    }
    else if (isType<polyModifyPoint>(action))
    {
        const polyModifyPoint& pmp = refCast<const polyModifyPoint>(action);

        modifyPoint
        (
            pmp.pointID(),
            pmp.newPoint(),
            pmp.zoneID(),
            pmp.inCell()
        );

        return -1;
    }
    else if (isType<polyRemovePoint>(action))
    {
        const polyRemovePoint& prp = refCast<const polyRemovePoint>(action);

        removePoint(prp.pointID(), prp.mergePointID());

        return -1;
    }
    else if (isType<polyAddFace>(action))
    {
        const polyAddFace& paf = refCast<const polyAddFace>(action);

        return addFace
        (
            paf.newFace(),
            paf.owner(),
            paf.neighbour(),
            paf.masterPointID(),
            paf.masterEdgeID(),
            paf.masterFaceID(),
            paf.flipFaceFlux(),
            paf.patchID(),
            paf.zoneID(),
            paf.zoneFlip()
        );
    }
    else if (isType<polyModifyFace>(action))
    {
        const polyModifyFace& pmf = refCast<const polyModifyFace>(action);

        modifyFace
        (
            pmf.newFace(),
            pmf.faceID(),
            pmf.owner(),
            pmf.neighbour(),
            pmf.flipFaceFlux(),
            pmf.patchID(),
            pmf.zoneID(),
            pmf.zoneFlip()
        );

        return -1;
    }
    else if (isType<polyRemoveFace>(action))
    {
        const polyRemoveFace& prf = refCast<const polyRemoveFace>(action);

        removeFace(prf.faceID(), prf.mergeFaceID());

        return -1;
    }
    else if (isType<polyAddCell>(action))
    {
        const polyAddCell& pac = refCast<const polyAddCell>(action);

        return addCell
        (
            pac.masterPointID(),
            pac.masterEdgeID(),
            pac.masterFaceID(),
            pac.masterCellID(),
            pac.zoneID()
        );
    }
    else if (isType<polyModifyCell>(action))
    {
        const polyModifyCell& pmc = refCast<const polyModifyCell>(action);

        if (pmc.removeFromZone())
        {
            modifyCell(pmc.cellID(), -1);
        }
        else
        {
            modifyCell(pmc.cellID(), pmc.zoneID());
        }

        return -1;
    }
    else if (isType<polyRemoveCell>(action))
    {
        const polyRemoveCell& prc = refCast<const polyRemoveCell>(action);

        removeCell(prc.cellID(), prc.mergeCellID());

        return -1;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown type of topoChange: " << action.type()
            << abort(FatalError);

        // Dummy return to keep compiler happy
        return -1;
    }
}


Foam::label Foam::polyTopoChange::addPoint
(
    const point& pt,
    const label masterPointID,
    const label zoneID,
    const bool inCell
)
{
    const label pointi = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointi);

    if (zoneID >= 0)
    {
        pointZone_.insert(pointi, zoneID);
    }
    // Delay extending the pointAdditionalZones

    if (!inCell)
    {
        retiredPoints_.insert(pointi);
    }

    return pointi;
}


Foam::label Foam::polyTopoChange::addPoint
(
    const point& pt,
    const label masterPointID,
    const labelUList& zoneIDs,
    const bool inCell
)
{
    const label pointi = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointi);

    if (zoneIDs.size())
    {
        const label minIndex = findMin(zoneIDs);

        pointZone_.set(pointi, zoneIDs[minIndex]);

        // Clear any additional storage
        if (pointi < pointAdditionalZones_.size())
        {
            pointAdditionalZones_[pointi].clear();
        }
        // (auto-vivify and) store additional zones
        forAll(zoneIDs, i)
        {
            if (i != minIndex)
            {
                if (zoneIDs[i] == zoneIDs[minIndex])
                {
                    FatalErrorInFunction << "Duplicates in zones "
                        << flatOutput(zoneIDs) << " for point " << pointi
                        << exit(FatalError);
                }
                pointAdditionalZones_(pointi).push_uniq(zoneIDs[i]);
            }
        }
    }

    if (!inCell)
    {
        retiredPoints_.insert(pointi);
    }

    return pointi;
}


void Foam::polyTopoChange::modifyPoint
(
    const label pointi,
    const point& pt,
    const label zoneID,
    const bool inCell,
    const bool multiZone
)
{
    if (pointi < 0 || pointi >= points_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(pointi) || pointMap_[pointi] == -1)
    {
        FatalErrorInFunction
            << "point " << pointi << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointi] = pt;

    if (!multiZone)
    {
        if (zoneID >= 0)
        {
            pointZone_.set(pointi, zoneID);
        }
        else
        {
            pointZone_.erase(pointi);
        }
        if (pointi < pointAdditionalZones_.size())
        {
            pointAdditionalZones_[pointi].clear();
        }
    }
    else
    {
        auto fnd = pointZone_.find(pointi);
        if (!fnd)
        {
            pointZone_.set(pointi, zoneID);

            // Check that no leftover additional zones
            if (pointAdditionalZones_(pointi).size())
            {
                FatalErrorInFunction
                    << "Additional zones for point:"
                    << pointAdditionalZones_(pointi)
                    << abort(FatalError);
            }
        }
        else if (zoneID == -1)
        {
            // Clear
            pointZone_.erase(fnd);
            if (pointi < pointAdditionalZones_.size())
            {
                pointAdditionalZones_[pointi].clear();
            }
        }
        else if (zoneID != fnd())
        {
            // New zone. Store in additional storage.
            pointAdditionalZones_(pointi).push_uniq(zoneID);
        }
    }

    if (inCell)
    {
        retiredPoints_.erase(pointi);
    }
    else
    {
        retiredPoints_.insert(pointi);
    }
}


void Foam::polyTopoChange::modifyPoint
(
    const label pointi,
    const point& pt,
    const labelUList& zoneIDs,
    const bool inCell
)
{
    if (pointi < 0 || pointi >= points_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(pointi) || pointMap_[pointi] == -1)
    {
        FatalErrorInFunction
            << "point " << pointi << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointi] = pt;

    if (zoneIDs.size())
    {
        if (zoneIDs.found(-1))
        {
            FatalErrorInFunction
                << "zones cannot contain -1 : " << flatOutput(zoneIDs)
                << " for point:" << pointi
                << abort(FatalError);
        }

        pointZone_.set(pointi, zoneIDs[0]);
        if (pointi < pointAdditionalZones_.size())
        {
            pointAdditionalZones_[pointi].clear();
        }
        for (label i = 1; i < zoneIDs.size(); ++i)
        {
            pointAdditionalZones_(pointi).push_uniq(zoneIDs[i]);
        }
    }
    else
    {
        pointZone_.erase(pointi);
        if (pointi < pointAdditionalZones_.size())
        {
            pointAdditionalZones_[pointi].clear();
        }
    }

    if (inCell)
    {
        retiredPoints_.erase(pointi);
    }
    else
    {
        retiredPoints_.insert(pointi);
    }
}

void Foam::polyTopoChange::movePoints(const pointField& newPoints)
{
    if (newPoints.size() != points_.size())
    {
        FatalErrorInFunction
            << "illegal pointField size." << endl
            << "Size:" << newPoints.size() << endl
            << "Points in mesh:" << points_.size()
            << abort(FatalError);
    }

    forAll(points_, pointi)
    {
        points_[pointi] = newPoints[pointi];
    }
}


void Foam::polyTopoChange::removePoint
(
    const label pointi,
    const label mergePointi
)
{
    if (pointi < 0 || pointi >= points_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (pointRemoved(pointi) || pointMap_[pointi] == -1)
    )
    {
        FatalErrorInFunction
            << "point " << pointi << " already marked for removal" << nl
            << "Point:" << points_[pointi] << " pointMap:" << pointMap_[pointi]
            << abort(FatalError);
    }

    if (pointi == mergePointi)
    {
        FatalErrorInFunction
            << "Cannot remove/merge point " << pointi << " onto itself."
            << abort(FatalError);
    }

    points_[pointi] = point::max;
    pointMap_[pointi] = -1;
    if (mergePointi >= 0)
    {
        reversePointMap_[pointi] = -mergePointi-2;
    }
    else
    {
        reversePointMap_[pointi] = -1;
    }
    pointZone_.erase(pointi);
    if (pointi < pointAdditionalZones_.size())
    {
        pointAdditionalZones_[pointi].clear();
    }
    retiredPoints_.erase(pointi);
}


Foam::label Foam::polyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, -1, own, nei, patchID, zoneID);
    }

    label facei = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterPointID >= 0)
    {
        faceMap_.append(-1);
        faceFromPoint_.insert(facei, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        faceMap_.append(-1);
        faceFromEdge_.insert(facei, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow inflate-from-nothing?
        //FatalErrorInFunction
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(facei);

    flipFaceFlux_.set(facei, flipFaceFlux);

    if (zoneID >= 0)
    {
        faceZone_.insert(facei, zoneID);
    }
    faceZoneFlip_.set(facei, zoneFlip);
    // Delay extending the faceAdditionalZones

    return facei;
}


Foam::label Foam::polyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const bool flipFaceFlux,
    const label patchID,
    const labelUList& zoneIDs,
    const UList<bool>& zoneFlips
)
{
    // Check validity
    if (debug)
    {
        // Note: no check on zones
        checkFace(f, -1, own, nei, patchID, -1);
    }

    label facei = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterPointID >= 0)
    {
        faceMap_.append(-1);
        faceFromPoint_.insert(facei, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        faceMap_.append(-1);
        faceFromEdge_.insert(facei, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow inflate-from-nothing?
        //FatalErrorInFunction
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(facei);

    flipFaceFlux_.set(facei, flipFaceFlux);

    if (zoneIDs.size())
    {
        const label minIndex = findMin(zoneIDs);

        faceZone_.set(facei, zoneIDs[minIndex]);
        faceZoneFlip_.set(facei, zoneFlips[minIndex]);

        // Clear any additional storage
        if (facei < faceAdditionalZones_.size())
        {
            faceAdditionalZones_[facei].clear();
        }
        // (auto-vivify and) store additional zones
        forAll(zoneIDs, i)
        {
            if (i != minIndex)
            {
                if (zoneIDs[i] == zoneIDs[minIndex])
                {
                    FatalErrorInFunction << "Duplicates in zones "
                        << flatOutput(zoneIDs) << " for face " << facei
                        << exit(FatalError);
                }
                const label zoneAndSign
                (
                    zoneFlips[i]
                  ? zoneIDs[i]+1
                  : -(zoneIDs[i]+1)
                );
                faceAdditionalZones_(facei).push_uniq(zoneAndSign);
            }
        }
    }

    return facei;
}


void Foam::polyTopoChange::modifyFace
(
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip,
    const bool multiZone
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, facei, own, nei, patchID, zoneID);
    }

    faces_[facei] = f;
    faceOwner_[facei] = own;
    faceNeighbour_[facei] = nei;
    region_[facei] = patchID;

    flipFaceFlux_.set(facei, flipFaceFlux);

    // Note: same logic as modifyCell except storage is now sparse instead
    //       of marked with -1

    if (!multiZone)
    {
        if (zoneID == -1)
        {
            faceZone_.erase(facei);
            faceZoneFlip_.set(facei, zoneFlip);
        }
        else
        {
            faceZone_.set(facei, zoneID);
            faceZoneFlip_.set(facei, zoneFlip);
        }
        if (facei < faceAdditionalZones_.size())
        {
            faceAdditionalZones_[facei].clear();
        }
    }
    else
    {
        auto fnd = faceZone_.find(facei);
        if (!fnd)
        {
            faceZone_.set(facei, zoneID);
            faceZoneFlip_.set(facei, zoneFlip);

            // Check that no leftover additional zones
            if (faceAdditionalZones_(facei).size())
            {
                FatalErrorInFunction
                    << "Additional zones for face:"
                    << faceAdditionalZones_(facei)
                    << abort(FatalError);
            }
        }
        else if (zoneID == -1)
        {
            // Clear
            faceZone_.erase(fnd);
            faceZoneFlip_.set(facei, false);
            if (facei < faceAdditionalZones_.size())
            {
                faceAdditionalZones_[facei].clear();
            }
        }
        else if (zoneID != fnd())
        {
            // New zone. Store in additional storage.
            const label zoneAndSign
            (
                zoneFlip
              ? zoneID+1
              : -(zoneID+1)
            );
            faceAdditionalZones_(facei).push_uniq(zoneAndSign);
        }
    }
}


void Foam::polyTopoChange::modifyFace
(
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label patchID,
    const labelUList& zoneIDs,
    const UList<bool>& zoneFlips
)
{
    // Check validity
    if (debug)
    {
        // Note: no check on zones
        checkFace(f, facei, own, nei, patchID, -1);
    }

    faces_[facei] = f;
    faceOwner_[facei] = own;
    faceNeighbour_[facei] = nei;
    region_[facei] = patchID;

    flipFaceFlux_.set(facei, flipFaceFlux);

    // Note: same logic as modifyCell except storage is now sparse instead
    //       of marked with -1

    if (zoneIDs.size())
    {
        if (zoneIDs.found(-1))
        {
            FatalErrorInFunction
                << "zones cannot contain -1 : " << flatOutput(zoneIDs)
                << " for face:" << facei
                << abort(FatalError);
        }

        faceZone_.set(facei, zoneIDs[0]);
        faceZoneFlip_.set(facei, zoneFlips[0]);
        if (facei < faceAdditionalZones_.size())
        {
            faceAdditionalZones_[facei].clear();
        }
        for (label i = 1; i < zoneIDs.size(); ++i)
        {
            const label zoneAndSign
            (
                zoneFlips[i]
              ? zoneIDs[i]+1
              : -(zoneIDs[i]+1)
            );
            faceAdditionalZones_(facei).push_uniq(zoneAndSign);
        }
    }
    else
    {
        faceZone_.erase(facei);
        faceZoneFlip_.set(facei, false);
        if (facei < faceAdditionalZones_.size())
        {
            faceAdditionalZones_[facei].clear();
        }
    }
}


void Foam::polyTopoChange::removeFace
(
    const label facei,
    const label mergeFacei
)
{
    if (facei < 0 || facei >= faces_.size())
    {
        FatalErrorInFunction
            << "illegal face label " << facei << endl
            << "Valid face labels are 0 .. " << faces_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (faceRemoved(facei) || faceMap_[facei] == -1)
    )
    {
        FatalErrorInFunction
            << "face " << facei
            << " already marked for removal"
            << abort(FatalError);
    }

    faces_[facei].clear();
    region_[facei] = -1;
    faceOwner_[facei] = -1;
    faceNeighbour_[facei] = -1;
    faceMap_[facei] = -1;
    if (mergeFacei >= 0)
    {
        reverseFaceMap_[facei] = -mergeFacei-2;
    }
    else
    {
        reverseFaceMap_[facei] = -1;
    }
    faceFromEdge_.erase(facei);
    faceFromPoint_.erase(facei);
    flipFaceFlux_.unset(facei);
    faceZoneFlip_.unset(facei);
    faceZone_.erase(facei);
    if (facei < faceAdditionalZones_.size())
    {
        faceAdditionalZones_[facei].clear();
    }
}


Foam::label Foam::polyTopoChange::addCell
(
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label masterCellID,
    const label zoneID
)
{
    label celli = cellMap_.size();

    if (masterPointID >= 0)
    {
        cellMap_.append(-1);
        cellFromPoint_.insert(celli, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        cellMap_.append(-1);
        cellFromEdge_.insert(celli, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        cellMap_.append(-1);
        cellFromFace_.insert(celli, masterFaceID);
    }
    else
    {
        cellMap_.append(masterCellID);
    }
    reverseCellMap_.append(celli);

    cellZone_.append(zoneID);
    // Delay extending the cellAdditionalZones

    return celli;
}


Foam::label Foam::polyTopoChange::addCell
(
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label masterCellID,
    const labelUList& zoneIDs
)
{
    label celli = cellMap_.size();

    if (masterPointID >= 0)
    {
        cellMap_.append(-1);
        cellFromPoint_.insert(celli, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        cellMap_.append(-1);
        cellFromEdge_.insert(celli, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        cellMap_.append(-1);
        cellFromFace_.insert(celli, masterFaceID);
    }
    else
    {
        cellMap_.append(masterCellID);
    }
    reverseCellMap_.append(celli);

    if (zoneIDs.size())
    {
        cellZone_.append(zoneIDs[0]);

        // Clear any additional storage
        if (celli < cellAdditionalZones_.size())
        {
            cellAdditionalZones_[celli].clear();
        }
        // (auto-vivify and) store additional zones
        for (label i = 1; i < zoneIDs.size(); ++i)
        {
            cellAdditionalZones_(celli).push_uniq(zoneIDs[i]);
        }
    }
    else
    {
        cellZone_.append(-1);
    }

    return celli;
}


void Foam::polyTopoChange::modifyCell
(
    const label celli,
    const label zoneID,
    const bool multiZone
)
{
    if (!multiZone)
    {
        cellZone_[celli] = zoneID;
        if (celli < cellAdditionalZones_.size())
        {
            cellAdditionalZones_[celli].clear();
        }
    }
    else
    {
        if (cellZone_[celli] == -1)
        {
            cellZone_[celli] = zoneID;

            // Check that no leftover additional zones
            if (cellAdditionalZones_(celli).size())
            {
                FatalErrorInFunction
                    << "Additional zones for cell:"
                    << cellAdditionalZones_(celli)
                    << abort(FatalError);
            }
        }
        else
        {
            if (zoneID == -1)
            {
                // Clear
                cellZone_[celli] = zoneID;
                if (celli < cellAdditionalZones_.size())
                {
                    cellAdditionalZones_[celli].clear();
                }
            }
            else if (zoneID != cellZone_[celli])
            {
                // New zone. Store in additional storage.
                cellAdditionalZones_(celli).push_uniq(zoneID);
            }
        }
    }
}


void Foam::polyTopoChange::modifyCell
(
    const label celli,
    const labelUList& zoneIDs
)
{
    if (celli < 0 || celli >= cellMap_.size())
    {
        FatalErrorInFunction
            << "illegal cell label " << celli << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }

    if (zoneIDs.size())
    {
        if (zoneIDs.found(-1))
        {
            FatalErrorInFunction
                << "zones cannot contain -1 : " << flatOutput(zoneIDs)
                << " for cell:" << celli
                << abort(FatalError);
        }

        cellZone_[celli] = zoneIDs[0];
        if (celli < cellAdditionalZones_.size())
        {
            cellAdditionalZones_[celli].clear();
        }
        for (label i = 1; i < zoneIDs.size(); ++i)
        {
            cellAdditionalZones_(celli).push_uniq(zoneIDs[i]);
        }
    }
    else
    {
        cellZone_[celli] = -1;
        if (celli < cellAdditionalZones_.size())
        {
            cellAdditionalZones_[celli].clear();
        }
    }
}


void Foam::polyTopoChange::removeCell
(
    const label celli,
    const label mergeCelli
)
{
    if (celli < 0 || celli >= cellMap_.size())
    {
        FatalErrorInFunction
            << "illegal cell label " << celli << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }

    if (strict_ && cellMap_[celli] == -2)
    {
        FatalErrorInFunction
            << "cell " << celli
            << " already marked for removal"
            << abort(FatalError);
    }

    cellMap_[celli] = -2;
    if (mergeCelli >= 0)
    {
        reverseCellMap_[celli] = -mergeCelli-2;
    }
    else
    {
        reverseCellMap_[celli] = -1;
    }
    cellFromPoint_.erase(celli);
    cellFromEdge_.erase(celli);
    cellFromFace_.erase(celli);
    cellZone_[celli] = -1;
    if (celli < cellAdditionalZones_.size())
    {
        cellAdditionalZones_[celli].clear();
    }
}


Foam::label Foam::polyTopoChange::pointZones
(
    const label pointi,
    DynamicList<label>& zones
) const
{
    if (pointi < 0 || pointi >= pointMap_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << pointMap_.size()-1
            << abort(FatalError);
    }
    zones.clear();

    const auto fnd = pointZone_.find(pointi);
    if (fnd)
    {
        zones.push_back(fnd());
        if (pointi < pointAdditionalZones_.size())
        {
            for (const label zonei : pointAdditionalZones_[pointi])
            {
                zones.push_uniq(zonei);
            }
        }
    }
    return zones.size();
}


Foam::label Foam::polyTopoChange::faceZones
(
    const label facei,
    DynamicList<label>& zones,
    DynamicList<bool>& flips
) const
{
    if (facei < 0 || facei >= faceMap_.size())
    {
        FatalErrorInFunction
            << "illegal face label " << facei << endl
            << "Valid face labels are 0 .. " << faceMap_.size()-1
            << abort(FatalError);
    }
    zones.clear();
    flips.clear();

    const auto fnd = faceZone_.find(facei);
    if (fnd)
    {
        zones.push_back(fnd());
        flips.push_back(flipFaceFlux_[facei]);
        if (facei < faceAdditionalZones_.size())
        {
            for (const label zoneAndSign : faceAdditionalZones_[facei])
            {
                const label zonei = mag(zoneAndSign)-1;
                if (!zones.found(zonei))
                {
                    zones.push_back(zonei);
                    flips.push_back((zoneAndSign > 0));
                }
            }
        }
    }
    return zones.size();
}


Foam::label Foam::polyTopoChange::cellZones
(
    const label celli,
    DynamicList<label>& zones
) const
{
    if (celli < 0 || celli >= cellMap_.size())
    {
        FatalErrorInFunction
            << "illegal cell label " << celli << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }
    zones.clear();
    if (cellZone_[celli] != -1)
    {
        zones.push_back(cellZone_[celli]);
    }
    if (celli < cellAdditionalZones_.size())
    {
        for (const label zonei : cellAdditionalZones_[celli])
        {
            zones.push_uniq(zonei);
        }
    }
    return zones.size();
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::changeMesh
(
    polyMesh& mesh,
    const labelUList& patchMap,
    const bool inflate,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (debug)
    {
        Pout<< "polyTopoChange::changeMesh"
            << "(polyMesh&, const bool, const bool, const bool, const bool)"
            << endl;
    }

    if (debug)
    {
        Pout<< "Old mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // new mesh points
    pointField newPoints;
    // number of internal points
    label nInternalPoints;
    // patch slicing
    labelList patchSizes;
    labelList patchStarts;
    // inflate maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromPoints;
    List<objectMap> facesFromEdges;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromPoints;
    List<objectMap> cellsFromEdges;
    List<objectMap> cellsFromFaces;
    List<objectMap> cellsFromCells;
    // old mesh info
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label>> oldFaceZoneMeshPointMaps;

    // Compact, reorder patch faces and calculate mesh/patch maps.
    compactAndReorder
    (
        mesh,
        patchMap,
        syncParallel,
        orderCells,
        orderPoints,

        nInternalPoints,
        newPoints,
        patchSizes,
        patchStarts,
        pointsFromPoints,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchStarts,
        oldFaceZoneMeshPointMaps
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Change the mesh
    // ~~~~~~~~~~~~~~~
    // This will invalidate any addressing so better make sure you have
    // all the information you need!!!

    if (inflate)
    {
        // Keep (renumbered) mesh points, store new points in map for inflation
        // (appended points (i.e. from nowhere) get value zero)
        pointField renumberedMeshPoints(newPoints.size());

        forAll(pointMap_, newPointi)
        {
            const label oldPointi = pointMap_[newPointi];

            if (oldPointi >= 0)
            {
                renumberedMeshPoints[newPointi] = mesh.points()[oldPointi];
            }
            else
            {
                renumberedMeshPoints[newPointi] = Zero;
            }
        }

        mesh.resetPrimitives
        (
            autoPtr<pointField>::New(std::move(renumberedMeshPoints)),
            autoPtr<faceList>::New(std::move(faces_)),
            autoPtr<labelList>::New(std::move(faceOwner_)),
            autoPtr<labelList>::New(std::move(faceNeighbour_)),
            patchSizes,
            patchStarts,
            syncParallel
        );

        mesh.topoChanging(true);
        // Note: could already set moving flag as well
        //       mesh.moving(true);
    }
    else
    {
        // Set new points.
        mesh.resetPrimitives
        (
            autoPtr<pointField>::New(std::move(newPoints)),
            autoPtr<faceList>::New(std::move(faces_)),
            autoPtr<labelList>::New(std::move(faceOwner_)),
            autoPtr<labelList>::New(std::move(faceNeighbour_)),
            patchSizes,
            patchStarts,
            syncParallel
        );
        mesh.topoChanging(true);
    }

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nAdd, nInflate, nMerge, nRemove;
        countMap(pointMap_, reversePointMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Points:"
            << "  added(from point):" << nAdd
            << "  added(from nothing):" << nInflate
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other cell):" << nMerge
            << "  removed:" << nRemove
            << nl
            << endl;
    }

    // Zones
    // ~~~~~

    // Inverse of point/face/cell zone addressing.
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, mesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // Clear zone info
    {
        pointZone_.clearStorage();
        faceZone_.clearStorage();
        faceZoneFlip_.clearStorage();
        cellZone_.clearStorage();
    }


    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(mesh.boundaryMesh().size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        patchMap,
        mesh.boundaryMesh(),
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(mesh.faceZones().size());
    calcFaceZonePointMap(mesh, oldFaceZoneMeshPointMaps, faceZonePointMap);

    labelHashSet flipFaceFluxSet(HashSetOps::used(flipFaceFlux_));

    return autoPtr<mapPolyMesh>::New
    (
        mesh,
        nOldPoints,
        nOldFaces,
        nOldCells,

        pointMap_,
        pointsFromPoints,

        faceMap_,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,

        cellMap_,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,

        reversePointMap_,
        reverseFaceMap_,
        reverseCellMap_,

        flipFaceFluxSet,

        patchPointMap,

        pointZoneMap,

        faceZonePointMap,
        faceZoneFaceMap,
        cellZoneMap,

        newPoints,          // if empty signals no inflation.
        oldPatchStarts,
        oldPatchNMeshPoints,

        oldCellVolumes,

        true                // steal storage.
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::changeMesh
(
    polyMesh& mesh,
    const bool inflate,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    return changeMesh
    (
        mesh,
        identity(mesh.boundaryMesh().size()),
        inflate,
        syncParallel,
        orderCells,
        orderPoints
    );
}


// ************************************************************************* //
