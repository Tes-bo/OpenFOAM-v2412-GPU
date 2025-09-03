/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "globalMeshData.H"
#include "globalIndex.H"
#include "CompactListList.H"
#include "processorPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// NOTE: the AgglomerationType is anything that behaves like a List with
// an operator[] and provides coverage in the (0-nCells) range.
// - identityOp() does this

template<class AgglomerationType>
static void calcCellCellsImpl
(
    const polyMesh& mesh,
    const AgglomerationType& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells,
    CompactListList<scalar>* cellCellWeightsPtr = nullptr
)
{
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Global cell numbers (agglomerated numbering)
    const label myProci = UPstream::myProcNo(UPstream::worldComm);
    const globalIndex globalAgglom(nLocalCoarse, UPstream::worldComm, parallel);

    // The agglomerated owner per boundary faces (global numbering)
    // from the other side (of coupled patches)

    labelList globalNeighbour;
    {
        const label myAgglomOffset = globalAgglom.localStart(myProci);

        const labelList::subList bndFaceOwner = pbm.faceOwner();

        const label nBoundaryFaces = bndFaceOwner.size();

        globalNeighbour.resize(nBoundaryFaces);

        for (label bfacei = 0; bfacei < nBoundaryFaces; ++bfacei)
        {
            label val = agglom[bndFaceOwner[bfacei]];
            if (val >= 0)
            {
                // Only offset 'real' (non-negative) agglomerations
                val += myAgglomOffset;
            }
            globalNeighbour[bfacei] = val;
        }
    }

    // Swap boundary neighbour information:
    // - cyclics and (optionally) processor
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour, parallel);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per coarse cell
    labelList nFacesPerCell(nLocalCoarse, Foam::zero{});

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        // Negative agglomeration (exclude from subset)
        if (own < 0 || nei < 0) continue;

        ++nFacesPerCell[own];
        ++nFacesPerCell[nei];
    }

    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            const label bndOffset = mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[facei-bndOffset];

                // Negative agglomeration (exclude from subset)
                if (own < 0 || globalNei < 0) continue;

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    ++nFacesPerCell[own];
                }
            }
        }
    }


    // Fill in offset and data
    // ~~~~~~~~~~~~~~~~~~~~~~~

    cellCells.resize_nocopy(nFacesPerCell);
    nFacesPerCell = 0;  // Restart the count


    // CSR indexing
    const labelList& offsets = cellCells.offsets();

    // CSR connections
    labelList& connect = cellCells.values();


    // CSR connection weights
    scalarList weights;
    if (cellCellWeightsPtr)
    {
        cellCellWeightsPtr->clear();
        weights.resize(cellCells.totalSize());
    }


    // For internal faces is just offsetted owner and neighbour
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        // Negative agglomeration (exclude from subset)
        if (own < 0 || nei < 0) continue;

        const label ownIndex = offsets[own] + nFacesPerCell[own]++;
        const label neiIndex = offsets[nei] + nFacesPerCell[nei]++;

        connect[ownIndex] = globalAgglom.toGlobal(myProci, nei);
        connect[neiIndex] = globalAgglom.toGlobal(myProci, own);

        if (!weights.empty())
        {
            weights[ownIndex] = Foam::mag(mesh.faceAreas()[facei]);
            weights[neiIndex] = weights[ownIndex];
        }
    }

    // For boundary faces is offsetted coupled neighbour
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            const label bndOffset = mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[facei-bndOffset];

                // Negative agglomeration (exclude from subset)
                if (own < 0 || globalNei < 0) continue;

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    const label ownIndex = offsets[own] + nFacesPerCell[own]++;

                    connect[ownIndex] = globalNei;

                    if (!weights.empty())
                    {
                        weights[ownIndex] = Foam::mag(mesh.faceAreas()[facei]);
                    }
                }
            }
        }
    }


    // Check for duplicates connections between cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done as postprocessing step since we now have cellCells.

    // NB: Because of agglomeration, self-connections will occur
    //     and must be filtered out.

    if (!cellCells.empty())
    {
        // Need non-const access to CSR indexing
        labelList& offsets = cellCells.offsets();

        label startIndex = offsets[0];
        label newIndex = 0;
        labelHashSet nbrCells;

        const label nCellCells = cellCells.size();

        for (label celli = 0; celli < nCellCells; ++celli)
        {
            // Could be done as combination of std::sort, std::unique_copy
            // except need to block out 'self' and std::copy_if/remove_if
            // are undefined for overlapping regions

            const label self = globalAgglom.toGlobal(myProci, celli);

            nbrCells.clear();
            nbrCells.insert(self);

            const label endIndex = offsets[celli+1];

            for (label i = startIndex; i < endIndex; ++i)
            {
                if (nbrCells.insert(connect[i]))
                {
                    connect[newIndex] = connect[i];

                    if (!weights.empty())
                    {
                        weights[newIndex] = weights[i];
                    }

                    ++newIndex;
                }
            }
            startIndex = endIndex;
            offsets[celli+1] = newIndex;
        }

        connect.resize(newIndex);
        if (!weights.empty())
        {
            weights.resize(newIndex);
        }
    }


    // CSR connection weights
    // - addressing is identical to the connections
    if (cellCellWeightsPtr)
    {
        cellCellWeightsPtr->offsets() = cellCells.offsets();
        cellCellWeightsPtr->values() = std::move(weights);
    }


    //forAll(cellCells, celli)
    //{
    //    Pout<< "Original: Coarse cell " << celli << endl;
    //    forAll(mesh.cellCells()[celli], i)
    //    {
    //        Pout<< "    nbr:" << mesh.cellCells()[celli][i] << endl;
    //    }
    //    Pout<< "Compacted: Coarse cell " << celli << endl;
    //    const labelUList cCells = cellCells[celli];
    //    forAll(cCells, i)
    //    {
    //        Pout<< "    nbr:" << cCells[i] << endl;
    //    }
    //}
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const labelUList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells
)
{
    calcCellCellsImpl
    (
        mesh,
        agglom,
        nLocalCoarse,
        parallel,
        cellCells
    );
}


void Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const labelUList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells,
    CompactListList<scalar>& cellCellWeights
)
{
    calcCellCellsImpl
    (
        mesh,
        agglom,
        nLocalCoarse,
        parallel,
        cellCells,
       &cellCellWeights
    );
}


// Convenience forms

void Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    CompactListList<label>& cellCells,
    const bool parallel
)
{
    calcCellCellsImpl
    (
        mesh,
        Foam::identityOp{},
        mesh.nCells(),
        parallel,
        cellCells
    );
}


Foam::labelList Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const bitSet& selectedCells,
    CompactListList<label>& cellCells,
    const bool parallel
)
{
    const label nCells = mesh.nCells();

    labelList agglom(nCells, -1);
    labelList cellMap;

    // First pass - sorted order without duplicates
    label nCompact = 0;
    for (const label celli : selectedCells)
    {
        // A bitSet has no negatives/duplicates, so just check the upper range
        if (celli >= nCells)
        {
            break;
        }
        else
        {
            agglom[celli] = celli;
            ++nCompact;
        }
    }

    // Second pass - finalize mappings
    if (nCompact)
    {
        cellMap.resize(nCompact);
        nCompact = 0;

        for (label& celli : agglom)
        {
            if (celli >= 0)
            {
                cellMap[nCompact] = celli;
                celli = nCompact;
                ++nCompact;

                if (nCompact == cellMap.size())
                {
                    break;  // Early termination
                }
            }
        }
    }

    globalMeshData::calcCellCells
    (
        mesh,
        agglom,
        nCompact,   // == cellMap.size()
        parallel,
        cellCells
    );

    return cellMap;
}


Foam::labelList Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const labelUList& selectedCells,
    CompactListList<label>& cellCells,
    const bool parallel
)
{
    const label nCells = mesh.nCells();

    labelList agglom(nCells, -1);
    labelList cellMap;

    // First pass - creates a sorted order without duplicates
    label nCompact = 0;
    for (const label celli : selectedCells)
    {
        // Check cell is in range, and squash out duplicates
        if (celli >= 0 && celli < nCells && agglom[celli] < 0)
        {
            agglom[celli] = celli;
            ++nCompact;
        }
    }

    // Second pass - finalize mappings
    if (nCompact)
    {
        cellMap.resize(nCompact);
        nCompact = 0;

        for (label& celli : agglom)
        {
            if (celli >= 0)
            {
                cellMap[nCompact] = celli;
                celli = nCompact;
                ++nCompact;

                if (nCompact == cellMap.size())
                {
                    break;  // Early termination
                }
            }
        }
    }

    globalMeshData::calcCellCells
    (
        mesh,
        agglom,
        nCompact,   // == cellMap.size()
        parallel,
        cellCells
    );

    return cellMap;
}


// ************************************************************************* //
