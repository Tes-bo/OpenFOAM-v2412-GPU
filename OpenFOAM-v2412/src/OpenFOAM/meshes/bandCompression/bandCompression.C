/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "bandCompression.H"
#include "bitSet.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "CircularBuffer.H"
#include "CompactListList.H"
#include "DynamicList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Process connections with the Cuthill-McKee algorithm.
// The connections are CompactListList<label> or a labelListList.
template<class ConnectionListListType>
Foam::labelList cuthill_mckee_algorithm
(
    const ConnectionListListType& cellCellAddressing
)
{
    using namespace Foam;

    const label nOldCells(cellCellAddressing.size());

    // Which cells are visited/unvisited
    bitSet unvisited(nOldCells, true);

    // The new output order
    labelList newOrder(nOldCells);


    // Various work arrays
    // ~~~~~~~~~~~~~~~~~~~

    // Neighbour cells
    DynamicList<label> nbrCells;

    // Neighbour ordering
    DynamicList<label> nbrOrder;

    // Corresponding weights for neighbour cells
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

                // Add in increasing order of connectivity

                // 1. Count neighbours of unvisited neighbours
                nbrCells.clear();
                weights.clear();

                // Find if the neighbours have been visited
                const auto& neighbours = cellCellAddressing[currCelli];

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

    // Debug:
    // - the peak capacity of queuedCells approximates the
    //   maximum intermediate bandwidth
    #if 0
    Pout<< "bandCompression: peak-capacity=" << queuedCells.capacity() << nl;
    #endif

    // Now we have new-to-old in newOrder.
    return newOrder;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelList Foam::meshTools::bandCompression
(
    const labelUList& cellCells,
    const labelUList& offsets
)
{
    // Protect against zero-sized offset list
    const label nOldCells = Foam::max(0, (offsets.size()-1));

    // Count number of neighbours
    labelList numNbrs(nOldCells, Foam::zero{});
    for (label celli = 0; celli < nOldCells; ++celli)
    {
        const label beg = offsets[celli];
        const label end = offsets[celli+1];

        for (label idx = beg; idx < end; ++idx)
        {
            ++numNbrs[celli];
            ++numNbrs[cellCells[idx]];
        }
    }


    // Which cells are visited/unvisited
    bitSet unvisited(nOldCells, true);

    // The new output order
    labelList newOrder(nOldCells);


    // Various work arrays
    // ~~~~~~~~~~~~~~~~~~~

    // Neighbour cells
    DynamicList<label> nbrCells;

    // Neighbour ordering
    DynamicList<label> nbrOrder;

    // Corresponding weights for neighbour cells
    DynamicList<label> weights;

    // FIFO buffer for renumbering.
    CircularBuffer<label> queuedCells(1024);


    label cellInOrder = 0;

    while (true)
    {
        // Find lowest connected cell that has not been visited yet
        label currCelli = -1;
        label minCount = labelMax;

        for (const label celli : unvisited)
        {
            const label nbrCount = numNbrs[celli];

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


        // Starting from currCellii - walk breadth-first

        queuedCells.push_back(currCelli);

        // loop through the nextCell list. Add the first cell into the
        // cell order if it has not already been visited and ask for its
        // neighbours. If the neighbour in question has not been visited,
        // add it to the end of the nextCell list

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

                // Add in increasing order of connectivity

                // 1. Count neighbours of unvisited neighbours
                nbrCells.clear();
                weights.clear();

                const label beg = offsets[currCelli];
                const label end = offsets[currCelli+1];

                for (label idx = beg; idx < end; ++idx)
                {
                    const label nbr = cellCells[idx];
                    const label nbrCount = numNbrs[nbr];

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

    // Debug:
    // - the peak capacity of queuedCells approximates the
    //   maximum intermediate bandwidth
    #if 0
    Pout<< "bandCompression: peak-capacity=" << queuedCells.capacity() << nl;
    #endif

    // Now we have new-to-old in newOrder.
    return newOrder;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelList Foam::meshTools::bandCompression(const polyMesh& mesh)
{
    // Local mesh connectivity
    CompactListList<label> cellCells;
    globalMeshData::calcCellCells(mesh, cellCells);

    return cuthill_mckee_algorithm(cellCells);
}


Foam::labelList Foam::meshTools::bandCompression
(
    const CompactListList<label>& cellCellAddressing
)
{
    return cuthill_mckee_algorithm(cellCellAddressing);
}


Foam::labelList Foam::meshTools::bandCompression
(
    const labelListList& cellCellAddressing
)
{
    return cuthill_mckee_algorithm(cellCellAddressing);
}


// ************************************************************************* //
