/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "primitiveMesh.H"
#include "cell.H"
#include "bitSet.H"
#include "DynamicList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcPointCells() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcPointCells() : "
            << "calculating pointCells"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate pointCells
    // if the pointer is already set
    if (pcPtr_)
    {
        FatalErrorInFunction
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else if (hasCellPoints())
    {
        // Invert cellPoints
        pcPtr_ = new labelListList(nPoints());
        invertManyToMany(nPoints(), cellPoints(), *pcPtr_);
    }
    else if (hasPointFaces())
    {
        // Calculate point-cell from point-face information

        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelListList& pFaces = pointFaces();

        // Tracking (only use each cell id once)
        bitSet usedCells(nCells());

        // Cell ids for the point currently being processed
        DynamicList<label> currCells(256);

        const label loopLen = nPoints();

        pcPtr_ = new labelListList(nPoints());
        auto& pointCellAddr = *pcPtr_;

        for (label pointi = 0; pointi < loopLen; ++pointi)
        {
            // Clear any previous contents
            usedCells.unset(currCells);
            currCells.clear();

            for (const label facei : pFaces[pointi])
            {
                // Owner cell - only allow one occurance
                if (usedCells.set(own[facei]))
                {
                    currCells.push_back(own[facei]);
                }

                // Neighbour cell - only allow one occurance
                if (facei < nInternalFaces())
                {
                    if (usedCells.set(nei[facei]))
                    {
                        currCells.push_back(nei[facei]);
                    }
                }
            }

            pointCellAddr[pointi] = currCells;  // NB: unsorted
        }
    }
    else
    {
        // Calculate point-cell topology

        const cellList& cellLst = cells();
        const faceList& faceLst = faces();

        // Tracking (only use each point id once)
        bitSet usedPoints(nPoints());

        // Which of usedPoints needs to be unset [faster]
        DynamicList<label> currPoints(256);

        const label loopLen = nCells();

        // Step 1: count number of cells per point

        labelList pointCount(nPoints(), Zero);

        for (label celli = 0; celli < loopLen; ++celli)
        {
            // Clear any previous contents
            usedPoints.unset(currPoints);
            currPoints.clear();

            for (const label facei : cellLst[celli])
            {
                for (const label pointi : faceLst[facei])
                {
                    // Only once for each point id
                    if (usedPoints.set(pointi))
                    {
                        currPoints.push_back(pointi);  // Needed for cleanup
                        ++pointCount[pointi];
                    }
                }
            }
        }


        // Step 2: set sizing, reset counters

        pcPtr_ = new labelListList(nPoints());
        auto& pointCellAddr = *pcPtr_;

        forAll(pointCellAddr, pointi)
        {
            pointCellAddr[pointi].resize_nocopy(pointCount[pointi]);
            pointCount[pointi] = 0;
        }


        // Step 3: fill in values. Logic as per step 1
        for (label celli = 0; celli < loopLen; ++celli)
        {
            // Clear any previous contents
            usedPoints.unset(currPoints);
            currPoints.clear();

            for (const label facei : cellLst[celli])
            {
                for (const label pointi : faceLst[facei])
                {
                    // Only once for each point id
                    if (usedPoints.set(pointi))
                    {
                        currPoints.push_back(pointi);  // Needed for cleanup
                        pointCellAddr[pointi][pointCount[pointi]++] = celli;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::pointCells() const
{
    if (!pcPtr_)
    {
        calcPointCells();
    }

    return *pcPtr_;
}


const Foam::labelList& Foam::primitiveMesh::pointCells
(
    const label pointi,
    DynamicList<label>& storage
) const
{
    if (hasPointCells())
    {
        return pointCells()[pointi];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelList& pFaces = pointFaces()[pointi];

        storage.clear();

        for (const label facei : pFaces)
        {
            // Owner cell
            storage.push_back(own[facei]);

            // Neighbour cell
            if (facei < nInternalFaces())
            {
                storage.push_back(nei[facei]);
            }
        }

        // Filter duplicates
        if (storage.size() > 1)
        {
            std::sort(storage.begin(), storage.end());
            auto last = std::unique(storage.begin(), storage.end());
            storage.resize(label(last - storage.begin()));
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::pointCells(const label pointi) const
{
    return pointCells(pointi, labels_);
}


// ************************************************************************* //
