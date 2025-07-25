/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "primitiveMesh.H"
#include "cell.H"
#include "bitSet.H"
#include "DynamicList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellPoints() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::cellCellPoints() : "
            << "calculating cellPoints" << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate cellPoints
    // if the pointer is already set
    if (cpPtr_)
    {
        FatalErrorInFunction
            << "cellPoints already calculated"
            << abort(FatalError);
    }
    else if (hasPointCells())
    {
        // Invert pointCells
        cpPtr_ = new labelListList(nCells());
        invertManyToMany(nCells(), pointCells(), *cpPtr_);
    }
    else
    {
        // Calculate cell-point topology

        cpPtr_ = new labelListList(nCells());
        auto& cellPointAddr = *cpPtr_;

        const cellList& cellLst = cells();
        const faceList& faceLst = faces();

        // Tracking (only use each point id once)
        bitSet usedPoints(nPoints());

        // Vertex labels for the current cell
        DynamicList<label> currPoints(256);

        const label loopLen = nCells();

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
                        currPoints.push_back(pointi);
                    }
                }
            }

            cellPointAddr[celli] = currPoints;  // NB: unsorted
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::cellPoints() const
{
    if (!cpPtr_)
    {
        calcCellPoints();
    }

    return *cpPtr_;
}


const Foam::labelList& Foam::primitiveMesh::cellPoints
(
    const label celli,
    labelHashSet& set,
    DynamicList<label>& storage
) const
{
    if (hasCellPoints())
    {
        return cellPoints()[celli];
    }

    const faceList& fcs = faces();
    const labelList& cFaces = cells()[celli];

    set.clear();

    for (const label facei : cFaces)
    {
        set.insert(fcs[facei]);
    }

    storage.clear();
    if (storage.capacity() < set.size())
    {
        storage.setCapacity(set.size());
    }

    for (const label pointi : set)
    {
        storage.push_back(pointi);
    }

    return storage;
}


const Foam::labelList& Foam::primitiveMesh::cellPoints(const label celli) const
{
    return cellPoints(celli, labelSet_, labels_);
}


// ************************************************************************* //
