/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Simple timing tests for some polyMesh primitives

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clockTime.H"
#include "Time.H"
#include "PDRblock.H"
#include "polyMesh.H"
#include "ListOps.H"

using namespace Foam;

void printAlloc(const polyMesh& mesh)
{
    Info<< "memory"
        << " hasCellPoints:" << mesh.hasCellPoints()
        << " hasPointCells:" << mesh.hasPointCells() << endl;
}


void printInfo(const polyMesh& mesh)
{
    Info<< "polyMesh"
        << " nPoints:" << mesh.nPoints()
        << " nInternalFaces:" << mesh.nInternalFaces()
        << " nFaces:" << mesh.nFaces()
        << " nCells:" << mesh.nCells() << endl;
}


// How point cells are calculated in OpenFOAM-v2212 and earlier
autoPtr<labelListList> pointCells_2212(const polyMesh& mesh)
{
    const cellList& cf = mesh.cells();

    // Count number of cells per point

    labelList npc(mesh.nPoints(), Zero);

    forAll(cf, celli)
    {
        const labelList curPoints = cf[celli].labels(mesh.faces());

        for (const label pointi : curPoints)
        {
            ++npc[pointi];
        }
    }


    // Size and fill cells per point

    auto pcPtr_ = autoPtr<labelListList>::New(npc.size());
    labelListList& pointCellAddr = *pcPtr_;

    forAll(pointCellAddr, pointi)
    {
        pointCellAddr[pointi].setSize(npc[pointi]);
        npc[pointi] = 0;
    }


    forAll(cf, celli)
    {
        const labelList curPoints = cf[celli].labels(mesh.faces());

        for (const label pointi : curPoints)
        {
            pointCellAddr[pointi][npc[pointi]++] = celli;
        }
    }

    return pcPtr_;
}


// Line cell::labels but with persistent storage
void cell_labels
(
    const cell& cFaces,
    const faceUList& meshFaces,
    DynamicList<label>& pointLabels
)
{
    // const labelList& cFaces = *this;

    label nVerts = 0;
    for (const label facei : cFaces)
    {
        nVerts += meshFaces[facei].size();
    }

    pointLabels.resize(pointLabels.capacity());  // Use full storage

    // The first face has no duplicates, can copy in values
    const labelList& firstFace = meshFaces[cFaces[0]];

    std::copy(firstFace.cbegin(), firstFace.cend(), pointLabels.begin());

    // Now already contains some vertices
    nVerts = firstFace.size();

    // For the rest of the faces. For each vertex, check if the point is
    // already inserted (up to nVerts, which now carries the number of real
    // points. If not, add it at the end of the list.

    for (label facei = 1; facei < cFaces.size(); ++facei)
    {
        for (const label curPoint : meshFaces[cFaces[facei]])
        {
            bool pointFound = false;

            for (label checki = 0; checki < nVerts; ++checki)
            {
                if (curPoint == pointLabels[checki])
                {
                    pointFound = true;
                    break;
                }
            }

            if (!pointFound)
            {
                pointLabels[nVerts] = curPoint;
                ++nVerts;
            }
        }
    }

    pointLabels.resize(nVerts);
}



// Like OpenFOAM-v2212, but with cell::labels unrolled to avoid allocations
autoPtr<labelListList> pointCells_2212mod(const polyMesh& mesh)
{
    const cellList& cf = mesh.cells();

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    // Count number of cells per point

    labelList npc(mesh.nPoints(), Zero);

    for (const cell& c : cf)
    {
        cell_labels(c, mesh.faces(), vertices);

        for (const label pointi : vertices)
        {
            ++npc[pointi];
        }
    }


    // Size and fill cells per point

    auto pcPtr_ = autoPtr<labelListList>::New(npc.size());
    labelListList& pointCellAddr = *pcPtr_;

    forAll(pointCellAddr, pointi)
    {
        pointCellAddr[pointi].resize(npc[pointi]);
        npc[pointi] = 0;
    }


    forAll(cf, celli)
    {
        cell_labels(cf[celli], mesh.faces(), vertices);

        for (const label pointi : vertices)
        {
            pointCellAddr[pointi][npc[pointi]++] = celli;
        }
    }

    return pcPtr_;
}


// How cells points are calculated in OpenFOAM-v2212 and earlier
autoPtr<labelListList> cellPoints_2212(const polyMesh& mesh)
{
    autoPtr<labelListList> pointCells = pointCells_2212(mesh);

    auto cpPtr_ = autoPtr<labelListList>::New(mesh.nCells());

    invertManyToMany(mesh.nCells(), pointCells(), *cpPtr_);

    return cpPtr_;
}


// Calculate with bitSet tracking and avoid cells::labels
autoPtr<labelListList> pointCells_bitSet(const polyMesh& mesh)
{
    // Calculate point-cell topology

    const cellList& cellLst = mesh.cells();
    const faceList& faceLst = mesh.faces();

    // For tracking (only use each point id once)
    bitSet usedPoints(mesh.nPoints());

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    const label loopLen = mesh.nCells();

    // Step 1: count number of cells per point

    labelList pointCount(mesh.nPoints(), Zero);

    for (label celli = 0; celli < loopLen; ++celli)
    {
        // Clear any previous contents
        usedPoints.unset(vertices);
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (usedPoints.set(pointi))
                {
                    vertices.push_back(pointi);
                    ++pointCount[pointi];
                }
            }
        }
    }


    // Step 2: set sizing, reset counters

    auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());
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
        usedPoints.unset(vertices);
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (usedPoints.set(pointi))
                {
                    vertices.push_back(pointi);
                    pointCellAddr[pointi][pointCount[pointi]++] = celli;
                }
            }
        }
    }

    return pcPtr_;
}


// Calculate with bitSet tracking and avoid cells::labels
autoPtr<labelListList> cellPoints_bitSet(const polyMesh& mesh)
{
    // Calculate cell-point topology

    auto cpPtr_ = autoPtr<labelListList>::New(mesh.nCells());
    auto& cellPointAddr = *cpPtr_;

    const cellList& cellLst = mesh.cells();
    const faceList& faceLst = mesh.faces();

    // For tracking (only use each point id once)
    bitSet usedPoints(mesh.nPoints());

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    const label loopLen = mesh.nCells();

    for (label celli = 0; celli < loopLen; ++celli)
    {
        // Clear any previous contents
        usedPoints.unset(vertices);
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (usedPoints.set(pointi))
                {
                    vertices.push_back(pointi);
                }
            }
        }

        cellPointAddr[celli] = vertices;  // unsorted
    }

    return cpPtr_;
}


// Calculate with linear lookup and avoid cells::labels
autoPtr<labelListList> pointCells_linear(const polyMesh& mesh)
{
    // Calculate point-cell topology

    const cellList& cellLst = mesh.cells();
    const faceList& faceLst = mesh.faces();

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    const label loopLen = mesh.nCells();

    // Step 1: count number of cells per point

    labelList pointCount(mesh.nPoints(), Zero);

    for (label celli = 0; celli < loopLen; ++celli)
    {
        // Clear any previous contents
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (!vertices.contains(pointi))
                {
                    vertices.push_back(pointi);
                    ++pointCount[pointi];
                }
            }
        }
    }


    // Step 2: set sizing, reset counters

    auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());
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
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (!vertices.contains(pointi))
                {
                    vertices.push_back(pointi);
                    pointCellAddr[pointi][pointCount[pointi]++] = celli;
                }
            }
        }
    }

    return pcPtr_;
}


// Calculate with linear lookup and avoid cells::labels
autoPtr<labelListList> cellPoints_linear(const polyMesh& mesh)
{
    // Calculate cell-point topology

    auto cpPtr_ = autoPtr<labelListList>::New(mesh.nCells());
    auto& cellPointAddr = *cpPtr_;

    const cellList& cellLst = mesh.cells();
    const faceList& faceLst = mesh.faces();

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    const label loopLen = mesh.nCells();

    for (label celli = 0; celli < loopLen; ++celli)
    {
        // Clear any previous contents
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                // Only once for each point id
                if (!vertices.contains(pointi))
                {
                    vertices.push_back(pointi);
                }
            }
        }

        cellPointAddr[celli] = vertices;  // unsorted
    }

    return cpPtr_;
}


// Calculate point-cell from point-face information
autoPtr<labelListList> pointCells_faces(const polyMesh& mesh)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const labelListList& pFaces = mesh.pointFaces();

    const label loopLen = mesh.nPoints();

    auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());
    auto& pointCellAddr = *pcPtr_;

    DynamicList<label> storage(256);

    for (label pointi = 0; pointi < loopLen; ++pointi)
    {
        // Clear any previous contents
        storage.clear();

        for (const label facei : pFaces[pointi])
        {
            // Owner cell
            storage.push_back(own[facei]);

            // Neighbour cell
            if (facei < mesh.nInternalFaces())
            {
                storage.push_back(nei[facei]);
            }
        }

        // Sort + unique to eliminate duplicates
        std::sort(storage.begin(), storage.end());
        auto last = std::unique(storage.begin(), storage.end());
        storage.resize(label(last - storage.begin()));

        pointCellAddr[pointi] = storage;
    }

    return pcPtr_;
}

// Calculate point-cell from point-face information
autoPtr<labelListList> pointCells_bitSet_faces(const polyMesh& mesh)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const labelListList& pFaces = mesh.pointFaces();

    const label loopLen = mesh.nPoints();

    auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());
    auto& pointCellAddr = *pcPtr_;

    // For tracking (only use each cell id once)
    bitSet usedCells(mesh.nCells());

    DynamicList<label> storage(256);

    for (label pointi = 0; pointi < loopLen; ++pointi)
    {
        // Clear any previous contents
        usedCells.unset(storage);
        storage.clear();

        for (const label facei : pFaces[pointi])
        {
            // Owner cell - only once
            if (usedCells.set(own[facei]))
            {
                storage.push_back(own[facei]);
            }

            // Neighbour cell
            if (facei < mesh.nInternalFaces() && usedCells.set(nei[facei]))
            {
                storage.push_back(nei[facei]);
            }
        }

        pointCellAddr[pointi] = storage;
    }

    return pcPtr_;
}


// Calculate point-cell from cell-point information
autoPtr<labelListList> pointCells_bitSet_alon(const polyMesh& mesh)
{
    autoPtr<labelListList> cellPoints = cellPoints_bitSet(mesh);

    auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());

    invertManyToMany(mesh.nPoints(), cellPoints(), *pcPtr_);

    return pcPtr_;
}


// Eliminate duplicates with sort+unique
autoPtr<labelListList> cellPoints_sorted(const polyMesh& mesh)
{
    // Calculate cell-point topology

    auto cpPtr_ = autoPtr<labelListList>::New(mesh.nCells());
    auto& cellPointAddr = *cpPtr_;

    const cellList& cellLst = mesh.cells();
    const faceList& faceLst = mesh.faces();

    // Vertex labels for the current cell
    DynamicList<label> vertices(256);

    const label loopLen = mesh.nCells();

    for (label celli = 0; celli < loopLen; ++celli)
    {
        // Clear any previous contents
        vertices.clear();

        for (const label facei : cellLst[celli])
        {
            for (const label pointi : faceLst[facei])
            {
                vertices.push_back(pointi);
            }
        }

        // Sort + unique to eliminate duplicates
        std::sort(vertices.begin(), vertices.end());
        auto last = std::unique(vertices.begin(), vertices.end());
        vertices.resize(label(last - vertices.begin()));

        cellPointAddr[celli] = vertices;
    }

    return cpPtr_;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addOption("nCells", "number", "The number of cells");

    #include "setRootCase.H"

    const scalar cellCount(args.getOrDefault<scalar>("nCells", 1000));

    const label nDivs(::round(::cbrt(cellCount)));

    PDRblock blkMesh(boundBox(zero_one{}), labelVector::uniform(nDivs));

    autoPtr<Time> dummyTimePtr(Time::New());

    Info<< "Requested " << cellCount
        << " cells, blockMesh with " << blkMesh.nCells() << " cells" << nl;

    autoPtr<polyMesh> meshPtr = blkMesh.innerMesh
    (
        IOobject
        (
            "Testing",
            dummyTimePtr->system(),
            *dummyTimePtr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto& mesh = meshPtr();

    printInfo(mesh);
    printAlloc(mesh);

    clockTime timing;

    // pointCells
    {
        mesh.clearOut();
        timing.resetTime();
        (void) mesh.pointCells();
        Info<< "pointCells (builtin): " << timing.elapsedTime() << " s" << nl;
    }

    // cellPoints
    {
        mesh.clearOut();
        timing.resetTime();
        (void) mesh.cellPoints();
        Info<< "cellPoints (builtin): " << timing.elapsedTime() << " s" << nl;
    }
    Info<< nl;


    // pointCells
    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_2212(mesh);
        Info<< "pointCells (2212): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_2212mod(mesh);
        Info<< "pointCells (2212mod): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_bitSet(mesh);
        Info<< "pointCells (bitSet): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_linear(mesh);
        Info<< "pointCells (linear): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_faces(mesh);
        Info<< "pointCells (faces): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_bitSet_faces(mesh);
        Info<< "pointCells (bitSet faces): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) pointCells_bitSet_alon(mesh);
        Info<< "pointCells (bitSet alon): " << timing.elapsedTime() << " s" << nl;
    }

    // cellPoints
    {
        mesh.clearOut();
        timing.resetTime();
        (void) cellPoints_2212(mesh);
        Info<< "cellPoints (2212): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) cellPoints_bitSet(mesh);
        Info<< "cellPoints (bitSet): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) cellPoints_linear(mesh);
        Info<< "cellPoints (linear): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) cellPoints_sorted(mesh);
        Info<< "cellPoints (sorted): " << timing.elapsedTime() << " s" << nl;
    }


    // With precalculated values
    {
        mesh.clearOut();
        const auto& cp = mesh.cellPoints();
        timing.resetTime();

        auto pcPtr_ = autoPtr<labelListList>::New(mesh.nPoints());

        invertManyToMany(mesh.nPoints(), cp, *pcPtr_);

        Info<< "pointCells (from cached cellPoints): " << timing.elapsedTime() << " s" << nl;
    }

    // With precalculated values
    {
        mesh.clearOut();
        (void)mesh.pointFaces();
        timing.resetTime();

        (void) pointCells_bitSet_faces(mesh);

        Info<< "pointCells (bitSet from cached pointFaces): " << timing.elapsedTime() << " s" << nl;
    }


    // With precalculated values
    {
        mesh.clearOut();
        const auto& pc = mesh.pointCells();
        timing.resetTime();

        auto cpPtr_ = autoPtr<labelListList>::New(mesh.nCells());

        invertManyToMany(mesh.nCells(), pc, *cpPtr_);

        Info<< "cellPoints (from cached pointCells): " << timing.elapsedTime() << " s" << nl;
    }


    // Re-measure timings
    Info<< nl;
    {
        mesh.clearOut();
        timing.resetTime();
        (void) mesh.pointCells();
        Info<< "pointCells (builtin): " << timing.elapsedTime() << " s" << nl;
    }

    {
        mesh.clearOut();
        timing.resetTime();
        (void) mesh.cellPoints();
        Info<< "cellPoints (builtin): " << timing.elapsedTime() << " s" << nl;
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
