/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "triangulatedPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triangulatedPatch::randomPoint
(
    Random& rnd,
    const scalar c,
    point& result,
    label& facei,
    label& celli
) const
{
    result = point::min;
    facei = -1;
    celli = -1;

    if (triWght_.empty() || c < triWght_.front() || c > triWght_.back())
    {
        return false;
    }

    // Find corresponding decomposed face triangle
    // Note: triWght_ is sized nTri+1 (zero added at start)
    //
    // TBD: binary search with findLower(triWght_, c) ??
    label trii = 0;
    for (label i = 0; i < triWght_.size() - 1; ++i)
    {
        if (c > triWght_[i] && c <= triWght_[i+1])
        {
            trii = i;
            break;
        }
    }

    // Find random point in triangle
    const pointField& points = patch_.points();

    result = triFace_[trii].tri(points).randomPoint(rnd);
    facei = triFace_[trii].index();
    celli = patch_.faceCells()[facei];

    if (perturbTol_ > 0)
    {
        const polyMesh& mesh = patch_.boundaryMesh().mesh();
        const point& cc = mesh.cellCentres()[celli];

        // Normal points out of domain => subtract correction
        const vector& n = patch_.faceNormals()[facei];
        result -= perturbTol_*n*mag(n & (cc - result));

        // Reset facei - point no longer resides on the face
        facei = -1;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triangulatedPatch::triangulatedPatch
(
    const polyPatch& patch,
    const scalar perturbTol
)
:
    patch_(patch),
    perturbTol_(perturbTol),
    triFace_(),
    triWght_()
{
    update();
}


Foam::triangulatedPatch::triangulatedPatch
(
    const polyMesh& mesh,
    const word& patchName,
    const scalar perturbTol
)
:
    triangulatedPatch(mesh.boundaryMesh()[patchName], perturbTol)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::triangulatedPatch::triangulate
(
    const polyPatch& pp,
    List<labelledTri>& tris
)
{
    const pointField& points = pp.points();

    // Triangulate the patch faces and create addressing
    label nTris = 0;
    for (const face& f : pp)
    {
        nTris += f.nTriangles();
    }

    DynamicList<labelledTri> dynTris(nTris);
    DynamicList<face> tfaces(8);  // work array

    label facei = 0;
    for (const face& f : pp)
    {
        tfaces.clear();
        f.triangles(points, tfaces);

        for (const auto& t : tfaces)
        {
            dynTris.emplace_back(t[0], t[1], t[2], facei);
        }
        ++facei;
    }

    tris.transfer(dynTris);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triangulatedPatch::update()
{
    triFace_.clear();
    triWght_.clear();

    triangulate(patch_, triFace_);

    const pointField& points = patch_.points();

    const label myProci = UPstream::myProcNo();
    const label numProc = UPstream::nProcs();

    // Calculate the cumulative triangle weights
    triWght_.resize_nocopy(triFace_.size()+1);

    auto iter = triWght_.begin();

    // Set zero value at the start of the tri area/weight list
    scalar patchArea = 0;
    *iter++ = patchArea;

    // Calculate cumulative and total area (processor-local at this point)
    for (const auto& t : triFace_)
    {
        patchArea += t.mag(points);
        *iter++ = patchArea;
    }

    scalarList procSumArea(numProc+1);
    procSumArea[0] = 0;

    {
        scalarList::subList slice(procSumArea, numProc, 1);
        slice[myProci] = patchArea;
        Pstream::allGatherList(slice);
    }

    // Convert to cumulative
    for (label i = 1; i < procSumArea.size(); ++i)
    {
        procSumArea[i] += procSumArea[i-1];
    }

    const scalar offset = procSumArea[myProci];
    const scalar totalArea = procSumArea.back();

    // Apply processor offset and normalise - for a global 0-1 interval
    for (scalar& w : triWght_)
    {
        w = (w + offset) / totalArea;
    }
}


bool Foam::triangulatedPatch::randomLocalPoint
(
    Random& rnd,
    point& result,
    label& facei,
    label& celli
) const
{
    if (triWght_.empty())
    {
        result = point::min;
        facei = -1;
        celli = -1;
        return false;
    }

    const scalar c = rnd.position<scalar>(triWght_.front(), triWght_.back());

    return randomPoint(rnd, c, result, facei, celli);
}


bool Foam::triangulatedPatch::randomGlobalPoint
(
    Random& rnd,
    point& result,
    label& facei,
    label& celli
) const
{
    if (UPstream::parRun())
    {
        const scalar c = rnd.sample01<scalar>();
        const bool ok = randomPoint(rnd, c, result, facei, celli);

        boolList valid(UPstream::listGatherValues(ok));

        // Select the first valid processor
        label proci = valid.find(true);
        Pstream::broadcast(proci);

        return (proci == UPstream::myProcNo());
    }
    else
    {
        return randomLocalPoint(rnd, result, facei, celli);
    }
}


// ************************************************************************* //
