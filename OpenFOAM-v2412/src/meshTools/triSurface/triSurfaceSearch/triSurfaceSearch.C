/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "PatchTools.H"
#include "volumeType.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurfaceSearch::checkUniqueHit
(
    const pointIndexHit& currHit,
    const UList<pointIndexHit>& hits,
    const vector& lineVec
) const
{
    // Classify the hit
    label nearType = -1;
    label nearLabel = -1;

    const labelledTri& f = surface()[currHit.index()];

    f.nearestPointClassify
    (
        currHit.hitPoint(),
        surface().points(),
        nearType,
        nearLabel
    );

    if (nearType == triPointRef::POINT)
    {
        // near point

        const label nearPointi = f[nearLabel];

        const labelList& pointFaces =
            surface().pointFaces()[surface().meshPointMap()[nearPointi]];

        forAll(pointFaces, pI)
        {
            const label pointFacei = pointFaces[pI];

            if (pointFacei != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == pointFacei)
                    {
                        return false;
                    }
                }
            }
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        // near edge
        // check if the other face of the edge is already hit

        const labelList& fEdges = surface().faceEdges()[currHit.index()];

        const label edgeI = fEdges[nearLabel];

        const labelList& edgeFaces = surface().edgeFaces()[edgeI];

        forAll(edgeFaces, fI)
        {
            const label edgeFacei = edgeFaces[fI];

            if (edgeFacei != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == edgeFacei)
                    {
                        // Check normals
                        const vector currHitNormal =
                            surface().faceNormals()[currHit.index()];

                        const vector existingHitNormal =
                            surface().faceNormals()[edgeFacei];

                        const label signCurrHit =
                            pos0(currHitNormal & lineVec);

                        const label signExistingHit =
                            pos0(existingHitNormal & lineVec);

                        if (signCurrHit == signExistingHit)
                        {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(nullptr)
{}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const dictionary& dict
)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(nullptr)
{
    // Have optional non-standard search tolerance for gappy surfaces.
    if (dict.readIfPresent("tolerance", tolerance_) && tolerance_ > 0)
    {
        Info<< "    using intersection tolerance " << tolerance_ << endl;
    }

    // Have optional non-standard tree-depth to limit storage.
    if (dict.readIfPresent("maxTreeDepth", maxTreeDepth_) && maxTreeDepth_ > 0)
    {
        Info<< "    using maximum tree depth " << maxTreeDepth_ << endl;
    }
}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const scalar tolerance,
    const label maxTreeDepth
)
:
    surface_(surface),
    tolerance_(tolerance),
    maxTreeDepth_(maxTreeDepth),
    treePtr_(nullptr)
{
    if (tolerance_ < 0)
    {
        tolerance_ = indexedOctree<treeDataTriSurface>::perturbTol();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::~triSurfaceSearch()
{
    clearOut();
}


void Foam::triSurfaceSearch::clearOut()
{
    treePtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::triSurfaceSearch::tree() const
{
    if (!treePtr_)
    {
        // Calculate bb without constructing local point numbering.
        treeBoundBox bb(point::zero);

        if (surface().size())
        {
            label nPoints;
            PatchTools::calcBounds(surface(), bb, nPoints);

            if (nPoints != surface().points().size())
            {
                WarningInFunction
                    << "Surface does not have compact point numbering."
                    << " Of " << surface().points().size()
                    << " only " << nPoints
                    << " are used. This might give problems in some routines."
                    << endl;
            }

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb.inflate(rndGen, 1e-4, ROOTVSMALL);
        }

        const scalar oldTol =
            indexedOctree<treeDataTriSurface>::perturbTol(tolerance_);

        treePtr_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(surface_, tolerance_),
                bb,
                maxTreeDepth_,  // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );

        indexedOctree<treeDataTriSurface>::perturbTol(oldTol);
    }

    return *treePtr_;
}


// Determine inside/outside for samples
Foam::boolList Foam::triSurfaceSearch::calcInside
(
    const pointField& samples
) const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().bb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if (tree().getVolumeType(sample) == volumeType::INSIDE)
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


void Foam::triSurfaceSearch::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const scalar oldTol =
        indexedOctree<treeDataTriSurface>::perturbTol(tolerance());

    const indexedOctree<treeDataTriSurface>& octree = tree();

    const treeDataTriSurface::findNearestOp fOp(octree);

    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = octree.findNearest
        (
            samples[i],
            nearestDistSqr[i],
            fOp
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol(oldTol);
}


Foam::pointIndexHit Foam::triSurfaceSearch::nearest
(
    const point& pt,
    const vector& span
)
const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    return tree().findNearest(pt, nearestDistSqr);
}


void Foam::triSurfaceSearch::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    const scalar oldTol =
        indexedOctree<treeDataTriSurface>::perturbTol(tolerance());

    forAll(start, i)
    {
        info[i] = octree.findLine(start[i], end[i]);
    }

    indexedOctree<treeDataTriSurface>::perturbTol(oldTol);
}


void Foam::triSurfaceSearch::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    const scalar oldTol =
        indexedOctree<treeDataTriSurface>::perturbTol(tolerance());

    forAll(start, i)
    {
        info[i] = octree.findLineAny(start[i], end[i]);
    }

    indexedOctree<treeDataTriSurface>::perturbTol(oldTol);
}


void Foam::triSurfaceSearch::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    const scalar oldTol =
        indexedOctree<treeDataTriSurface>::perturbTol(tolerance());

    // Work array
    DynamicList<pointIndexHit> hits;

    DynamicList<label> shapeMask;

    treeDataTriSurface::findAllIntersectOp allIntersectOp(octree, shapeMask);

    forAll(start, pointi)
    {
        hits.clear();
        shapeMask.clear();

        while (true)
        {
            // See if any intersection between pt and end
            pointIndexHit inter = octree.findLine
            (
                start[pointi],
                end[pointi],
                allIntersectOp
            );

            if (inter.hit())
            {
                const vector lineVec = normalised(end[pointi] - start[pointi]);

                if
                (
                    checkUniqueHit
                    (
                        inter,
                        hits,
                        lineVec
                    )
                )
                {
                    hits.append(inter);
                }

                shapeMask.append(inter.index());
            }
            else
            {
                break;
            }
        }

        info[pointi].transfer(hits);
    }

    indexedOctree<treeDataTriSurface>::perturbTol(oldTol);
}


// ************************************************************************* //
