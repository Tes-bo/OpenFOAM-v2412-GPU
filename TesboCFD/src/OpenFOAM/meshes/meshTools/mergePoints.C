/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "ListOps.H"  // sortedOrder

// * * * * * * * * * * * * * * * Implementation  * * * * * * * * * * * * * * //

template<class PointList, class IndexerOp>
Foam::label Foam::Detail::mergePoints
(
    const PointList& points,
    const IndexerOp& indexer,
    const label nSubPoints,
    labelList& pointToUnique,
    labelList& uniquePoints,
    const scalar mergeTol,
    const bool verbose
)
{
    const label nTotPoints = points.size();

    if (!nTotPoints || !nSubPoints)
    {
        // Nothing to do
        pointToUnique = identity(nTotPoints);
        uniquePoints = pointToUnique;
        return 0;  // No points removed
    }

    // Properly size for old to new mapping array
    pointToUnique.resize_nocopy(nTotPoints);


    // Use the boundBox minimum as the reference point. This
    // stretches distances with fewer collisions than a mid-point
    // reference would.

    auto comparePoint(points[indexer(0)]);
    for (label pointi = 1; pointi < nSubPoints; ++pointi)
    {
        comparePoint = min(comparePoint, points[indexer(pointi)]);
    }

    // We're comparing distance squared to reference point first.
    // Say if starting from two close points:
    //     x, y, z
    //     x+mergeTol, y+mergeTol, z+mergeTol
    // Then the magSqr of both will be
    //     x^2+y^2+z^2
    //     x^2+y^2+z^2 + 2*mergeTol*(x+z+y) + mergeTol^2*...
    // so the difference will be 2*mergeTol*(x+y+z)

    const scalar mergeTolSqr(magSqr(mergeTol));

    // Use magSqr distance for the points sort order
    List<scalar> sqrDist(nSubPoints);
    for (label pointi = 0; pointi < nSubPoints; ++pointi)
    {
        const auto& p = points[indexer(pointi)];

        sqrDist[pointi] =
        (
            // Use scalar precision
            magSqr(scalar(p.x() - comparePoint.x()))
          + magSqr(scalar(p.y() - comparePoint.y()))
          + magSqr(scalar(p.z() - comparePoint.z()))
        );
    }
    labelList order(Foam::sortedOrder(sqrDist));

    List<scalar> sortedTol(nSubPoints);
    forAll(order, sorti)
    {
        const auto& p = points[indexer(order[sorti])];

        sortedTol[sorti] =
        (
            2*mergeTol*
            (
                // Use scalar precision
                mag(scalar(p.x() - comparePoint.x()))
              + mag(scalar(p.y() - comparePoint.y()))
              + mag(scalar(p.z() - comparePoint.z()))
            )
        );
    }


    // Bookkeeping parameters
    // ~~~~~~~~~~~~~~~~~~~~~~

    // Will only be working on a subset of the points
    // Can use a slice of pointToUnique (full length not needed until later).
    SubList<label> subPointMap(pointToUnique, nSubPoints);

    // Track number of unique points - this will form an offsets table
    labelList newPointCounts(nSubPoints, Zero);

    label nNewPoints = 0;
    for (label sorti = 0; sorti < order.size(); ++sorti)
    {
        // The (sub)point index
        const label pointi = order[sorti];
        const scalar currDist = sqrDist[order[sorti]];
        const auto& currPoint = points[indexer(pointi)];

        // Compare to previous points to find equal one
        // - automatically a no-op for sorti == 0 (the first point)

        bool matched = false;

        for
        (
            label prevSorti = sorti - 1;
            (
                prevSorti >= 0
             && (mag(sqrDist[order[prevSorti]] - currDist) <= sortedTol[sorti])
            );
            --prevSorti
        )
        {
            const label prevPointi = order[prevSorti];
            const auto& prevPoint = points[indexer(prevPointi)];

            // Matched within tolerance?
            matched =
            (
                (
                    // Use scalar precision
                    magSqr(scalar(currPoint.x() - prevPoint.x()))
                  + magSqr(scalar(currPoint.y() - prevPoint.y()))
                  + magSqr(scalar(currPoint.z() - prevPoint.z()))
                ) <= mergeTolSqr
            );

            if (matched)
            {
                // Both pointi and prevPointi have similar coordinates.
                // Map to the same new point.
                subPointMap[pointi] = subPointMap[prevPointi];

                if (verbose)
                {
                    Pout<< "Foam::mergePoints : [" << subPointMap[pointi]
                        << "] Point " << pointi << " duplicate of "
                        << prevPointi << " : coordinates:" << currPoint
                        << " and " << prevPoint << endl;
                }
                break;
            }
        }

        if (!matched)
        {
            // Differs. Store new point.
            subPointMap[pointi] = nNewPoints++;

            /// Too verbose
            ///if (verbose)
            ///{
            ///    Pout<< "Foam::mergePoints : [" << subPointMap[pointi]
            ///        << "] Point " << pointi << endl;
            ///}
        }
        ++newPointCounts[subPointMap[pointi]];
    }

    const label nDupPoints(nSubPoints - nNewPoints);
    const label nUniqPoints(nTotPoints - nDupPoints);

    if (verbose)
    {
        Pout<< "Foam::mergePoints : "
            << "Merging removed " << nDupPoints << '/'
            << nTotPoints << " points" << endl;
    }

    if (!nDupPoints)
    {
        // Nothing to do
        pointToUnique = identity(nTotPoints);
        uniquePoints = pointToUnique;
        return 0;  // No points removed
    }


    // The subPointMap now contains a mapping of the sub-selection
    // to the list of (sorted) merged points.
    // Get its sort order to bundle according to the merged point target.
    // This is in effect an adjacent list of graph edges to mapping back
    // to the merged points, but in compact form.
    // Use the previously obtained newPointCounts for the offsets list.

    labelList lookupMerged(std::move(order));
    Foam::sortedOrder(subPointMap, lookupMerged);

    // Remap inplace to the original points ids
    for (label& idx : lookupMerged)
    {
        idx = indexer(idx);
    }
    // The subPointMap slice is not needed beyond here


    // Setup initial identity +1 mapping for pointToUnique
    // The +1 allows negatives to mark duplicates

    Foam::identity(pointToUnique, 1);

    // The newPointCounts is an offsets table that we use to walk
    // across the adjacency list (lookupMerged), picking the original
    // point with the lowest id as the one to retain (master).
    {
        label beg = 0;
        for (const label len : newPointCounts)
        {
            if (!len) continue;  // Can be empty

            const label end = (beg + len);

            // Pass 1:
            // Find the 'master' (lowest point id)

            label masterPointi = lookupMerged[beg];

            for (label iter = beg + 1; iter < end; ++iter)
            {
                const label origPointi = lookupMerged[iter];

                if (masterPointi > origPointi)
                {
                    masterPointi = origPointi;
                }
            }

            // Pass 2:
            // Markup duplicate points, encoding information about master
            for (label iter = beg; iter < end; ++iter)
            {
                const label origPointi = lookupMerged[iter];

                if (masterPointi != origPointi)
                {
                    // Encode the originating 'master' point
                    pointToUnique[origPointi] = (-masterPointi-1);
                }
            }

            beg = end;
        }
    }

    // Now have all the information needed

    uniquePoints.resize_nocopy(nUniqPoints);
    {
        label uniquei = 0;

        forAll(pointToUnique, pointi)
        {
            const label origPointi = pointToUnique[pointi];

            if (origPointi > 0)
            {
                // Subtract one to align addressing
                uniquePoints[uniquei] = (origPointi - 1);
                pointToUnique[pointi] = uniquei;
                ++uniquei;
            }
            else
            {
                // A duplicate point. Also guaranteed that the 'master' point
                // has a lower index and thus already been seen.
                const label masterPointi = mag(origPointi) - 1;
                pointToUnique[pointi] = pointToUnique[masterPointi];
            }
        }
    }

    return nDupPoints;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PointList>
Foam::label Foam::mergePoints
(
    const PointList& points,
    labelList& pointToUnique,
    labelList& uniquePoints,
    const scalar mergeTol,
    const bool verbose
)
{
    const label nTotPoints = points.size();

    if (!nTotPoints)
    {
        // Nothing to do
        pointToUnique.clear();
        uniquePoints.clear();
        return 0;  // No points removed
    }

    return Foam::Detail::mergePoints
    (
        points,
        identityOp(),   // identity indexer
        nTotPoints,     // == nSubPoints
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );
}


template<class PointList>
Foam::label Foam::mergePoints
(
    const PointList& points,
    const labelUList& selection,
    labelList& pointToUnique,
    labelList& uniquePoints,
    const scalar mergeTol,
    const bool verbose
)
{
    const label nTotPoints = points.size();
    const label nSubPoints = selection.size();

    if (!nTotPoints || !nSubPoints)
    {
        // Nothing to do
        pointToUnique.clear();
        uniquePoints.clear();
        return 0;  // No points removed
    }

    const auto indexer = [&](const label i) -> label { return selection[i]; };

    return Foam::Detail::mergePoints
    (
        points,
        indexer,
        nSubPoints,
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );
}


template<class PointList>
Foam::label Foam::mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
)
{
    labelList uniquePoints;
    const label nChanged = Foam::mergePoints
    (
        points,
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );

    // Number of unique points
    return (points.size() - nChanged);
}


template<class PointList>
Foam::label Foam::inplaceMergePoints
(
    PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
)
{
    labelList uniquePoints;
    const label nChanged = Foam::mergePoints
    (
        points,
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );

    if (nChanged)
    {
        // Overwrite
        points = List<typename PointList::value_type>(points, uniquePoints);
    }
    else
    {
        // TDB:
        // pointToUnique.clear();
    }

    return nChanged;
}


template<class PointList>
Foam::label Foam::inplaceMergePoints
(
    PointList& points,
    const labelUList& selection,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
)
{
    labelList uniquePoints;
    const label nChanged = Foam::mergePoints
    (
        points,
        selection,
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );

    if (nChanged)
    {
        // Overwrite
        points = List<typename PointList::value_type>(points, uniquePoints);
    }
    else
    {
        // TDB:
        // pointToUnique.clear();
    }

    return nChanged;
}


template<class PointList>
Foam::label Foam::mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique,
    List<typename PointList::value_type>& newPoints
)
{
    const label nTotPoints = points.size();

    if (!nTotPoints)
    {
        // Nothing to do
        pointToUnique.clear();
        newPoints.clear();
        return 0;  // No points removed
    }

    labelList uniquePoints;
    const label nChanged = Foam::mergePoints
    (
        points,
        pointToUnique,
        uniquePoints,
        mergeTol,
        verbose
    );

    if (nChanged)
    {
        newPoints = List<typename PointList::value_type>(points, uniquePoints);
    }
    else
    {
        // TDB:
        // pointToUnique.clear();
        newPoints = points;
    }

    return nChanged;
}


// ************************************************************************* //
