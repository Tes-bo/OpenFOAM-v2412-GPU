/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "linearDistance.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearDistance, 0);
    addToRunTimeSelectionTable(cellSizeFunction, linearDistance, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearDistance::linearDistance
(
    const dictionary& initialPointsDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize,
    const labelList regionIndices
)
:
    cellSizeFunction
    (
        typeName,
        initialPointsDict,
        surface,
        defaultCellSize,
        regionIndices
    ),
    distanceCellSize_
    (
        coeffsDict().get<scalar>("distanceCellSizeCoeff") * defaultCellSize
    ),
    distance_
    (
        coeffsDict().get<scalar>("distanceCoeff") * defaultCellSize
    ),
    distanceSqr_(sqr(distance_))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::linearDistance::sizeFunction
(
    const point& pt,
    scalar d,
    label index
) const
{
    const scalar interpolatedSize
        = surfaceCellSizeFunction_().interpolate(pt, index);

    const scalar gradient
        = (distanceCellSize_ - interpolatedSize)
          /distance_;

    scalar size = gradient*d + interpolatedSize;

    return size;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::linearDistance::sizeLocations
(
    const pointIndexHit& hitPt,
    const vector& n,
    pointField& shapePts,
    scalarField& shapeSizes
) const
{
    const Foam::point& pt = hitPt.hitPoint();

    if (sideMode_ == rmBothsides)
    {
        shapePts.resize(2);
        shapeSizes.resize(2);

        shapePts[0] = pt - n*distance_;
        shapePts[1] = pt + n*distance_;

        shapeSizes[0] = distanceCellSize_;
        shapeSizes[1] = distanceCellSize_;
    }
    else if (sideMode_ == smInside)
    {
        shapePts.resize(1);
        shapeSizes.resize(1);

        shapePts[0] = pt - n*distance_;
        shapeSizes[0] = distanceCellSize_;
    }
    else if (sideMode_ == smOutside)
    {
        shapePts.resize(1);
        shapeSizes.resize(1);

        shapePts[0] = pt + n*distance_;
        shapeSizes[0] = distanceCellSize_;
    }

    return false;
}


bool Foam::linearDistance::cellSize(const point& pt, scalar& size) const
{
    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, distanceSqr_),
        regionIndices_,
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    if (hitInfo.hit())
    {
        const point& hitPt = hitInfo.point();
        const label hitIndex = hitInfo.index();

        const scalar dist = hitPt.dist(pt);

        if (sideMode_ == rmBothsides)
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            return true;
        }

        // If the nearest point is essentially on the surface, do not do a
        // getVolumeType calculation, as it will be prone to error.
        if (dist < snapToSurfaceTol_)
        {
            size = sizeFunction(hitPt, 0, hitIndex);

            return true;
        }

        pointField ptF(1, pt);
        List<volumeType> vTL;

        surface_.getVolumeType(ptF, vTL);

        bool functionApplied = false;

        if
        (
            sideMode_ == smInside
         && vTL[0] == volumeType::INSIDE
        )
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            functionApplied = true;
        }
        else if
        (
            sideMode_ == smOutside
         && vTL[0] == volumeType::OUTSIDE
        )
        {
            size = sizeFunction(hitPt, dist, hitIndex);

            functionApplied = true;
        }

        return functionApplied;
    }

    return false;
}


bool Foam::linearDistance::setCellSize(const pointField& pts)
{
//    labelHashSet surfaceAlreadyHit(surfaceCellSize_.size());

//    for (const Foam::point& pt : pts)
//    {
//        List<pointIndexHit> hits;

//        surface_.findNearest
//        (
//            pointField(1, pt),
//            scalarField(1, distanceSqr_),
//            regionIndices_,
//            hits
//        );

//        const label surfHitI = hits[0].index();

//        if
//        (
//            hits[0].hit()
//         && !surfaceAlreadyHit.found(surfHitI)
//        )
//        {
//            // Halving cell size is arbitrary
//            surfaceCellSizeFunction_().refineSurfaceSize(surfHitI);

//            surfaceAlreadyHit.insert(surfHitI);
//        }
//    }

//    return true;
    return false;
}


// ************************************************************************* //
