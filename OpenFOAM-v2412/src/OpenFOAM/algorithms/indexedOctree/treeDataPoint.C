/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "treeDataPoint.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(treeDataPoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataPoint::treeDataPoint(const pointField& points)
:
    points_(points),
    useSubset_(false)
{}


Foam::treeDataPoint::treeDataPoint
(
    const pointField& points,
    const labelUList& pointLabels,
    const bool useSubsetPoints
)
:
    points_(points),
    pointLabels_(pointLabels),
    useSubset_(useSubsetPoints)
{}


Foam::treeDataPoint::treeDataPoint
(
    const pointField& points,
    labelList&& pointLabels,
    const bool useSubsetPoints
)
:
    points_(points),
    pointLabels_(std::move(pointLabels)),
    useSubset_(useSubsetPoints)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeDataPoint::bounds(const labelUList& indices) const
{
    if (useSubset_)
    {
        treeBoundBox bb;
        for (const label index : indices)
        {
            bb.add(points_[pointLabels_[index]]);
        }
        return bb;
    }

    return treeBoundBox(points_, indices);
}


Foam::tmp<Foam::pointField> Foam::treeDataPoint::centres() const
{
    if (useSubset_)
    {
        return tmp<pointField>::New(points_, pointLabels_);
    }

    return points_;
}


Foam::volumeType Foam::treeDataPoint::getVolumeType
(
    const indexedOctree<treeDataPoint>& oc,
    const point& sample
) const
{
    return volumeType::UNKNOWN;
}


bool Foam::treeDataPoint::overlaps
(
    const label index,
    const treeBoundBox& searchBox
) const
{
    return searchBox.contains(centre(index));
}


bool Foam::treeDataPoint::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    return (centre.distSqr(this->centre(index)) <= radiusSqr);
}


// * * * * * * * * * * * * * * * * Searching * * * * * * * * * * * * * * * * //

Foam::treeDataPoint::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataPoint>& tree
)
:
    tree_(tree)
{}


Foam::treeDataPoint::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataPoint>& tree
)
{}


void Foam::treeDataPoint::findNearest
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    for (const label index : indices)
    {
        const point& pt = centre(index);

        const scalar distSqr = sample.distSqr(pt);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = pt;
        }
    }
}


void Foam::treeDataPoint::findNearestOp::operator()
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    tree_.shapes().findNearest
    (
        indices,
        sample,
        nearestDistSqr,
        minIndex,
        nearestPoint
    );
}


void Foam::treeDataPoint::findNearestOp::operator()
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    const treeDataPoint& shape = tree_.shapes();

    const treeBoundBox lnBb(ln.box());

    // Best so far
    scalar nearestDistSqr = GREAT;
    if (minIndex >= 0)
    {
        nearestDistSqr = linePoint.distSqr(nearestPoint);
    }

    for (const label index : indices)
    {
        const point& pt = shape.centre(index);

        if (tightest.contains(pt))
        {
            // Nearest point on line
            pointHit pHit = ln.nearestDist(pt);
            const scalar distSqr = sqr(pHit.distance());

            if (distSqr < nearestDistSqr)
            {
                nearestDistSqr = distSqr;
                minIndex = index;
                linePoint = pHit.point();
                nearestPoint = pt;

                tightest = lnBb;
                tightest.grow(pHit.distance());
            }
        }
    }
}


bool Foam::treeDataPoint::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& result
) const
{
    NotImplemented;
    return false;
}


// ************************************************************************* //
