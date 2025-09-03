/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "treeDataEdge.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(treeDataEdge);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Bound boxes corresponding to specified edges
template<class ElementIds>
static treeBoundBoxList boxesImpl
(
    const edgeList& edges,
    const pointField& points,
    const ElementIds& elemIds
)
{
    treeBoundBoxList bbs(elemIds.size());

    std::transform
    (
        elemIds.cbegin(),
        elemIds.cend(),
        bbs.begin(),
        [&](label edgei)
        {
            return treeBoundBox(edges[edgei].box(points));
        }
    );

    return bbs;
}


// Overall bound box for specified edges
template<class ElementIds>
static treeBoundBox boundsImpl
(
    const edgeList& edges,
    const pointField& points,
    const ElementIds& elemIds
)
{
    treeBoundBox bb;

    for (const label edgei : elemIds)
    {
        const edge& e = edges[edgei];

        bb.add(points[e.first()], points[e.second()]);
    }

    return bb;
}

}  // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::treeBoundBoxList
Foam::treeDataEdge::boxes
(
    const edgeList& edges,
    const pointField& points
)
{
    treeBoundBoxList bbs(edges.size());

    std::transform
    (
        edges.cbegin(),
        edges.cend(),
        bbs.begin(),
        [&](const edge& e) { return treeBoundBox(e.box(points)); }
    );

    return bbs;
}


Foam::treeBoundBoxList
Foam::treeDataEdge::boxes
(
    const edgeList& edges,
    const pointField& points,
    const labelRange& range
)
{
    return boxesImpl(edges, points, range);
}


Foam::treeBoundBoxList
Foam::treeDataEdge::boxes
(
    const edgeList& edges,
    const pointField& points,
    const labelUList& edgeIds
)
{
    return boxesImpl(edges, points, edgeIds);
}


Foam::treeBoundBox
Foam::treeDataEdge::bounds
(
    const edgeList& edges,
    const pointField& points
)
{
    treeBoundBox bb;

    for (const edge& e : edges)
    {
        bb.add(points[e.first()], points[e.second()]);
    }

    return bb;
}


Foam::treeBoundBox
Foam::treeDataEdge::bounds
(
    const edgeList& edges,
    const pointField& points,
    const labelRange& range
)
{
    return boundsImpl(edges, points, range);
}


Foam::treeBoundBox
Foam::treeDataEdge::bounds
(
    const edgeList& edges,
    const pointField& points,
    const labelUList& edgeIds
)
{
    return boundsImpl(edges, points, edgeIds);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::treeDataEdge::update()
{
    if (cacheBb_)
    {
        if (useSubset_)
        {
            bbs_ = treeDataEdge::boxes(edges_, points_, edgeLabels_);
        }
        else
        {
            bbs_ = treeDataEdge::boxes(edges_, points_);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points
)
:
    points_(points),
    edges_(edges),
    edgeLabels_(),
    useSubset_(false),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points,
    const labelRange& range
)
:
    points_(points),
    edges_(edges),
    edgeLabels_(identity(range)),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points,
    const labelUList& edgeLabels
)
:
    points_(points),
    edges_(edges),
    edgeLabels_(edgeLabels),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataEdge::treeDataEdge
(
    const bool cacheBb,
    const edgeList& edges,
    const pointField& points,
    labelList&& edgeLabels
)
:
    points_(points),
    edges_(edges),
    edgeLabels_(std::move(edgeLabels)),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeDataEdge::bounds(const labelUList& indices) const
{
    if (useSubset_)
    {
        treeBoundBox bb;

        for (const label index : indices)
        {
            const edge& e = edges_[edgeLabels_[index]];

            bb.add(points_[e.first()], points_[e.second()]);
        }

        return bb;
    }

    return treeDataEdge::bounds(edges_, points_, indices);
}


Foam::tmp<Foam::pointField> Foam::treeDataEdge::centres() const
{
    tmp<pointField> tpts;

    if (useSubset_)
    {
        tpts = tmp<pointField>::New(edgeLabels_.size());

        std::transform
        (
            edgeLabels_.cbegin(),
            edgeLabels_.cend(),
            tpts.ref().begin(),
            [&](label edgei) { return edges_[edgei].centre(points_); }
        );
    }
    else
    {
        tpts = tmp<pointField>::New(edges_.size());

        std::transform
        (
            edges_.cbegin(),
            edges_.cend(),
            tpts.ref().begin(),
            [&](const edge& e) { return e.centre(points_); }
        );
    }

    return tpts;
}


Foam::volumeType Foam::treeDataEdge::getVolumeType
(
    const indexedOctree<treeDataEdge>& oc,
    const point& sample
) const
{
    return volumeType::UNKNOWN;
}


bool Foam::treeDataEdge::overlaps
(
    const label index,
    const treeBoundBox& searchBox
) const
{
    point intersect;
    return searchBox.intersects(this->line(index), intersect);
}


bool Foam::treeDataEdge::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    const pointHit nearHit = this->line(index).nearestDist(centre);

    return (sqr(nearHit.distance()) <= radiusSqr);
}


// * * * * * * * * * * * * * * * * Searching * * * * * * * * * * * * * * * * //

Foam::treeDataEdge::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataEdge>& tree
)
:
    tree_(tree)
{}


Foam::treeDataEdge::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataEdge>& tree
)
{}


void Foam::treeDataEdge::findNearest
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
        pointHit nearHit = this->line(index).nearestDist(sample);

        const scalar distSqr = sqr(nearHit.distance());

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = nearHit.point();
        }
    }
}


void Foam::treeDataEdge::findNearestOp::operator()
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


void Foam::treeDataEdge::findNearestOp::operator()
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    const treeDataEdge& shape = tree_.shapes();

    const treeBoundBox lnBb(ln.box());

    // Best so far
    scalar nearestDistSqr = linePoint.distSqr(nearestPoint);

    for (const label index : indices)
    {
        // Note: could do bb test ? Worthwhile?

        // Nearest point on line
        point ePoint, lnPt;
        const scalar dist = shape.line(index).nearestDist(ln, ePoint, lnPt);

        const scalar distSqr = sqr(dist);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            linePoint = lnPt;
            nearestPoint = ePoint;

            tightest = lnBb;
            tightest.grow(dist);
        }
    }
}


bool Foam::treeDataEdge::findIntersectOp::operator()
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
