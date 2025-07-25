/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "treeDataPrimitivePatch.H"
#include "indexedOctree.H"
#include "triangle.H"
#include "triSurfaceTools.H"
#include "triFace.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class PatchType>
Foam::treeBoundBoxList
Foam::treeDataPrimitivePatch<PatchType>::boxes(const PatchType& pp)
{
    const auto& points = pp.points();

    treeBoundBoxList bbs(pp.size());

    // Like std::transform with [&](const auto& f)
    // which may not work with older C++ versions

    {
        auto iter = bbs.begin();

        for (const auto& f : pp)
        {
            *iter = treeBoundBox(points, f);
            ++iter;
        }
    }

    return bbs;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PatchType>
inline Foam::treeBoundBox
Foam::treeDataPrimitivePatch<PatchType>::getBounds(const label patchFacei) const
{
    return treeBoundBox(patch_.points(), patch_[patchFacei]);
}


template<class PatchType>
void Foam::treeDataPrimitivePatch<PatchType>::update()
{
    if (cacheBb_)
    {
        bbs_ = treeDataPrimitivePatch<PatchType>::boxes(patch_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PatchType>
Foam::treeDataPrimitivePatch<PatchType>::treeDataPrimitivePatch
(
    const bool cacheBb,
    const PatchType& patch,
    const scalar planarTol
)
:
    patch_(patch),
    cacheBb_(cacheBb),
    planarTol_(planarTol)
{
    update();
}


template<class PatchType>
Foam::treeDataPrimitivePatch<PatchType>::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree
)
:
    tree_(tree)
{}


template<class PatchType>
Foam::treeDataPrimitivePatch<PatchType>::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree
)
:
    tree_(tree)
{}


template<class PatchType>
Foam::treeDataPrimitivePatch<PatchType>::findAllIntersectOp::findAllIntersectOp
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
    DynamicList<label>& shapeMask
)
:
    tree_(tree),
    shapeMask_(shapeMask)
{}


template<class PatchType>
Foam::treeDataPrimitivePatch<PatchType>::
findSelfIntersectOp::findSelfIntersectOp
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
    const label edgeID
)
:
    tree_(tree),
    edgeID_(edgeID)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchType>
Foam::treeBoundBox
Foam::treeDataPrimitivePatch<PatchType>::bounds
(
    const labelUList& indices
) const
{
    treeBoundBox bb;

    for (const label patchFacei : indices)
    {
        bb.add(patch_.points(), patch_[patchFacei]);
    }

    return bb;
}


template<class PatchType>
Foam::volumeType Foam::treeDataPrimitivePatch<PatchType>::getVolumeType
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& oc,
    const point& sample
) const
{
    // Need to determine whether sample is 'inside' or 'outside'
    // Done by finding nearest face. This gives back a face which is
    // guaranteed to contain nearest point. This point can be
    // - in interior of face: compare to face normal
    // - on edge of face: compare to edge normal
    // - on point of face: compare to point normal
    // Unfortunately the octree does not give us back the intersection point
    // or where on the face it has hit so we have to recreate all that
    // information.


    // Find nearest face to sample
    pointIndexHit info = oc.findNearest(sample, sqr(GREAT));

    if (info.index() == -1)
    {
        FatalErrorInFunction
            << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }

    // Get actual intersection point on face
    label facei = info.index();

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << facei;
    }

    const auto& localF = patch_.localFaces()[facei];
    const auto& f = patch_[facei];
    const auto& points = patch_.points();

    // Retest to classify where on face info is. Note: could be improved. We
    // already have point.

    pointHit curHit = f.nearestPoint(sample, points);
    const vector area = f.areaNormal(points);
    const point& curPt = curHit.point();

    //
    // 1] Check whether sample is above face
    //

    if (curHit.hit())
    {
        // Nearest point inside face. Compare to face normal.

        if (debug & 2)
        {
            Pout<< " -> face hit:" << curPt
                << " comparing to face normal " << area << endl;
        }
        return indexedOctree<treeDataPrimitivePatch>::getSide
        (
            area,
            sample - curPt
        );
    }

    if (debug & 2)
    {
        Pout<< " -> face miss:" << curPt;
    }

    //
    // 2] Check whether intersection is on one of the face vertices or
    //    face centre
    //

    const scalar typDimSqr = mag(area) + VSMALL;


    forAll(f, fp)
    {
        const scalar relDistSqr = (magSqr(points[f[fp]] - curPt)/typDimSqr);

        if (relDistSqr < planarTol_)
        {
            // Face intersection point equals face vertex fp

            // Calculate point normal (wrong: uses face normals instead of
            // triangle normals)

            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                patch_.pointNormals()[localF[fp]],
                sample - curPt
            );
        }
    }

    const point fc(f.centre(points));

    const scalar relDistSqr = (magSqr(fc - curPt)/typDimSqr);

    if (relDistSqr < planarTol_)
    {
        // Face intersection point equals face centre. Normal at face centre
        // is already average of face normals

        if (debug & 2)
        {
            Pout<< " -> centre hit:" << fc
                << " distance:" << relDistSqr << endl;
        }

        return indexedOctree<treeDataPrimitivePatch>::getSide
        (
            area,
            sample - curPt
        );
    }


    //
    // 3] Get the 'real' edge the face intersection is on
    //

    for (const label edgei : patch_.faceEdges()[facei])
    {
        pointHit edgeHit =
            patch_.meshEdge(edgei).line(points).nearestDist(sample);

        const scalar relDistSqr = (edgeHit.point().distSqr(curPt)/typDimSqr);

        if (relDistSqr < planarTol_)
        {
            // Face intersection point lies on edge e

            // Calculate edge normal (wrong: uses face normals instead of
            // triangle normals)
            vector edgeNormal(Zero);

            for (const label eFacei : patch_.edgeFaces()[edgei])
            {
                edgeNormal += patch_.faceNormals()[eFacei];
            }

            if (debug & 2)
            {
                Pout<< " -> real edge hit point:" << edgeHit.point()
                    << " comparing to edge normal:" << edgeNormal
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                edgeNormal,
                sample - curPt
            );
        }
    }


    //
    // 4] Get the internal edge the face intersection is on
    //

    forAll(f, fp)
    {
        pointHit edgeHit = linePointRef(points[f[fp]], fc).nearestDist(sample);

        const scalar relDistSqr = (edgeHit.point().distSqr(curPt)/typDimSqr);

        if (relDistSqr < planarTol_)
        {
            // Face intersection point lies on edge between two face triangles

            // Calculate edge normal as average of the two triangle normals
            vector e = points[f[fp]] - fc;
            vector ePrev = points[f[f.rcIndex(fp)]] - fc;
            vector eNext = points[f[f.fcIndex(fp)]] - fc;

            vector nLeft = normalised(ePrev ^ e);
            vector nRight = normalised(e ^ eNext);

            if (debug & 2)
            {
                Pout<< " -> internal edge hit point:" << edgeHit.point()
                    << " comparing to edge normal "
                    << 0.5*(nLeft + nRight)
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                0.5*(nLeft + nRight),
                sample - curPt
            );
        }
    }

    if (debug & 2)
    {
        Pout<< "Did not find sample " << sample
            << " anywhere related to nearest face " << facei << endl
            << "Face:";

        forAll(f, fp)
        {
            Pout<< "    vertex:" << f[fp]
                << "  coord:" << points[f[fp]]
                << endl;
        }
    }

    // Can't determine status of sample with respect to nearest face.
    // Either
    // - tolerances are wrong. (if e.g. face has zero area)
    // - or (more likely) surface is not closed.

    return volumeType::UNKNOWN;
}


// Check if any point on shape is inside searchBox.
template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::overlaps
(
    const label index,
    const treeBoundBox& searchBox
) const
{
    // 1. Quick rejection: bb does not intersect face bb at all
    if
    (
        cacheBb_
      ? !searchBox.overlaps(bbs_[index])
      : !searchBox.overlaps(getBounds(index))
    )
    {
        return false;
    }


    // 2. Check if one or more face points inside

    const auto& points = patch_.points();
    const auto& f = patch_[index];

    if (f.size() == 3)
    {
        const triPointRef tri(points[f[0]], points[f[1]], points[f[2]]);

        return searchBox.intersects(tri);
    }

    if (searchBox.containsAny(points, f))
    {
        return true;
    }

    // 3. Difficult case: all points are outside but connecting edges might
    // go through cube. Use triangle-bounding box intersection.

    const point fc = f.centre(points);

    forAll(f, fp)
    {
        const triPointRef tri
        (
            points[f.thisLabel(fp)], points[f.nextLabel(fp)], fc
        );

        if (searchBox.intersects(tri))
        {
            return true;
        }
    }

    return false;
}


// Check if any point on shape is inside sphere.
template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    // 1. Quick rejection: sphere does not intersect face bb at all
    if
    (
        cacheBb_
      ? !bbs_[index].overlaps(centre, radiusSqr)
      : !getBounds(index).overlaps(centre, radiusSqr)
    )
    {
        return false;
    }

    const auto& points = patch_.points();
    const auto& f = patch_[index];

    pointHit nearHit = f.nearestPoint(centre, points);

    // If the distance to the nearest point on the face from the sphere centres
    // is within the radius, then the sphere touches the face.
    if (sqr(nearHit.distance()) < radiusSqr)
    {
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Searching * * * * * * * * * * * * * * * * //

template<class PatchType>
void Foam::treeDataPrimitivePatch<PatchType>::findNearest
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    const auto& points = patch_.points();

    for (const label index : indices)
    {
        const auto& f = patch_[index];

        const pointHit nearHit = f.nearestPoint(sample, points);
        const scalar distSqr = sqr(nearHit.distance());

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = nearHit.point();
        }
    }
}


template<class PatchType>
void Foam::treeDataPrimitivePatch<PatchType>::findNearestOp::operator()
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


template<class PatchType>
void Foam::treeDataPrimitivePatch<PatchType>::findNearestOp::operator()
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    NotImplemented;
}


template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    return findIntersection(tree_, index, start, end, intersectionPoint);
}


template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::findAllIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    if (shapeMask_.found(index))
    {
        return false;
    }

    return findIntersection(tree_, index, start, end, intersectionPoint);
}


template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::findSelfIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    if (edgeID_ == -1)
    {
        FatalErrorInFunction
            << "EdgeID not set. Please set edgeID to the index of"
            << " the edge you are testing"
            << exit(FatalError);
    }

    const treeDataPrimitivePatch<PatchType>& shape = tree_.shapes();
    const PatchType& patch = shape.patch();

    const auto& f = patch.localFaces()[index];
    const edge& e = patch.edges()[edgeID_];

    if (!f.found(e[0]) && !f.found(e[1]))
    {
        return findIntersection(tree_, index, start, end, intersectionPoint);
    }

    return false;
}


template<class PatchType>
bool Foam::treeDataPrimitivePatch<PatchType>::findIntersection
(
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
)
{
    const treeDataPrimitivePatch<PatchType>& shape = tree.shapes();
    const PatchType& patch = shape.patch();

    const auto& points = patch.points();
    const auto& f = patch[index];

    // Do quick rejection test
    if (shape.cacheBb_)
    {
        const treeBoundBox& faceBb = shape.bbs_[index];

        if ((faceBb.posBits(start) & faceBb.posBits(end)) != 0)
        {
            // start and end in same block outside of faceBb.
            return false;
        }
    }

    const vector dir(end - start);
    pointHit inter;

    if (f.size() == 3)
    {
        inter = triPointRef
        (
            points[f[0]],
            points[f[1]],
            points[f[2]]
        ).intersection(start, dir, intersection::HALF_RAY, shape.planarTol_);
    }
    else
    {
        const pointField& faceCentres = patch.faceCentres();

        inter = f.intersection
        (
            start,
            dir,
            faceCentres[index],
            points,
            intersection::HALF_RAY,
            shape.planarTol_
        );
    }

    if (inter.hit() && inter.distance() <= 1)
    {
        // Note: no extra test on whether intersection is in front of us
        // since using half_ray
        intersectionPoint = inter.point();

        return true;
    }

    return false;
}


// ************************************************************************* //
