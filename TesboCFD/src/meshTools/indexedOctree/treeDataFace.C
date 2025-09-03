/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "treeDataFace.H"
#include "polyMesh.H"
#include "triangle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(treeDataFace, 0);

    scalar treeDataFace::tolSqr = sqr(1e-6);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Bound boxes corresponding to specified faces
template<class ElementIds>
static treeBoundBoxList boxesImpl
(
    const primitiveMesh& mesh,
    const ElementIds& elemIds
)
{
    treeBoundBoxList bbs(elemIds.size());

    std::transform
    (
        elemIds.cbegin(),
        elemIds.cend(),
        bbs.begin(),
        [&](label facei)
        {
            return treeBoundBox(mesh.points(), mesh.faces()[facei]);
        }
    );

    return bbs;
}


// Overall bound box for specified faces
template<class ElementIds>
static treeBoundBox boundsImpl
(
    const primitiveMesh& mesh,
    const ElementIds& elemIds
)
{
    treeBoundBox bb;

    for (const label facei : elemIds)
    {
        bb.add(mesh.points(), mesh.faces()[facei]);
    }

    return bb;
}

}  // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::treeBoundBoxList
Foam::treeDataFace::boxes(const primitiveMesh& mesh)
{
    // All faces
    return boxesImpl(mesh, labelRange(mesh.nFaces()));
}


Foam::treeBoundBoxList
Foam::treeDataFace::boxes
(
    const primitiveMesh& mesh,
    const labelRange& range
)
{
    return boxesImpl(mesh, range);
}


Foam::treeBoundBoxList
Foam::treeDataFace::boxes
(
    const primitiveMesh& mesh,
    const labelUList& faceIds
)
{
    return boxesImpl(mesh, faceIds);
}


Foam::treeBoundBox
Foam::treeDataFace::bounds
(
    const primitiveMesh& mesh,
    const labelRange& range
)
{
    return boundsImpl(mesh, range);
}


Foam::treeBoundBox
Foam::treeDataFace::bounds
(
    const primitiveMesh& mesh,
    const labelUList& faceIds
)
{
    return boundsImpl(mesh, faceIds);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool
Foam::treeDataFace::usesFace(const label facei) const
{
    return (!useSubset_ || isTreeFace_.test(facei));
}


inline Foam::treeBoundBox
Foam::treeDataFace::getBounds(const label index) const
{
    const label facei = objectIndex(index);
    return treeBoundBox(mesh_.points(), mesh_.faces()[facei]);
}


void Foam::treeDataFace::update()
{
    if (cacheBb_)
    {
        if (useSubset_)
        {
            bbs_ = treeDataFace::boxes(mesh_, faceLabels_);
        }
        else
        {
            bbs_ = treeDataFace::boxes(mesh_);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataFace::treeDataFace
(
    const bool cacheBb,
    const primitiveMesh& mesh
)
:
    mesh_(mesh),
    faceLabels_(),
    isTreeFace_(),
    useSubset_(false),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataFace::treeDataFace
(
    const bool cacheBb,
    const primitiveMesh& mesh,
    const labelRange& range
)
:
    mesh_(mesh),
    faceLabels_(identity(range)),
    isTreeFace_(range),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataFace::treeDataFace
(
    const bool cacheBb,
    const polyPatch& patch
)
:
    treeDataFace(cacheBb, patch.boundaryMesh().mesh(), patch.range())
{}


Foam::treeDataFace::treeDataFace
(
    const bool cacheBb,
    const primitiveMesh& mesh,
    const labelUList& faceLabels
)
:
    mesh_(mesh),
    faceLabels_(faceLabels),
    isTreeFace_(mesh_.nFaces(), faceLabels_),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


Foam::treeDataFace::treeDataFace
(
    const bool cacheBb,
    const primitiveMesh& mesh,
    labelList&& faceLabels
)
:
    mesh_(mesh),
    faceLabels_(std::move(faceLabels)),
    isTreeFace_(mesh_.nFaces(), faceLabels_),
    useSubset_(true),
    cacheBb_(cacheBb)
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox
Foam::treeDataFace::bounds(const labelUList& indices) const
{
    if (useSubset_)
    {
        treeBoundBox bb;

        for (const label index : indices)
        {
            const label facei = faceLabels_[index];

            bb.add(mesh_.points(), mesh_.faces()[facei]);
        }

        return bb;
    }

    return treeDataFace::bounds(mesh_, indices);
}


Foam::tmp<Foam::pointField> Foam::treeDataFace::centres() const
{
    if (useSubset_)
    {
        return tmp<pointField>::New(mesh_.faceCentres(), faceLabels_);
    }

    return mesh_.faceCentres();
}


Foam::volumeType Foam::treeDataFace::getVolumeType
(
    const indexedOctree<treeDataFace>& oc,
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

    const label index = info.index();

    if (index == -1)
    {
        FatalErrorInFunction
            << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }


    // Get actual intersection point on face
    const label facei = objectIndex(index);

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << facei;
    }

    const pointField& points = mesh_.points();

    // Retest to classify where on face info is. Note: could be improved. We
    // already have point.

    const face& f = mesh_.faces()[facei];
    const vector& area = mesh_.faceAreas()[facei];
    const point& fc = mesh_.faceCentres()[facei];

    pointHit curHit = f.nearestPoint(sample, points);
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
        return indexedOctree<treeDataFace>::getSide(area, sample - curPt);
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

    for (const label fp : f)
    {
        const scalar relDistSqr = (magSqr(points[fp] - curPt)/typDimSqr);

        if (relDistSqr < tolSqr)
        {
            // Face intersection point equals face vertex

            // Calculate point normal (wrong: uses face normals instead of
            // triangle normals)

            vector pointNormal(Zero);

            for (const label ptFacei : mesh_.pointFaces()[fp])
            {
                if (usesFace(ptFacei))
                {
                    pointNormal += normalised(mesh_.faceAreas()[ptFacei]);
                }
            }

            if (debug & 2)
            {
                    Pout<< " -> face point hit :" << points[fp]
                        << " point normal:" << pointNormal
                        << " distance:" << relDistSqr << endl;
            }
            return indexedOctree<treeDataFace>::getSide
            (
                pointNormal,
                sample - curPt
            );
        }
    }

    const scalar relDistSqr = (magSqr(fc - curPt)/typDimSqr);

    if (relDistSqr < tolSqr)
    {
        // Face intersection point equals face centre. Normal at face centre
        // is already average of face normals

        if (debug & 2)
        {
            Pout<< " -> centre hit:" << fc
                << " distance:" << relDistSqr << endl;
        }

        return indexedOctree<treeDataFace>::getSide(area,  sample - curPt);
    }


    //
    // 3] Get the 'real' edge the face intersection is on
    //

    for (const label edgei : mesh_.faceEdges()[facei])
    {
        pointHit edgeHit =
            mesh_.edges()[edgei].line(points).nearestDist(sample);

        const scalar relDistSqr = edgeHit.point().distSqr(curPt)/typDimSqr;

        if (relDistSqr < tolSqr)
        {
            // Face intersection point lies on edge e

            // Calculate edge normal (wrong: uses face normals instead of
            // triangle normals)

            vector edgeNormal(Zero);

            for (const label eFacei : mesh_.edgeFaces()[edgei])
            {
                if (usesFace(eFacei))
                {
                    edgeNormal += normalised(mesh_.faceAreas()[eFacei]);
                }
            }

            if (debug & 2)
            {
                Pout<< " -> real edge hit point:" << edgeHit.point()
                    << " comparing to edge normal:" << edgeNormal
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return indexedOctree<treeDataFace>::getSide
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

        const scalar relDistSqr = edgeHit.point().distSqr(curPt)/typDimSqr;

        if (relDistSqr < tolSqr)
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
            return indexedOctree<treeDataFace>::getSide
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
            Pout<< "    vertex:" << f[fp] << "  coord:" << points[f[fp]]
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
bool Foam::treeDataFace::overlaps
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

    const pointField& points = mesh_.points();

    // 2. Check if one or more face points inside
    const label facei = objectIndex(index);

    const face& f = mesh_.faces()[facei];

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

    const point& fc = mesh_.faceCentres()[facei];

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
bool Foam::treeDataFace::overlaps
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

    const label facei = objectIndex(index);

    const face& f = mesh().faces()[facei];

    pointHit nearHit = f.nearestPoint(centre, mesh().points());

    if (sqr(nearHit.distance()) < radiusSqr)
    {
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Searching * * * * * * * * * * * * * * * * //

Foam::treeDataFace::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataFace>& tree
)
:
    tree_(tree)
{}


Foam::treeDataFace::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataFace>& tree
)
:
    tree_(tree)
{}


void Foam::treeDataFace::findNearest
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
        const label facei = objectIndex(index);

        const face& f = mesh().faces()[facei];

        pointHit nearHit = f.nearestPoint(sample, mesh().points());
        scalar distSqr = sqr(nearHit.distance());

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = nearHit.point();
        }
    }
}


void Foam::treeDataFace::findNearestOp::operator()
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


void Foam::treeDataFace::findNearestOp::operator()
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


bool Foam::treeDataFace::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    const treeDataFace& shape = tree_.shapes();

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

    const label facei = shape.objectIndex(index);

    const vector dir(end - start);

    pointHit inter = shape.mesh_.faces()[facei].intersection
    (
        start,
        dir,
        shape.mesh_.faceCentres()[facei],
        shape.mesh_.points(),
        intersection::HALF_RAY
    );

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
