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

#include "treeDataCell.H"
#include "indexedOctree.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(treeDataCell);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Bound boxes corresponding to specified cells
template<class ElementIds>
static treeBoundBoxList boxesImpl
(
    const primitiveMesh& mesh,
    const ElementIds& elemIds
)
{
    treeBoundBoxList bbs(elemIds.size());

    if (mesh.hasCellPoints())
    {
        std::transform
        (
            elemIds.cbegin(),
            elemIds.cend(),
            bbs.begin(),
            [&](label celli)
            {
                return treeBoundBox(mesh.points(), mesh.cellPoints(celli));
            }
        );
    }
    else
    {
        std::transform
        (
            elemIds.cbegin(),
            elemIds.cend(),
            bbs.begin(),
            [&](label celli)
            {
                return treeBoundBox(mesh.cells()[celli].box(mesh));
            }
        );
    }

    return bbs;
}

}  // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::treeBoundBoxList
Foam::treeDataCell::boxes(const primitiveMesh& mesh)
{
    // All cells
    return boxesImpl(mesh, labelRange(mesh.nCells()));
}


Foam::treeBoundBoxList
Foam::treeDataCell::boxes
(
    const primitiveMesh& mesh,
    const labelUList& cellIds
)
{
    return boxesImpl(mesh, cellIds);
}


Foam::treeBoundBox
Foam::treeDataCell::bounds
(
    const primitiveMesh& mesh,
    const labelUList& cellIds
)
{
    treeBoundBox bb;

    if (mesh.hasCellPoints())
    {
        for (const label celli : cellIds)
        {
            bb.add(mesh.points(), mesh.cellPoints(celli));
        }
    }
    else
    {
        for (const label celli : cellIds)
        {
            bb.add(mesh.cells()[celli].box(mesh));
        }
    }

    return bb;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::treeBoundBox
Foam::treeDataCell::getBounds(const label index) const
{
    return treeBoundBox(mesh_.cellBb(objectIndex(index)));
}


void Foam::treeDataCell::update()
{
    if (cacheBb_)
    {
        if (useSubset_)
        {
            bbs_ = treeDataCell::boxes(mesh_, cellLabels_);
        }
        else
        {
            bbs_ = treeDataCell::boxes(mesh_);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(),
    useSubset_(false),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    const labelUList& cellLabels,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(cellLabels),
    useSubset_(true),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


Foam::treeDataCell::treeDataCell
(
    const bool cacheBb,
    const polyMesh& mesh,
    labelList&& cellLabels,
    const polyMesh::cellDecomposition decompMode
)
:
    mesh_(mesh),
    cellLabels_(std::move(cellLabels)),
    useSubset_(true),
    cacheBb_(cacheBb),
    decompMode_(decompMode)
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox
Foam::treeDataCell::bounds(const labelUList& indices) const
{
    if (useSubset_)
    {
        treeBoundBox bb;

        if (mesh_.hasCellPoints())
        {
            for (const label index : indices)
            {
                const label celli = cellLabels_[index];

                bb.add(mesh_.points(), mesh_.cellPoints(celli));
            }
        }
        else
        {
            for (const label index : indices)
            {
                const label celli = cellLabels_[index];

                bb.add(mesh_.cells()[celli].box(mesh_));
            }
        }

        return bb;
    }

    return treeDataCell::bounds(mesh_, indices);
}


Foam::tmp<Foam::pointField> Foam::treeDataCell::centres() const
{
    if (useSubset_)
    {
        return tmp<pointField>::New(mesh_.cellCentres(), cellLabels_);
    }

    return mesh_.cellCentres();
}


bool Foam::treeDataCell::overlaps
(
    const label index,
    const treeBoundBox& searchBox
) const
{
    return
    (
        cacheBb_
      ? searchBox.overlaps(bbs_[index])
      : searchBox.overlaps(getBounds(index))
    );
}


bool Foam::treeDataCell::contains
(
    const label index,
    const point& sample
) const
{
    return mesh_.pointInCell(sample, objectIndex(index), decompMode_);
}


// * * * * * * * * * * * * * * * * Searching * * * * * * * * * * * * * * * * //

Foam::treeDataCell::findNearestOp::findNearestOp
(
    const indexedOctree<treeDataCell>& tree
)
:
    tree_(tree)
{}


Foam::treeDataCell::findIntersectOp::findIntersectOp
(
    const indexedOctree<treeDataCell>& tree
)
:
    tree_(tree)
{}


void Foam::treeDataCell::findNearest
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


void Foam::treeDataCell::findNearestOp::operator()
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


void Foam::treeDataCell::findNearestOp::operator()
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


bool Foam::treeDataCell::findIntersectOp::operator()
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    const treeDataCell& shape = tree_.shapes();

    // Do quick rejection test
    if (shape.cacheBb_)
    {
        const treeBoundBox& cellBb = shape.bbs_[index];

        if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
        {
            // start and end in same block outside of cellBb.
            return false;
        }
    }
    else
    {
        const treeBoundBox cellBb = shape.getBounds(index);

        if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
        {
            // start and end in same block outside of cellBb.
            return false;
        }
    }


    // Do intersection with all faces of cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Disable picking up intersections behind us.
    const scalar oldTol = intersection::setPlanarTol(0);

    const vector dir(end - start);
    scalar minDistSqr = magSqr(dir);
    bool hasMin = false;

    const label celli = shape.objectIndex(index);

    for (const label facei : shape.mesh().cells()[celli])
    {
        const face& f = shape.mesh().faces()[facei];

        pointHit inter = f.ray
        (
            start,
            dir,
            shape.mesh().points(),
            intersection::HALF_RAY
        );

        if (inter.hit() && sqr(inter.distance()) <= minDistSqr)
        {
            // Note: no extra test on whether intersection is in front of us
            // since using half_ray AND zero tolerance. (note that tolerance
            // is used to look behind us)
            minDistSqr = sqr(inter.distance());
            intersectionPoint = inter.point();
            hasMin = true;
        }
    }

    // Restore picking tolerance
    intersection::setPlanarTol(oldTol);

    return hasMin;
}


// ************************************************************************* //
