/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "triangle.H"
#include "plane.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
template<class AboveOp, class BelowOp>
inline void Foam::triangle<Point, PointRef>::triSliceWithPlane
(
    const plane& pln,
    const triPoints& tri,
    AboveOp& aboveOp,
    BelowOp& belowOp
)
{
    // distance to cutting plane
    FixedList<scalar, 3> d;

    // determine how many of the points are above the cutting plane
    label nPos = 0;
    label posI = -1;
    label negI = -1;
    forAll(tri, i)
    {
        d[i] = pln.signedDistance(tri[i]);

        if (d[i] > 0)
        {
            nPos++;
            posI = i;
        }
        else
        {
            negI = i;
        }
    }

    if (nPos == 3)
    {
        aboveOp(tri);
    }
    else if (nPos == 2)
    {
        // point under the plane
        label i0 = negI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        aboveOp(triPoints(tri[i1], tri[i2], p02));
        aboveOp(triPoints(tri[i1], p02, p01));
        belowOp(triPoints(tri[i0], p01, p02));
    }
    else if (nPos == 1)
    {
        // point above the plane
        label i0 = posI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        belowOp(triPoints(tri[i1], tri[i2], p02));
        belowOp(triPoints(tri[i1], p02, p01));
        aboveOp(triPoints(tri[i0], p01, p02));
    }
    else
    {
        // All below
        belowOp(tri);
    }
}


template<class Point, class PointRef>
template<class AboveOp, class BelowOp>
inline void Foam::triangle<Point, PointRef>::sliceWithPlane
(
    const plane& pl,
    AboveOp& aboveOp,
    BelowOp& belowOp
) const
{
    triSliceWithPlane(pl, triPoints(a_, b_, c_), aboveOp, belowOp);
}


template<class Point, class PointRef>
template<class InsideOp, class OutsideOp>
inline void Foam::triangle<Point, PointRef>::triangleOverlap
(
    const vector& n,
    const triangle<Point, PointRef>& tgt,
    InsideOp& insideOp,
    OutsideOp& outsideOp
) const
{
    // There are two possibilities with this algorithm - we either cut
    // the outside triangles with all the edges or not (and keep them
    // as disconnected triangles). In the first case
    // we cannot do any evaluation short cut so we've chosen not to re-cut
    // the outside triangles.


    triIntersectionList insideTrisA;
    label nInsideA = 0;
    storeOp insideOpA(insideTrisA, nInsideA);

    triIntersectionList outsideTrisA;
    label nOutsideA = 0;
    storeOp outsideOpA(outsideTrisA, nOutsideA);


    const triPoints thisTri(a_, b_, c_);


    // Cut original triangle with tgt edge 0.
    // From *this to insideTrisA, outsideTrisA.
    {
        scalar s = Foam::mag(tgt.b() - tgt.a());
        const plane pl0(tgt.a(), tgt.b(), tgt.b() + s*n);
        triSliceWithPlane(pl0, thisTri, insideOpA, outsideOpA);
    }

    // Shortcut if nothing cut

    if (insideOpA.nTris_ == 0)
    {
        outsideOp(thisTri);
        return;
    }

    if (outsideOpA.nTris_ == 0)
    {
        insideOp(thisTri);
        return;
    }


    // Cut all triangles with edge 1.
    // From insideTrisA to insideTrisB, outsideTrisA

    triIntersectionList insideTrisB;
    label nInsideB = 0;
    storeOp insideOpB(insideTrisB, nInsideB);

    //triIntersectionList outsideTrisB;
    //label nOutsideB = 0;
    //storeOp outsideOpB(outsideTrisB, nOutsideB);

    {
        scalar s = Foam::mag(tgt.c() - tgt.b());
        const plane pl0(tgt.b(), tgt.c(), tgt.c() + s*n);

        for (label i = 0; i < insideOpA.nTris_; i++)
        {
            const triPoints& tri = insideOpA.tris_[i];
            triSliceWithPlane(pl0, tri, insideOpB, outsideOpA);
        }

        //// Recut outside triangles (not necessary if only interested in
        //// intersection properties)
        //for (label i = 0; i < outsideOpA.nTris_; i++)
        //{
        //    const triPoints& tri = outsideOpA.tris_[i];
        //    triSliceWithPlane(pl0, tri, outsideOpB, outsideOpB);
        //}
    }


    // Cut all triangles with edge 2.
    // From insideTrisB to insideTrisA, outsideTrisA
    {
        scalar s = Foam::mag(tgt.a() - tgt.c());
        const plane pl0(tgt.c(), tgt.a(), tgt.a() + s*n);

        insideOpA.nTris_ = 0;
        //outsideOpA.nTris_ = 0;
        for (label i = 0; i < insideOpB.nTris_; i++)
        {
            const triPoints& tri = insideOpB.tris_[i];
            triSliceWithPlane(pl0, tri, insideOpA, outsideOpA);
        }

        //// Recut outside triangles (not necessary if only interested in
        //// intersection properties)
        //for (label i = 0; i < outsideOpB.nTris_; i++)
        //{
        //    const triPoints& tri = outsideOpB.tris_[i];
        //    triSliceWithPlane(pl0, tri, outsideOpA, outsideOpA);
        //}
    }

    // Transfer from A to argument
    for (label i = 0; i < insideOpA.nTris_; i++)
    {
        insideOp(insideOpA.tris_[i]);
    }

    for (label i = 0; i < outsideOpA.nTris_; i++)
    {
        outsideOp(outsideOpA.tris_[i]);
    }
}


// ************************************************************************* //
