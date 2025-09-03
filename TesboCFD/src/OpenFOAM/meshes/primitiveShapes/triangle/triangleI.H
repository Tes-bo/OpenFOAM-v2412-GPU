/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "IOstreams.H"
#include "pointHit.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

inline Foam::triPoints::triPoints
(
    const point& p0,
    const point& p1,
    const point& p2
)
{
    a() = p0;
    b() = p1;
    c() = p2;
}


inline Foam::triPoints::triPoints(const triPointRef& pts)
{
    a() = pts.a();
    b() = pts.b();
    c() = pts.c();
}


inline Foam::triPoints::triPoints(const FixedList<point, 3>& pts)
:
    FixedList<point, 3>(pts)
{}


inline Foam::triPoints::triPoints
(
    const UList<point>& points,
    const FixedList<label, 3>& indices
)
:
    FixedList<point, 3>(points, indices)
{}


inline Foam::triPoints::triPoints
(
    const UList<point>& points,
    const label p0,
    const label p1,
    const label p2
)
{
    a() = points[p0];
    b() = points[p1];
    c() = points[p2];
}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const Point& p0,
    const Point& p1,
    const Point& p2
)
:
    a_(p0),
    b_(p1),
    c_(p2)
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const FixedList<Point, 3>& pts
)
:
    a_(pts.template get<0>()),
    b_(pts.template get<1>()),
    c_(pts.template get<2>())
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const UList<Point>& points,
    const FixedList<label, 3>& indices
)
:
    a_(points[indices.template get<0>()]),
    b_(points[indices.template get<1>()]),
    c_(points[indices.template get<2>()])
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const UList<Point>& points,
    const label p0,
    const label p1,
    const label p2
)
:
    a_(points[p0]),
    b_(points[p1]),
    c_(points[p2])
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::centre
(
    const Point& p0,
    const Point& p1,
    const Point& p2
)
{
    return (1.0/3.0)*(p0 + p1 + p2);
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::centre() const
{
    return (1.0/3.0)*(a_ + b_ + c_);
}


inline Foam::point Foam::triPoints::centre() const
{
    return triPointRef::centre(a(), b(), c());
}


template<class Point, class PointRef>
inline Foam::vector Foam::triangle<Point, PointRef>::areaNormal
(
    const Point& p0,
    const Point& p1,
    const Point& p2
)
{
    return 0.5*((p1 - p0)^(p2 - p0));
}


template<class Point, class PointRef>
inline Foam::vector Foam::triangle<Point, PointRef>::areaNormal() const
{
    return 0.5*((b_ - a_)^(c_ - a_));
}


inline Foam::vector Foam::triPoints::areaNormal() const
{
    return triPointRef::areaNormal(a(), b(), c());
}


template<class Point, class PointRef>
inline Foam::vector Foam::triangle<Point, PointRef>::unitNormal
(
    const Point& p0,
    const Point& p1,
    const Point& p2
)
{
    vector n(areaNormal(p0, p1, p2));
    (void) n.normalise(ROOTVSMALL);
    return n;
}


template<class Point, class PointRef>
inline Foam::vector Foam::triangle<Point, PointRef>::unitNormal() const
{
    return triangle<Point, PointRef>::unitNormal(a_, b_, c_);
}


inline Foam::vector Foam::triPoints::unitNormal() const
{
    return triPointRef::unitNormal(a(), b(), c());
}


template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::triangle<Point, PointRef>::box
(
    const Point& p0,
    const Point& p1,
    const Point& p2
)
{
    return Pair<Point>
    (
        min(p0, min(p1, p2)),
        max(p0, max(p1, p2))
    );
}


template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::triangle<Point, PointRef>::box() const
{
    return triangle<Point, PointRef>::box(a_, b_, c_);
}


inline Foam::Pair<Foam::point> Foam::triPoints::box() const
{
    return triPointRef::box(a(), b(), c());
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::vecA() const
{
    return Point(c_ - b_);
}

template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::vecB() const
{
    return Point(a_ - c_);
}

template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::vecC() const
{
    return Point(b_ - a_);
}


inline Foam::vector Foam::triPoints::vecA() const
{
    return vector(c() - b());
}


inline Foam::vector Foam::triPoints::vecB() const
{
    return vector(a() - c());
}


inline Foam::vector Foam::triPoints::vecC() const
{
    return vector(b() - a());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::triPointRef Foam::triPoints::tri() const
{
    return triPointRef(a(), b(), c());
}


inline void Foam::triPoints::flip()
{
    // swap pt1 <-> pt2
    point t(b());
    b() = c();
    c() = t;
}


inline Foam::scalar Foam::triPoints::mag() const
{
    return ::Foam::mag(areaNormal());
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::mag() const
{
    return ::Foam::mag(areaNormal());
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::circumCentre() const
{
    scalar d1 =  (c_ - a_) & (b_ - a_);
    scalar d2 = -(c_ - b_) & (b_ - a_);
    scalar d3 =  (c_ - a_) & (c_ - b_);

    scalar c1 = d2*d3;
    scalar c2 = d3*d1;
    scalar c3 = d1*d2;

    scalar c = c1 + c2 + c3;

    if (Foam::mag(c) < ROOTVSMALL)
    {
        // Degenerate triangle, returning centre instead of circumCentre.

        return centre();
    }

    return
    (
        ((c2 + c3)*a_ + (c3 + c1)*b_ + (c1 + c2)*c_)/(2*c)
    );
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::circumRadius() const
{
    const scalar d1 =  (c_ - a_) & (b_ - a_);
    const scalar d2 = -(c_ - b_) & (b_ - a_);
    const scalar d3 =  (c_ - a_) & (c_ - b_);

    const scalar denom = d2*d3 + d3*d1 + d1*d2;

    if (Foam::mag(denom) < VSMALL)
    {
        // Degenerate triangle, returning GREAT for circumRadius.

        return GREAT;
    }

    const scalar a = (d1 + d2)*(d2 + d3)*(d3 + d1) / denom;
    return 0.5*Foam::sqrt(clamp(a, scalar(0), scalar(GREAT)));
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::quality() const
{
    scalar c = circumRadius();

    if (c < ROOTVSMALL)
    {
        // zero circumRadius, something has gone wrong.
        return SMALL;
    }

    return mag()/(Foam::sqr(c)*3.0*sqrt(3.0)/4.0);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::sweptVol
(
    const triangle& t
) const
{
    return (1.0/12.0)*
    (
        ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((t.b_ - b_) & ((c_ - b_)^(t.a_ - b_)))
      + ((c_ - t.c_) & ((t.b_ - t.c_)^(t.a_ - t.c_)))

      + ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((b_ - t.b_) & ((t.a_ - t.b_)^(t.c_ - t.b_)))
      + ((c_ - t.c_) & ((b_ - t.c_)^(t.a_ - t.c_)))
    );
}


template<class Point, class PointRef>
inline Foam::tensor Foam::triangle<Point, PointRef>::inertia
(
    PointRef refPt,
    scalar density
) const
{
    Point aRel = a_ - refPt;
    Point bRel = b_ - refPt;
    Point cRel = c_ - refPt;

    tensor V
    (
        aRel.x(), aRel.y(), aRel.z(),
        bRel.x(), bRel.y(), bRel.z(),
        cRel.x(), cRel.y(), cRel.z()
    );

    scalar a = Foam::mag((b_ - a_)^(c_ - a_));

    tensor S = 1/24.0*(tensor::one + I);

    return
    (
        a*I/24.0*
        (
            (aRel & aRel)
          + (bRel & bRel)
          + (cRel & cRel)
          + ((aRel + bRel + cRel) & (aRel + bRel + cRel))
        )
      - a*(V.T() & S & V)
    )
   *density;
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::randomPoint(Random& rndGen) const
{
    return barycentricToPoint(barycentric2D01(rndGen));
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::barycentricToPoint
(
    const barycentric2D& bary
) const
{
    return bary[0]*a_ + bary[1]*b_ + bary[2]*c_;
}


template<class Point, class PointRef>
inline Foam::barycentric2D Foam::triangle<Point, PointRef>::pointToBarycentric
(
    const point& pt
) const
{
    barycentric2D bary;
    pointToBarycentric(pt, bary);
    return bary;
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::pointToBarycentric
(
    const point& pt,
    barycentric2D& bary
) const
{
    // Reference:
    // Real-time collision detection, Christer Ericson, 2005, p47-48

    vector v0 = b_ - a_;
    vector v1 = c_ - a_;
    vector v2 = pt - a_;

    scalar d00 = v0 & v0;
    scalar d01 = v0 & v1;
    scalar d11 = v1 & v1;
    scalar d20 = v2 & v0;
    scalar d21 = v2 & v1;

    scalar denom = d00*d11 - d01*d01;

    if (Foam::mag(denom) < SMALL)
    {
        // Degenerate triangle, returning 1/3 barycentric coordinates.

        bary = barycentric2D(1.0/3.0, 1.0/3.0, 1.0/3.0);

        return denom;
    }

    bary[1] = (d11*d20 - d01*d21)/denom;
    bary[2] = (d00*d21 - d01*d20)/denom;
    bary[0] = 1.0 - bary[1] - bary[2];

    return denom;
}


template<class Point, class PointRef>
inline bool Foam::triangle<Point, PointRef>::intersects
(
    const point& origin,
    const vector& normal
) const
{
    // Check points cut below(1) or above(2). Mix of above/below == 3
    // Cf. plane::whichSide()
    return
    (
        (
            (((a_ - origin) & normal) < 0 ? 1 : 2)
          | (((b_ - origin) & normal) < 0 ? 1 : 2)
          | (((c_ - origin) & normal) < 0 ? 1 : 2)
        ) == 3
    );
}


template<class Point, class PointRef>
inline bool Foam::triangle<Point, PointRef>::intersects
(
    const point& origin,
    const vector::components axis
) const
{
    // Direct check of points cut below(1) or above(2)
    // Cf. plane::whichSide()
    return
    (
        (
            (a_[axis] < origin[axis] ? 1 : 2)
          | (b_[axis] < origin[axis] ? 1 : 2)
          | (c_[axis] < origin[axis] ? 1 : 2)
        ) == 3
    );
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::ray
(
    const point& p,
    const vector& q,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Express triangle in terms of baseVertex (point a_) and
    // two edge vectors
    vector E0 = b_ - a_;
    vector E1 = c_ - a_;

    // Initialize intersection to miss.
    pointHit inter(p);

    vector n(0.5*(E0 ^ E1));
    scalar area = Foam::mag(n);

    if (area < VSMALL)
    {
        // Ineligible miss.
        inter.setMiss(false);

        // The miss point is the nearest point on the triangle. Take any one.
        inter.setPoint(a_);

        // The distance to the miss is the distance between the
        // original point and plane of intersection. No normal so use
        // distance from p to a_
        inter.setDistance(Foam::mag(a_ - p));

        return inter;
    }

    vector q1 = q/Foam::mag(q);

    if (dir == intersection::CONTACT_SPHERE)
    {
        n /= area;

        return ray(p, q1 - n, alg, intersection::VECTOR);
    }

    // Intersection point with triangle plane
    point pInter;
    // Is intersection point inside triangle
    bool hit;
    {
        // Reuse the fast ray intersection routine below in FULL_RAY
        // mode since the original intersection routine has rounding problems.
        pointHit fastInter = intersection(p, q1, intersection::FULL_RAY);
        hit = fastInter.hit();

        if (hit)
        {
            pInter = fastInter.point();
        }
        else
        {
            // Calculate intersection of ray with triangle plane
            vector v = a_ - p;
            pInter = p + (q1&v)*q1;
        }
    }

    // Distance to intersection point
    scalar dist = q1 & (pInter - p);

    const scalar planarPointTol =
        Foam::min
        (
            Foam::min
            (
                Foam::mag(E0),
                Foam::mag(E1)
            ),
            Foam::mag(c_ - b_)
        )*intersection::planarTol();

    bool eligible =
        alg == intersection::FULL_RAY
     || (alg == intersection::HALF_RAY && dist > -planarPointTol)
     || (
            alg == intersection::VISIBLE
         && ((q1 & areaNormal()) < -VSMALL)
        );

    if (hit && eligible)
    {
        // Hit. Set distance to intersection.
        inter.hitPoint(pInter);
        inter.setDistance(dist);
    }
    else
    {
        // Miss or ineligible hit.
        inter.setMiss(eligible);

        // The miss point is the nearest point on the triangle
        inter.setPoint(nearestPoint(p).point());

        // The distance to the miss is the distance between the
        // original point and plane of intersection
        inter.setDistance(Foam::mag(pInter - p));
    }


    return inter;
}


// From "Fast, Minimum Storage Ray/Triangle Intersection"
// Moeller/Trumbore.
template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::intersection
(
    const point& orig,
    const vector& dir,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    const vector edge1 = b_ - a_;
    const vector edge2 = c_ - a_;

    // begin calculating determinant - also used to calculate U parameter
    const vector pVec = dir ^ edge2;

    // if determinant is near zero, ray lies in plane of triangle
    const scalar det = edge1 & pVec;

    // Initialise to miss
    pointHit intersection(false, Zero, GREAT, false);

    if (alg == intersection::VISIBLE)
    {
        // Culling branch
        if (det < ROOTVSMALL)
        {
            // Ray on wrong side of triangle. Return miss
            return intersection;
        }
    }
    else if (alg == intersection::HALF_RAY || alg == intersection::FULL_RAY)
    {
        // Non-culling branch
        if (det > -ROOTVSMALL && det < ROOTVSMALL)
        {
            // Ray parallel to triangle. Return miss
            return intersection;
        }
    }

    const scalar inv_det = 1.0 / det;

    /* calculate distance from a_ to ray origin */
    const vector tVec = orig-a_;

    /* calculate U parameter and test bounds */
    const scalar u = (tVec & pVec)*inv_det;

    if (u < -tol || u > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* prepare to test V parameter */
    const vector qVec = tVec ^ edge1;

    /* calculate V parameter and test bounds */
    const scalar v = (dir & qVec) * inv_det;

    if (v < -tol || u + v > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* calculate t, scale parameters, ray intersects triangle */
    const scalar t = (edge2 & qVec) * inv_det;

    if (alg == intersection::HALF_RAY && t < -tol)
    {
        // Wrong side of orig. Return miss
        return intersection;
    }

    intersection.hitPoint(a_ + u*edge1 + v*edge2);
    intersection.setDistance(t);

    return intersection;
}


template<class Point, class PointRef>
Foam::pointHit Foam::triangle<Point, PointRef>::nearestPointClassify
(
    const point& p,
    label& nearType,
    label& nearLabel
) const
{
    // Adapted from:
    // Real-time collision detection, Christer Ericson, 2005, p136-142

    // Check if P in vertex region outside A
    vector ab = b_ - a_;
    vector ac = c_ - a_;
    vector ap = p - a_;

    scalar d1 = ab & ap;
    scalar d2 = ac & ap;

    if (d1 <= 0.0 && d2 <= 0.0)
    {
        // barycentric coordinates (1, 0, 0)

        nearType = POINT;
        nearLabel = 0;
        return pointHit(false, a_, Foam::mag(a_ - p), true);
    }

    // Check if P in vertex region outside B
    vector bp = p - b_;
    scalar d3 = ab & bp;
    scalar d4 = ac & bp;

    if (d3 >= 0.0 && d4 <= d3)
    {
        // barycentric coordinates (0, 1, 0)

        nearType = POINT;
        nearLabel = 1;
        return pointHit(false, b_, Foam::mag(b_ - p), true);
    }

    // Check if P in edge region of AB, if so return projection of P onto AB
    scalar vc = d1*d4 - d3*d2;

    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        if ((d1 - d3) < ROOTVSMALL)
        {
            // Degenerate triangle, for d1 = d3, a_ and b_ are likely coincident
            nearType = POINT;
            nearLabel = 0;
            return pointHit(false, a_, Foam::mag(a_ - p), true);
        }

        // barycentric coordinates (1-v, v, 0)
        scalar v = d1/(d1 - d3);

        point nearPt =  a_ + v*ab;
        nearType = EDGE;
        nearLabel = 0;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // Check if P in vertex region outside C
    vector cp = p - c_;
    scalar d5 = ab & cp;
    scalar d6 = ac & cp;

    if (d6 >= 0.0 && d5 <= d6)
    {
        // barycentric coordinates (0, 0, 1)

        nearType = POINT;
        nearLabel = 2;
        return pointHit(false, c_, Foam::mag(c_ - p), true);
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    scalar vb = d5*d2 - d1*d6;

    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        if ((d2 - d6) < ROOTVSMALL)
        {
            // Degenerate triangle, for d2 = d6, a_ and c_ are likely coincident
            nearType = POINT;
            nearLabel = 0;
            return pointHit(false, a_, Foam::mag(a_ - p), true);
        }

        // barycentric coordinates (1-w, 0, w)
        scalar w = d2/(d2 - d6);

        point nearPt = a_ + w*ac;
        nearType = EDGE;
        nearLabel = 2;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    scalar va = d3*d6 - d5*d4;

    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        if (((d4 - d3) + (d5 - d6)) < ROOTVSMALL)
        {
            // Degenerate triangle, for (d4 - d3) = (d6 - d5), b_ and c_ are
            // likely coincident
            nearType = POINT;
            nearLabel = 1;
            return pointHit(false, b_, Foam::mag(b_ - p), true);
        }

        // barycentric coordinates (0, 1-w, w)
        scalar w = (d4 - d3)/((d4 - d3) + (d5 - d6));

        point nearPt = b_ + w*(c_ - b_);
        nearType = EDGE;
        nearLabel = 1;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // P inside face region. Compute Q through its barycentric
    // coordinates (u, v, w)

    if ((va + vb + vc) < ROOTVSMALL)
    {
        // Degenerate triangle, return the centre because no edge or points are
        // closest
        point nearPt = centre();
        nearType = NONE,
        nearLabel = -1;
        return pointHit(true, nearPt, Foam::mag(nearPt - p), false);
    }

    scalar denom = 1.0/(va + vb + vc);
    scalar v = vb * denom;
    scalar w = vc * denom;

    // = u*a + v*b + w*c, u = va*denom = 1.0 - v - w

    point nearPt = a_ + ab*v + ac*w;
    nearType = NONE,
    nearLabel = -1;
    return pointHit(true, nearPt, Foam::mag(nearPt - p), false);
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::nearestPoint
(
    const point& p
) const
{
    // Dummy labels
    label nearType = -1;
    label nearLabel = -1;

    return nearestPointClassify(p, nearType, nearLabel);
}


template<class Point, class PointRef>
inline bool Foam::triangle<Point, PointRef>::classify
(
    const point& p,
    label& nearType,
    label& nearLabel
) const
{
    return nearestPointClassify(p, nearType, nearLabel).hit();
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::nearestPoint
(
    const linePointRef& ln,
    pointHit& lnInfo
) const
{
    vector q = ln.vec();
    pointHit triInfo
    (
        triangle<Point, PointRef>::intersection
        (
            ln.start(),
            q,
            intersection::FULL_RAY
        )
    );

    if (triInfo.hit())
    {
        // Line hits triangle. Find point on line.
        if (triInfo.distance() > 1)
        {
            // Hit beyond endpoint
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.end());
            scalar dist = triInfo.point().dist(lnInfo.point());
            lnInfo.setDistance(dist);
            triInfo.setMiss(true);
            triInfo.setDistance(dist);
        }
        else if (triInfo.distance() < 0)
        {
            // Hit beyond startpoint
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.start());
            scalar dist = triInfo.point().dist(lnInfo.point());
            lnInfo.setDistance(dist);
            triInfo.setMiss(true);
            triInfo.setDistance(dist);
        }
        else
        {
            // Hit on line
            lnInfo.hitPoint(triInfo.point());
            lnInfo.setDistance(0);
            triInfo.setDistance(0);
        }
    }
    else
    {
        // Line skips triangle. See which triangle edge it gets closest to

        point nearestEdgePoint;
        point nearestLinePoint;
        //label minEdgeIndex = 0;
        scalar minDist = ln.nearestDist
        (
            linePointRef(a_, b_),
            nearestLinePoint,
            nearestEdgePoint
        );

        {
            point linePoint;
            point triEdgePoint;
            scalar dist = ln.nearestDist
            (
                linePointRef(b_, c_),
                linePoint,
                triEdgePoint
            );
            if (dist < minDist)
            {
                minDist = dist;
                nearestEdgePoint = triEdgePoint;
                nearestLinePoint = linePoint;
                //minEdgeIndex = 1;
            }
        }

        {
            point linePoint;
            point triEdgePoint;
            scalar dist = ln.nearestDist
            (
                linePointRef(c_, a_),
                linePoint,
                triEdgePoint
            );
            if (dist < minDist)
            {
                minDist = dist;
                nearestEdgePoint = triEdgePoint;
                nearestLinePoint = linePoint;
                //minEdgeIndex = 2;
            }
        }

        lnInfo.setDistance(minDist);
        triInfo.setDistance(minDist);
        triInfo.setMiss(false);
        triInfo.setPoint(nearestEdgePoint);

        // Convert point on line to pointHit
        if (Foam::mag(nearestLinePoint-ln.start()) < SMALL)
        {
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.start());
        }
        else if (Foam::mag(nearestLinePoint-ln.end()) < SMALL)
        {
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.end());
        }
        else
        {
            lnInfo.hitPoint(nearestLinePoint);
        }
    }
    return triInfo;
}


template<class Point, class PointRef>
inline int Foam::triangle<Point, PointRef>::sign
(
    const point& p,
    const scalar tol
) const
{
    const scalar dist = ((p - a_) & unitNormal());

    return ((dist < -tol) ? -1 : (dist > tol) ? +1 : 0);
}


template<class Point, class PointRef>
inline void Foam::triangle<Point, PointRef>::dummyOp::operator()
(
    const triPoints&
)
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::sumAreaOp::sumAreaOp()
:
    area_(0.0)
{}


template<class Point, class PointRef>
inline void Foam::triangle<Point, PointRef>::sumAreaOp::operator()
(
    const triPoints& tri
)
{
    area_ += triangle<Point, const Point&>(tri).mag();
}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::storeOp::storeOp
(
    triIntersectionList& tris,
    label& nTris
)
:
    tris_(tris),
    nTris_(nTris)
{}


template<class Point, class PointRef>
inline void Foam::triangle<Point, PointRef>::storeOp::operator()
(
    const triPoints& tri
)
{
    tris_[nTris_++] = tri;
}


template<class Point, class PointRef>
inline Foam::point Foam::triangle<Point, PointRef>::planeIntersection
(
    const FixedList<scalar, 3>& d,
    const triPoints& t,
    const label negI,
    const label posI
)
{
    return (d[posI]*t[negI] - d[negI]*t[posI])/(-d[negI] + d[posI]);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    triangle<Point, PointRef>& t
)
{
    is.readBegin("triangle");
    is  >> t.a_ >> t.b_ >> t.c_;
    is.readEnd("triangle");

    is.check(FUNCTION_NAME);
    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const triangle<Point, PointRef>& t
)
{
    // Format like FixedList
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE
        << t.b_ << token::SPACE
        << t.c_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
