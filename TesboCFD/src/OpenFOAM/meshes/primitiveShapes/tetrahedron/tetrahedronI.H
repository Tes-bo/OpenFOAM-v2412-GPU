/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "triangle.H"
#include "IOstreams.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::tetPoints::tetPoints
(
    const point& p0,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    a() = p0;
    b() = p1;
    c() = p2;
    d() = p3;
}


inline Foam::tetPoints::tetPoints(const tetPointRef& pts)
{
    a() = pts.a();
    b() = pts.b();
    c() = pts.c();
    d() = pts.d();
}


inline Foam::tetPoints::tetPoints(const FixedList<point, 4>& pts)
:
    FixedList<point, 4>(pts)
{}


inline Foam::tetPoints::tetPoints
(
    const UList<point>& points,
    const FixedList<label, 4>& indices
)
:
    FixedList<point, 4>(points, indices)
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron
(
    const Point& p0,
    const Point& p1,
    const Point& p2,
    const Point& p3
)
:
    a_(p0),
    b_(p1),
    c_(p2),
    d_(p3)
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron
(
    const FixedList<Point, 4>& pts
)
:
    a_(pts.template get<0>()),
    b_(pts.template get<1>()),
    c_(pts.template get<2>()),
    d_(pts.template get<3>())
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron
(
    const UList<Point>& points,
    const FixedList<label, 4>& indices
)
:
    a_(points[indices.template get<0>()]),
    b_(points[indices.template get<1>()]),
    c_(points[indices.template get<2>()]),
    d_(points[indices.template get<3>()])
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::tetPointRef Foam::tetPoints::tet() const
{
    return tetPointRef(a(), b(), c(), d());
}


inline void Foam::tetPoints::flip()
{
    // swap pt2 <-> pt3
    point t(c());
    c() = d();
    d() = t;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::tetrahedron<Point, PointRef>::box
(
    const Point& p0,
    const Point& p1,
    const Point& p2,
    const Point& p3
)
{
    return Pair<Point>
    (
        min(p0, min(p1, min(p2, p3))),
        max(p0, max(p1, max(p2, p3)))
    );
}


template<class Point, class PointRef>
inline Foam::Pair<Point> Foam::tetrahedron<Point, PointRef>::box() const
{
    return tetrahedron<Point, PointRef>::box(a_, b_, c_, d_);
}


inline Foam::Pair<Foam::point> Foam::tetPoints::box() const
{
    return tetPointRef::box(a(), b(), c(), d());
}


template<class Point, class PointRef>
inline Foam::treeBoundBox Foam::tetrahedron<Point, PointRef>::bounds() const
{
    return treeBoundBox(box());
}


inline Foam::treeBoundBox Foam::tetPoints::bounds() const
{
    return treeBoundBox(box());
}


// Warning. Ordering of faces needs to be the same for a tetrahedron class,
// tetrahedron cell shape model and a tetCell

template<class Point, class PointRef>
inline Foam::triPointRef Foam::tetrahedron<Point, PointRef>::tri
(
    const label facei
) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell

    if (facei == 0)
    {
        return triPointRef(b_, c_, d_);
    }
    else if (facei == 1)
    {
        return triPointRef(a_, d_, c_);
    }
    else if (facei == 2)
    {
        return triPointRef(a_, b_, d_);
    }
    else if (facei == 3)
    {
        return triPointRef(a_, c_, b_);
    }

    FatalErrorInFunction
        << "Face index (" << facei << ") out of range 0..3\n"
        << abort(FatalError);
    return triPointRef(b_, c_, d_);
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sa() const
{
    return triangle<Point, PointRef>::areaNormal(b_, c_, d_);
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sb() const
{
    return triangle<Point, PointRef>::areaNormal(a_, d_, c_);
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sc() const
{
    return triangle<Point, PointRef>::areaNormal(a_, b_, d_);
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sd() const
{
    return triangle<Point, PointRef>::areaNormal(a_, c_, b_);
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::centre() const
{
    return 0.25*(a_ + b_ + c_ + d_);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::mag() const
{
    return (1.0/6.0)*(((b_ - a_) ^ (c_ - a_)) & (d_ - a_));
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::circumCentre() const
{
    vector a = b_ - a_;
    vector b = c_ - a_;
    vector c = d_ - a_;

    scalar lambda = magSqr(c) - (a & c);
    scalar mu = magSqr(b) - (a & b);

    vector ba = b ^ a;
    vector ca = c ^ a;

    vector num = lambda*ba - mu*ca;
    scalar denom = (c & ba);

    if (Foam::mag(denom) < ROOTVSMALL)
    {
        // Degenerate tetrahedron, returning centre instead of circumCentre.

        return centre();
    }

    return a_ + 0.5*(a + num/denom);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::circumRadius() const
{
    vector a = b_ - a_;
    vector b = c_ - a_;
    vector c = d_ - a_;

    scalar lambda = magSqr(c) - (a & c);
    scalar mu = magSqr(b) - (a & b);

    vector ba = b ^ a;
    vector ca = c ^ a;

    vector num = lambda*ba - mu*ca;
    scalar denom = (c & ba);

    if (Foam::mag(denom) < ROOTVSMALL)
    {
        // Degenerate tetrahedron, returning GREAT for circumRadius.
        return GREAT;
    }

    return Foam::mag(0.5*(a + num/denom));
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::quality() const
{
    return
        mag()
       /(
            8.0/(9.0*sqrt(3.0))
           *pow3(Foam::min(circumRadius(), GREAT))
          + ROOTVSMALL
        );
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::randomPoint
(
    Random& rndGen
) const
{
    return barycentricToPoint(barycentric01(rndGen));
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::barycentricToPoint
(
    const barycentric& bary
) const
{
    return Point(bary.a()*a_ + bary.b()*b_ + bary.c()*c_ + bary.d()*d_);
}


template<class Point, class PointRef>
inline Foam::barycentric Foam::tetrahedron<Point, PointRef>::pointToBarycentric
(
    const point& pt
) const
{
    barycentric bary;
    pointToBarycentric(pt, bary);
    return bary;
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::pointToBarycentric
(
    const point& pt,
    barycentric& bary
) const
{
    // Reference:
    // http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)

    vector e0(a_ - d_);
    vector e1(b_ - d_);
    vector e2(c_ - d_);

    tensor t
    (
        e0.x(), e1.x(), e2.x(),
        e0.y(), e1.y(), e2.y(),
        e0.z(), e1.z(), e2.z()
    );

    scalar detT = det(t);

    if (Foam::mag(detT) < SMALL)
    {
        // Degenerate tetrahedron, returning 1/4 barycentric coordinates

        bary = barycentric(0.25, 0.25, 0.25, 0.25);

        return detT;
    }

    vector res = inv(t, detT) & (pt - d_);

    bary[0] = res.x();
    bary[1] = res.y();
    bary[2] = res.z();
    bary[3] = 1 - bary[0] - bary[1] - bary[2];

    return detT;
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::tetrahedron<Point, PointRef>::nearestPoint
(
    const point& p
) const
{
    // Adapted from:
    // Real-time collision detection, Christer Ericson, 2005, p142-144

    // Assuming initially that the point is inside all of the
    // halfspaces, so closest to itself.

    point closestPt = p;

    scalar minOutsideDistance = VGREAT;

    bool inside = true;

    // Side a
    {
        const triangle<Point, PointRef> tria(b_, c_, d_);

        if (((p - b_) & tria.areaNormal()) >= 0)
        {
            // p is outside halfspace plane of tri
            pointHit info = tria.nearestPoint(p);

            inside = false;

            if (info.distance() < minOutsideDistance)
            {
                closestPt = info.point();

                minOutsideDistance = info.distance();
            }
        }
    }

    // Side b
    {
        const triangle<Point, PointRef> tria(a_, d_, c_);

        if (((p - a_) & tria.areaNormal()) >= 0)
        {
            // p is outside halfspace plane of tri
            pointHit info = tria.nearestPoint(p);

            inside = false;

            if (info.distance() < minOutsideDistance)
            {
                closestPt = info.point();

                minOutsideDistance = info.distance();
            }
        }
    }

    // Side c
    {
        const triangle<Point, PointRef> tria(a_, b_, d_);

        if (((p - a_) & tria.areaNormal()) >= 0)
        {
            // p is outside halfspace plane of tri
            pointHit info = tria.nearestPoint(p);

            inside = false;

            if (info.distance() < minOutsideDistance)
            {
                closestPt = info.point();

                minOutsideDistance = info.distance();
            }
        }
    }

    // Side c
    {
        const triangle<Point, PointRef> tria(a_, c_, b_);

        if (((p - a_) & tria.areaNormal()) >= 0)
        {
            // p is outside halfspace plane of tri
            pointHit info = tria.nearestPoint(p);

            inside = false;

            if (info.distance() < minOutsideDistance)
            {
                closestPt = info.point();

                minOutsideDistance = info.distance();
            }
        }
    }

    // If the point is inside, then the distance to the closest point
    // is zero
    if (inside)
    {
        minOutsideDistance = 0;
    }

    return pointHit
    (
        inside,
        closestPt,
        minOutsideDistance,
        !inside
    );
}


template<class Point, class PointRef>
bool Foam::tetrahedron<Point, PointRef>::inside(const point& p) const
{
    // For robustness, assuming that the point is in the tet unless
    // "definitively" shown otherwise by obtaining a positive dot
    // product greater than a tolerance of SMALL.

    // The tet is defined: tet(Cc, tetBasePt, pA, pB) where the area normal
    // vectors and base points for the half-space planes are:
    // area[0] = Sa();
    // area[1] = Sb();
    // area[2] = Sc();
    // area[3] = Sd();
    // planeBase[0] = tetBasePt = b_
    // planeBase[1] = ptA       = c_
    // planeBase[2] = tetBasePt = b_
    // planeBase[3] = tetBasePt = b_

    // Side a
    {
        const triangle<Point, PointRef> tria(b_, c_, d_);

        if (((p - b_) & tria.unitNormal()) > SMALL)
        {
            return false;
        }
    }

    // Side b
    {
        const triangle<Point, PointRef> tria(a_, d_, c_);

        if (((p - a_) & tria.unitNormal()) > SMALL)
        {
            return false;
        }
    }

    // Side c
    {
        const triangle<Point, PointRef> tria(a_, b_, d_);

        if (((p - a_) & tria.unitNormal()) > SMALL)
        {
            return false;
        }
    }

    // Side d
    {
        const triangle<Point, PointRef> tria(a_, c_, b_);

        if (((p - a_) & tria.unitNormal()) > SMALL)
        {
            return false;
        }
    }

    return true;
}


template<class Point, class PointRef>
inline void Foam::tetrahedron<Point, PointRef>::dummyOp::operator()
(
    const tetPoints&
)
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::sumVolOp::sumVolOp()
:
    vol_(0.0)
{}


template<class Point, class PointRef>
inline void Foam::tetrahedron<Point, PointRef>::sumVolOp::operator()
(
    const tetPoints& tet
)
{
    vol_ += tet.tet().mag();
}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::storeOp::storeOp
(
    tetIntersectionList& tets,
    label& nTets
)
:
    tets_(tets),
    nTets_(nTets)
{}


template<class Point, class PointRef>
inline void Foam::tetrahedron<Point, PointRef>::storeOp::operator()
(
    const tetPoints& tet
)
{
    tets_[nTets_++] = tet;
}


template<class Point, class PointRef>
inline Foam::point Foam::tetrahedron<Point, PointRef>::planeIntersection
(
    const FixedList<scalar, 4>& d,
    const tetPoints& t,
    const label negI,
    const label posI
)
{
    return
        (d[posI]*t[negI] - d[negI]*t[posI])
      / (-d[negI]+d[posI]);
}


template<class Point, class PointRef>
template<class TetOp>
inline void Foam::tetrahedron<Point, PointRef>::decomposePrism
(
    const FixedList<point, 6>& points,
    TetOp& op
)
{
    op(tetPoints(points[1], points[3], points[2], points[0]));
    op(tetPoints(points[1], points[2], points[3], points[4]));
    op(tetPoints(points[4], points[2], points[3], points[5]));
}


template<class Point, class PointRef>
template<class AboveTetOp, class BelowTetOp>
inline void Foam::tetrahedron<Point, PointRef>::tetSliceWithPlane
(
    const plane& pln,
    const tetPoints& tet,
    AboveTetOp& aboveOp,
    BelowTetOp& belowOp
)
{
    // Distance to plane
    FixedList<scalar, 4> d;
    label nPos = 0;
    forAll(tet, i)
    {
        d[i] = pln.signedDistance(tet[i]);
        if (d[i] > 0)
        {
            ++nPos;
        }
    }

    if (nPos == 4)
    {
        aboveOp(tet);
    }
    else if (nPos == 3)
    {
        // Sliced into below tet and above prism.
        // Prism gets split into two tets.

        // Find the below tet
        label i0 = -1;
        forAll(d, i)
        {
            if (d[i] <= 0)
            {
                i0 = i;
                break;
            }
        }

        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);
        label i3 = d.fcIndex(i2);

        point p01(planeIntersection(d, tet, i0, i1));
        point p02(planeIntersection(d, tet, i0, i2));
        point p03(planeIntersection(d, tet, i0, i3));

        // i0 = tetCell vertex 0: p01,p02,p03 outwards pointing triad
        //          ,,         1 :     ,,     inwards pointing triad
        //          ,,         2 :     ,,     outwards pointing triad
        //          ,,         3 :     ,,     inwards pointing triad

        //Pout<< "Split 3pos tet " << tet << " d:" << d << " into" << nl;

        if (i0 == 0 || i0 == 2)
        {
            tetPoints t(tet[i0], p01, p02, p03);
            //Pout<< "    belowtet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 3, belowTet i0==0 or 2");
            belowOp(t);

            // Prism
            FixedList<point, 6> p
            (
                {
                    tet[i1],
                    tet[i3],
                    tet[i2],
                    p01,
                    p03,
                    p02
                }
            );

            //Pout<< "    aboveprism:" << p << endl;
            decomposePrism(p, aboveOp);
        }
        else
        {
            tetPoints t(p01, p02, p03, tet[i0]);
            //Pout<< "    belowtet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 3, belowTet i0==1 or 3");
            belowOp(t);

            // Prism
            FixedList<point, 6> p
            (
                {
                    tet[i3],
                    tet[i1],
                    tet[i2],
                    p03,
                    p01,
                    p02
                }
            );
            //Pout<< "    aboveprism:" << p << endl;
            decomposePrism(p, aboveOp);
        }
    }
    else if (nPos == 2)
    {
        // Tet cut into two prisms. Determine the positive one.
        label pos0 = -1;
        label pos1 = -1;
        forAll(d, i)
        {
            if (d[i] > 0)
            {
                if (pos0 == -1)
                {
                    pos0 = i;
                }
                else
                {
                    pos1 = i;
                }
            }
        }

        //Pout<< "Split 2pos tet " << tet << " d:" << d
        //    << " around pos0:" << pos0 << " pos1:" << pos1
        //    << " neg0:" << neg0 << " neg1:" << neg1 << " into" << nl;

        const edge posEdge(pos0, pos1);

        if (posEdge == edge(0, 1))
        {
            point p02(planeIntersection(d, tet, 0, 2));
            point p03(planeIntersection(d, tet, 0, 3));
            point p12(planeIntersection(d, tet, 1, 2));
            point p13(planeIntersection(d, tet, 1, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[0],
                        p02,
                        p03,
                        tet[1],
                        p12,
                        p13
                    }
                );

                //Pout<< "    01 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[2],
                        p02,
                        p12,
                        tet[3],
                        p03,
                        p13
                    }
                );

                //Pout<< "    01 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(1, 2))
        {
            point p01(planeIntersection(d, tet, 0, 1));
            point p13(planeIntersection(d, tet, 1, 3));
            point p02(planeIntersection(d, tet, 0, 2));
            point p23(planeIntersection(d, tet, 2, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[1],
                        p01,
                        p13,
                        tet[2],
                        p02,
                        p23
                    }
                );

                //Pout<< "    12 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[3],
                        p23,
                        p13,
                        tet[0],
                        p02,
                        p01
                    }
                );

                //Pout<< "    12 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(2, 0))
        {
            point p01(planeIntersection(d, tet, 0, 1));
            point p03(planeIntersection(d, tet, 0, 3));
            point p12(planeIntersection(d, tet, 1, 2));
            point p23(planeIntersection(d, tet, 2, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[2],
                        p12,
                        p23,
                        tet[0],
                        p01,
                        p03
                    }
                );

                //Pout<< "    20 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[1],
                        p12,
                        p01,
                        tet[3],
                        p23,
                        p03
                    }
                );

                //Pout<< "    20 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(0, 3))
        {
            point p01(planeIntersection(d, tet, 0, 1));
            point p02(planeIntersection(d, tet, 0, 2));
            point p13(planeIntersection(d, tet, 1, 3));
            point p23(planeIntersection(d, tet, 2, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[3],
                        p23,
                        p13,
                        tet[0],
                        p02,
                        p01
                    }
                );

                //Pout<< "    03 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[2],
                        p23,
                        p02,
                        tet[1],
                        p13,
                        p01
                    }
                );

                //Pout<< "    03 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(1, 3))
        {
            point p01(planeIntersection(d, tet, 0, 1));
            point p12(planeIntersection(d, tet, 1, 2));
            point p03(planeIntersection(d, tet, 0, 3));
            point p23(planeIntersection(d, tet, 2, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[1],
                        p12,
                        p01,
                        tet[3],
                        p23,
                        p03
                    }
                );

                //Pout<< "    13 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[2],
                        p12,
                        p23,
                        tet[0],
                        p01,
                        p03
                    }
                );

                //Pout<< "    13 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(2, 3))
        {
            point p02(planeIntersection(d, tet, 0, 2));
            point p12(planeIntersection(d, tet, 1, 2));
            point p03(planeIntersection(d, tet, 0, 3));
            point p13(planeIntersection(d, tet, 1, 3));
            // Split the resulting prism
            {
                FixedList<point, 6> p
                (
                    {
                        tet[2],
                        p02,
                        p12,
                        tet[3],
                        p03,
                        p13
                    }
                );

                //Pout<< "    23 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p
                (
                    {
                        tet[0],
                        p02,
                        p03,
                        tet[1],
                        p12,
                        p13
                    }
                );

                //Pout<< "    23 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Missed edge:" << posEdge
                << abort(FatalError);
        }
    }
    else if (nPos == 1)
    {
        // Find the positive tet
        label i0 = -1;
        forAll(d, i)
        {
            if (d[i] > 0)
            {
                i0 = i;
                break;
            }
        }

        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);
        label i3 = d.fcIndex(i2);

        point p01(planeIntersection(d, tet, i0, i1));
        point p02(planeIntersection(d, tet, i0, i2));
        point p03(planeIntersection(d, tet, i0, i3));

        //Pout<< "Split 1pos tet " << tet << " d:" << d << " into" << nl;

        if (i0 == 0 || i0 == 2)
        {
            tetPoints t(tet[i0], p01, p02, p03);
            //Pout<< "    abovetet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 1, aboveTets i0==0 or 2");
            aboveOp(t);

            // Prism
            FixedList<point, 6> p
            (
                {
                    tet[i1],
                    tet[i3],
                    tet[i2],
                    p01,
                    p03,
                    p02
                }
            );

            //Pout<< "    belowprism:" << p << endl;
            decomposePrism(p, belowOp);
        }
        else
        {
            tetPoints t(p01, p02, p03, tet[i0]);
            //Pout<< "    abovetet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 1, aboveTets i0==1 or 3");
            aboveOp(t);

            // Prism
            FixedList<point, 6> p
            (
                {
                    tet[i3],
                    tet[i1],
                    tet[i2],
                    p03,
                    p01,
                    p02
                }
            );

            //Pout<< "    belowprism:" << p << endl;
            decomposePrism(p, belowOp);
        }
    }
    else    // nPos == 0
    {
        belowOp(tet);
    }
}


template<class Point, class PointRef>
template<class AboveTetOp, class BelowTetOp>
inline void Foam::tetrahedron<Point, PointRef>::sliceWithPlane
(
    const plane& pl,
    AboveTetOp& aboveOp,
    BelowTetOp& belowOp
) const
{
    tetSliceWithPlane(pl, tetPoints(a_, b_, c_, d_), aboveOp, belowOp);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    tetrahedron<Point, PointRef>& t
)
{
    is.readBegin("tetrahedron");
    is  >> t.a_ >> t.b_ >> t.c_ >> t.d_;
    is.readEnd("tetrahedron");

    is.check(FUNCTION_NAME);

    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tetrahedron<Point, PointRef>& t
)
{
    // Format like FixedList
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE
        << t.b_ << token::SPACE
        << t.c_ << token::SPACE
        << t.d_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
