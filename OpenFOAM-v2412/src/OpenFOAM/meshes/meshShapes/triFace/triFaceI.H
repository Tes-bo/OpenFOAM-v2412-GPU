/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "IOstreams.H"
#include "face.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::triFace::compare(const triFace& a, const triFace& b)
{
    if
    (
        (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
     || (a[0] == b[1] && a[1] == b[2] && a[2] == b[0])
     || (a[0] == b[2] && a[1] == b[0] && a[2] == b[1])
    )
    {
        // identical
        return 1;
    }
    else if
    (
        (a[0] == b[2] && a[1] == b[1] && a[2] == b[0])
     || (a[0] == b[1] && a[1] == b[0] && a[2] == b[2])
     || (a[0] == b[0] && a[1] == b[2] && a[2] == b[1])
    )
    {
        // same face, but reversed orientation
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::triFace::triFace()
:
    FixedList<label, 3>(-1)
{}


inline Foam::triFace::triFace
(
    const label p0,
    const label p1,
    const label p2
) noexcept
{
    a() = p0;
    b() = p1;
    c() = p2;
}


inline Foam::triFace::triFace(std::initializer_list<label> list)
:
    FixedList<label, 3>(list)
{}


inline Foam::triFace::triFace(const labelUList& list)
:
    FixedList<label, 3>(list)
{}


inline Foam::triFace::triFace
(
    const labelUList& list,
    const FixedList<label, 3>& indices
)
:
    FixedList<label, 3>(list, indices)
{}


inline Foam::triFace::triFace(Istream& is)
:
    FixedList<label, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::triFace::good() const noexcept
{
    return
    (
        a() >= 0 && a() != b()
     && b() >= 0 && b() != c()
     && c() >= 0 && c() != a()
    );
}


inline Foam::label Foam::triFace::collapse()
{
    // Cannot resize FixedList, so mark duplicates with '-1'
    // (the lower vertex is retained)
    // catch any '-1' (eg, if called multiple times)

    label n = 3;
    if (operator[](0) == operator[](1) || operator[](1) == -1)
    {
        operator[](1) = -1;
        n--;
    }
    else if (operator[](1) == operator[](2) || operator[](2) == -1)
    {
        operator[](2) = -1;
        n--;
    }
    if (operator[](0) == operator[](2))
    {
        operator[](2) = -1;
        n--;
    }

    return n;
}


inline void Foam::triFace::flip()
{
    std::swap(get<1>(), get<2>());
}


inline Foam::pointField Foam::triFace::points(const UList<point>& pts) const
{
    pointField p(3);

    p[0] = pts[a()];
    p[1] = pts[b()];
    p[2] = pts[c()];

    return p;
}


inline Foam::face Foam::triFace::triFaceFace() const
{
    return Foam::face(*this);
}


inline Foam::triPointRef Foam::triFace::tri(const UList<point>& points) const
{
    return triPointRef(points[a()], points[b()], points[c()]);
}


inline Foam::point Foam::triFace::centre(const UList<point>& points) const
{
    return triPointRef::centre(points[a()], points[b()], points[c()]);
}


inline Foam::vector Foam::triFace::areaNormal(const UList<point>& points) const
{
    return triPointRef::areaNormal(points[a()], points[b()], points[c()]);
}


inline Foam::vector Foam::triFace::unitNormal(const UList<point>& points) const
{
    return triPointRef::unitNormal(points[a()], points[b()], points[c()]);
}


inline Foam::scalar Foam::triFace::mag(const UList<point>& points) const
{
    return ::Foam::mag(areaNormal(points));
}


inline Foam::Pair<Foam::point>
Foam::triFace::box(const UList<point>& points) const
{
    return triPointRef::box(points[a()], points[b()], points[c()]);
}


inline Foam::label Foam::triFace::nTriangles() const noexcept
{
    return 1;
}


inline Foam::triFace Foam::triFace::reverseFace() const
{
    // The starting points of the original and reverse face are identical.
    return triFace(a(), c(), b());
}


inline Foam::label Foam::triFace::which(const label vertex) const
{
    return FixedList<label, 3>::find(vertex);
}


inline Foam::label Foam::triFace::thisLabel(const label i) const
{
    return operator[](i);
}


inline Foam::label Foam::triFace::nextLabel(const label i) const
{
    return operator[]((i == 2 ? 0 : i+1));
}


inline Foam::label Foam::triFace::prevLabel(const label i) const
{
    return operator[]((i ? i-1 : 2));
}


inline Foam::scalar Foam::triFace::sweptVol
(
    const UList<point>& opts,
    const UList<point>& npts
) const
{
    return (1.0/6.0)*
    (
        (
            (npts[operator[](0)] - opts[operator[](0)])
          & (
                (opts[operator[](1)] - opts[operator[](0)])
              ^ (opts[operator[](2)] - opts[operator[](0)])
            )
        )
      + (
            (npts[operator[](1)] - opts[operator[](1)])
          & (
                (opts[operator[](2)] - opts[operator[](1)])
              ^ (npts[operator[](0)] - opts[operator[](1)])
            )
        )
      + (
            (opts[operator[](2)] - npts[operator[](2)])
          & (
                (npts[operator[](1)] - npts[operator[](2)])
              ^ (npts[operator[](0)] - npts[operator[](2)])
            )
        )
    );
}


Foam::tensor Foam::triFace::inertia
(
    const UList<point>& points,
    const point& refPt,
    scalar density
) const
{
    // a triangle, do a direct calculation
    return this->tri(points).inertia(refPt, density);
}


inline Foam::pointHit Foam::triFace::ray
(
    const point& p,
    const vector& q,
    const UList<point>& points,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    return this->tri(points).ray(p, q, alg, dir);
}



inline Foam::pointHit Foam::triFace::intersection
(
    const point& p,
    const vector& q,
    const UList<point>& points,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    return this->tri(points).intersection(p, q, alg, tol);
}


inline Foam::pointHit Foam::triFace::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const UList<point>& points,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    return intersection(p, q, points, alg, tol);
}


inline Foam::pointHit Foam::triFace::nearestPoint
(
    const point& p,
    const UList<point>& points
) const
{
    return this->tri(points).nearestPoint(p);
}


inline Foam::pointHit Foam::triFace::nearestPointClassify
(
    const point& p,
    const UList<point>& points,
    label& nearType,
    label& nearLabel
) const
{
    return this->tri(points).nearestPointClassify(p, nearType, nearLabel);
}


inline int Foam::triFace::sign
(
    const point& p,
    const UList<point>& points,
    const scalar tol
) const
{
    return this->tri(points).sign(p, tol);
}


inline Foam::label Foam::triFace::nEdges() const noexcept
{
    return 3;
}


inline Foam::edge Foam::triFace::edge(const label edgei) const
{
    return Foam::edge(thisLabel(edgei), nextLabel(edgei));
}


inline Foam::vector Foam::triFace::edge
(
    const label edgei,
    const UList<point>& pts
) const
{
    return vector(pts[nextLabel(edgei)] - pts[thisLabel(edgei)]);
}


inline Foam::edge Foam::triFace::rcEdge(const label edgei) const
{
    // Edge 0 (forward and reverse) always starts at [0]
    // for consistency with face flipping
    const label pointi = edgei ? (3 - edgei) : 0;
    return Foam::edge(thisLabel(pointi), prevLabel(pointi));
}


inline Foam::vector Foam::triFace::rcEdge
(
    const label edgei,
    const UList<point>& pts
) const
{
    // Edge 0 (forward and reverse) always starts at [0]
    // for consistency with face flipping
    const label pointi = edgei ? (3 - edgei) : 0;
    return vector(pts[prevLabel(pointi)] - pts[thisLabel(pointi)]);
}


inline Foam::edgeList Foam::triFace::edges() const
{
    edgeList theEdges(3);

    theEdges[0].first()  = a();
    theEdges[0].second() = b();

    theEdges[1].first()  = b();
    theEdges[1].second() = c();

    theEdges[2].first()  = c();
    theEdges[2].second() = a();

    return theEdges;
}


inline Foam::edgeList Foam::triFace::rcEdges() const
{
    edgeList theEdges(3);

    theEdges[0].first()  = a();
    theEdges[0].second() = c();

    theEdges[1].first()  = c();
    theEdges[1].second() = b();

    theEdges[2].first()  = b();
    theEdges[2].second() = a();

    return theEdges;
}


inline int Foam::triFace::edgeDirection(const Foam::edge& e) const
{
    if (e.first() == a())
    {
        if (e.second() == b()) return 1;    // Forward edge 0 (encoded 1)
        if (e.second() == c()) return -1;   // Reverse edge 2 (encoded -3)
    }
    if (e.first() == b())
    {
        if (e.second() == c()) return 1;    // Forward edge 1 (encoded 2)
        if (e.second() == a()) return -1;   // Reverse edge 0 (encoded -1)
    }
    if (e.first() == c())
    {
        if (e.second() == a()) return 1;    // Forward edge 2 (encoded 3)
        if (e.second() == b()) return -1;   // Reverse edge 1 (encoded -2)
    }

    return 0;  // Not found
}


inline Foam::label Foam::triFace::find(const Foam::edge& e) const
{
    if (e.first() == a())
    {
        if (e.second() == b()) return 0;    // Forward edge 0
        if (e.second() == c()) return 2;    // Reverse edge 2
    }
    if (e.first() == b())
    {
        if (e.second() == c()) return 1;    // Forward edge 1
        if (e.second() == a()) return 0;    // Reverse edge 0
    }
    if (e.first() == c())
    {
        if (e.second() == a()) return 2;    // Forward edge 2
        if (e.second() == b()) return 1;    // Reverse edge 1
    }

    return -1;  // Not found
}


inline bool Foam::triFace::contains(const Foam::edge& e) const
{
    // or (find(e) >= 0)
    return (edgeDirection(e) != 0);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::triFace::operator+=(const label vertexOffset)
{
    if (vertexOffset)
    {
        (*this)[0] += vertexOffset;
        (*this)[1] += vertexOffset;
        (*this)[2] += vertexOffset;
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const triFace& a, const triFace& b)
{
    return triFace::compare(a,b) != 0;
}


inline bool Foam::operator!=(const triFace& a, const triFace& b)
{
    return triFace::compare(a,b) == 0;
}


// ************************************************************************* //
