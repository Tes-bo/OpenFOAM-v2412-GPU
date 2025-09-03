/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::face::face(const label sz)
:
    labelList(sz, -1)
{}


inline Foam::face::face(const labelUList& list)
:
    labelList(list)
{}


inline Foam::face::face(labelList&& list)
:
    labelList(std::move(list))
{}


inline Foam::face::face(std::initializer_list<label> list)
:
    labelList(list)
{}


template<unsigned N>
inline Foam::face::face(const FixedList<label, N>& list)
:
    labelList(list)
{}


inline Foam::face::face(const labelUList& list, const labelUList& indices)
:
    labelList(list, indices)
{}


template<unsigned N>
inline Foam::face::face
(
    const labelUList& list,
    const FixedList<label, N>& indices
)
:
    labelList(list, indices)
{}


inline Foam::face::face(Istream& is)
:
    labelList(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::face::points(const UList<point>& pts) const
{
    // There are as many points as there are labels for them
    pointField p(size());

    auto iter = p.begin();

    for (const label pointi : *this)
    {
        *iter = pts[pointi];
        ++iter;
    }

    return p;
}


inline Foam::vector Foam::face::unitNormal(const UList<point>& p) const
{
    vector n(areaNormal(p));
    (void) n.normalise(ROOTVSMALL);
    return n;
}


inline Foam::scalar Foam::face::mag(const UList<point>& p) const
{
    return ::Foam::mag(areaNormal(p));
}


inline Foam::Pair<Foam::point>
Foam::face::box(const UList<point>& pts) const
{
    Pair<point> bb(point::rootMax, point::rootMin);

    for (const label pointi : *this)
    {
        bb.first()  = min(bb.first(),  pts[pointi]);
        bb.second() = max(bb.second(), pts[pointi]);
    }

    return bb;
}


inline Foam::label Foam::face::nEdges() const noexcept
{
    // for a closed polygon a number of edges is the same as number of points
    return size();
}


inline Foam::edge Foam::face::edge(const label edgei) const
{
    return Foam::edge(thisLabel(edgei), nextLabel(edgei));
}


inline Foam::vector Foam::face::edge
(
    const label edgei,
    const UList<point>& pts
) const
{
    return vector(pts[nextLabel(edgei)] - pts[thisLabel(edgei)]);
}


inline Foam::edge Foam::face::rcEdge(const label edgei) const
{
    // Edge 0 (forward and reverse) always starts at [0]
    // for consistency with face flipping
    const label pointi = edgei ? (nEdges() - edgei) : 0;
    return Foam::edge(thisLabel(pointi), prevLabel(pointi));
}


inline Foam::vector Foam::face::rcEdge
(
    const label edgei,
    const UList<point>& pts
) const
{
    // Edge 0 (forward and reverse) always starts at [0]
    // for consistency with face flipping
    const label pointi = edgei ? (nEdges() - edgei) : 0;
    return vector(pts[prevLabel(pointi)] - pts[thisLabel(pointi)]);
}


inline Foam::label Foam::face::which(const label vertex) const
{
    return labelList::find(vertex);
}


inline Foam::label Foam::face::thisLabel(const label i) const
{
    return labelList::operator[](i);
}


inline Foam::label Foam::face::nextLabel(const label i) const
{
    return labelList::fcValue(i);
}


inline Foam::label Foam::face::prevLabel(const label i) const
{
    return labelList::rcValue(i);
}


inline Foam::label Foam::face::nTriangles() const noexcept
{
    return labelList::size() - 2;
}


inline bool Foam::face::connected(const labelUList& other) const
{
    for (const label pointi : *this)
    {
        if (other.contains(pointi))
        {
            return true;
        }
    }
    return false;
}


template<unsigned N>
inline bool Foam::face::connected(const FixedList<label, N>& other) const
{
    for (const label pointi : *this)
    {
        if (other.contains(pointi))
        {
            return true;
        }
    }
    return false;
}


inline bool Foam::face::contains(const Foam::edge& e) const
{
    // or (find(e) >= 0)
    return (edgeDirection(e) != 0);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::face::operator+=(const label vertexOffset)
{
    if (vertexOffset)
    {
        for (label& vrt : static_cast<labelList&>(*this))
        {
            vrt += vertexOffset;
        }
    }
}


// * * * * * * * * * * * * * * * Global Operators   * * * * * * * * * * * * * //

inline bool Foam::operator==(const face& a, const face& b)
{
    return face::compare(a,b) != 0;
}

inline bool Foam::operator!=(const face& a, const face& b)
{
    return face::compare(a,b) == 0;
}


// ************************************************************************* //
