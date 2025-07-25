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

#include "IOstreams.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::edge::compare(const edge& a, const edge& b)
{
    return labelPair::compare(a, b);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::edge::edge()
:
    labelPair(-1, -1)
{}


inline Foam::edge::edge(const label from, const label to)
:
    labelPair(from, to)
{}


inline Foam::edge::edge(const labelPair& pair)
:
    labelPair(pair.first(), pair.second())
{}


inline Foam::edge::edge(const FixedList<label, 2>& list)
:
    labelPair(list.get<0>(), list.get<1>())
{}


inline Foam::edge::edge(const label from, const label to, const bool doSort)
:
    labelPair(from, to, doSort)
{}


inline Foam::edge::edge(const FixedList<label, 2>& list, const bool doSort)
:
    labelPair(list, doSort)
{}


inline Foam::edge::edge
(
    const labelUList& list,
    const FixedList<label, 2>& indices
)
:
    labelPair(list[indices.get<0>()], list[indices.get<1>()])
{}


inline Foam::edge::edge(Istream& is)
:
    labelPair(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::edge::minVertex() const noexcept
{
    return (first() < second() ? first() : second());
}


inline Foam::label Foam::edge::maxVertex() const noexcept
{
    return (second() < first() ? first() : second());
}


inline bool Foam::edge::good() const noexcept
{
    return (first() != second() && first() >= 0 && second() >= 0);
}


inline bool Foam::edge::contains(const label vertex) const noexcept
{
    // -1: always false
    return
    (
        vertex >= 0
     && (vertex == first() || vertex == second())
    );
}


inline Foam::label Foam::edge::which(const label vertex) const
{
    // -1: always false
    if (vertex >= 0)
    {
        if (vertex == first())
        {
            return 0;
        }
        if (vertex == second())
        {
            return 1;
        }
    }

    return -1;
}


inline bool Foam::edge::connected(const edge& other) const
{
    return (other.contains(first()) || other.contains(second()));
}


inline Foam::label Foam::edge::commonVertex(const edge& other) const
{
    if (other.contains(first()))
    {
        return first();
    }
    if (other.contains(second()))
    {
        return second();
    }

    // No shared vertex.
    return -1;
}


inline Foam::label Foam::edge::otherVertex(const label vertex) const
{
    if (vertex == first())
    {
        return second();
    }
    if (vertex == second())
    {
        return first();
    }

    // The given vertex is not on the edge in the first place.
    return -1;
}


inline Foam::label Foam::edge::collapse()
{
    // Cannot resize FixedList, so mark duplicates with '-1'
    // (the lower vertex is retained)
    // catch any '-1' (eg, if called multiple times)

    label n = 2;
    if (first() == second() || second() < 0)
    {
        second() = -1;
        --n;
    }
    if (first() < 0)
    {
        --n;
    }

    return n;
}


inline Foam::edge Foam::edge::reverseEdge() const
{
    return Foam::edge(second(), first());
}


inline void Foam::edge::clear()
{
    first()  = -1;
    second() = -1;
}


inline Foam::label Foam::edge::count() const
{
    label n = 2;
    if (first() == second() || second() < 0)
    {
        --n;
    }
    if (first() < 0)
    {
        --n;
    }

    return n;
}


inline bool Foam::edge::empty() const noexcept
{
    return (first() < 0 && second() < 0);
}


inline bool Foam::edge::insert(const label vertex)
{
    if (vertex < 0)
    {
        // Cannot insert invalid vertex labels (use direct assignment for that)
        return false;
    }

    if (first() < 0)
    {
        // Store at first, if not duplicate of second
        if (vertex != second())
        {
            first() = vertex;
            return true;
        }
    }
    else if (second() < 0)
    {
        // Store at second, if not duplicate of first
        if (vertex != first())
        {
            second() = vertex;
            return true;
        }
    }

    return false;
}


template<class InputIterator>
inline Foam::label Foam::edge::insert
(
    InputIterator begIter,
    InputIterator endIter
)
{
    // Available slots.
    // Don't use count() since it has special treatment for duplicates
    const int maxChange = ((first() < 0 ? 1 : 0) + (second() < 0 ? 1 : 0));

    int changed = 0;
    for (; changed < maxChange && begIter != endIter; ++begIter)
    {
        if (insert(*begIter))
        {
            ++changed;
        }
    }

    return changed;
}


inline Foam::label Foam::edge::insert(std::initializer_list<label> list)
{
    return insert(list.begin(), list.end());
}


template<unsigned N>
inline Foam::label Foam::edge::insert(const FixedList<label, N>& list)
{
    return insert(list.begin(), list.end());
}


inline Foam::label Foam::edge::insert(const labelUList& list)
{
    return insert(list.begin(), list.end());
}


inline Foam::label Foam::edge::erase(const label vertex)
{
    if (vertex < 0)
    {
        // Can never remove invalid point labels!
        return 0;
    }

    label n = 0;
    if (vertex == first())
    {
        first() = -1;
        ++n;
    }

    // Automatically handle duplicates, which should not have been there anyhow
    if (vertex == second())
    {
        second() = -1;
        ++n;
    }

    return n;
}


template<class InputIterator>
inline Foam::label Foam::edge::erase
(
    InputIterator begIter,
    InputIterator endIter
)
{
    // Occupied slots.
    // Don't use count() since it has special treatment for duplicates
    const int maxChange = ((first() >= 0 ? 1 : 0) + (second() >= 0 ? 1 : 0));

    int changed = 0;
    for (; changed < maxChange && begIter != endIter; ++begIter)
    {
        changed += erase(*begIter);
    }

    return changed;
}


inline Foam::label Foam::edge::erase(std::initializer_list<label> list)
{
    return erase(list.begin(), list.end());
}


template<unsigned N>
inline Foam::label Foam::edge::erase(const FixedList<label, N>& list)
{
    return erase(list.begin(), list.end());
}


inline Foam::label Foam::edge::erase(const labelUList& list)
{
    return erase(list.begin(), list.end());
}


// Geometric

inline Foam::point Foam::edge::centre(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (first() < 0 || second() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return 0.5*(pts[first()] + pts[second()]);
}


inline Foam::vector Foam::edge::vec(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (first() < 0 || second() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return pts[second()] - pts[first()];
}


inline Foam::vector Foam::edge::unitVec(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (first() < 0 || second() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    vector v = (pts[second()] - pts[first()]);
    (void) v.normalise(ROOTVSMALL);
    return v;
}


inline Foam::scalar Foam::edge::mag(const UList<point>& pts) const
{
    return pts[first()].dist(pts[second()]);
}


inline Foam::scalar Foam::edge::magSqr(const UList<point>& pts) const
{
    return pts[first()].distSqr(pts[second()]);
}


inline Foam::Pair<Foam::point>
Foam::edge::box(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (first() < 0 || second() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return linePointRef::box(pts[first()], pts[second()]);
}


inline Foam::linePointRef Foam::edge::line(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (first() < 0 || second() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return linePointRef(pts[first()], pts[second()]);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label& Foam::edge::operator[](const label i)
{
    #ifdef FULLDEBUG
    if (i < 0 || i > 1)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range [0,1]" << abort(FatalError);
    }
    #endif
    return (i ? second() : first());
}


inline const Foam::label& Foam::edge::operator[](const label i) const
{
    #ifdef FULLDEBUG
    if (i < 0 || i > 1)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range [0,1]" << abort(FatalError);
    }
    #endif
    return (i ? second() : first());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool Foam::operator==(const edge& a, const edge& b)
{
    return edge::compare(a,b) != 0;
}


inline bool Foam::operator!=(const edge& a, const edge& b)
{
    return edge::compare(a,b) == 0;
}


// ************************************************************************* //
