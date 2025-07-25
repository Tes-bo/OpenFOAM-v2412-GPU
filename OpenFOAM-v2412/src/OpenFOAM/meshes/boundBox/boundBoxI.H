/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "boundBox.H"

// * * * * * * * * * * * * * Geometrical Information * * * * * * * * * * * * //

namespace Foam
{

// Box corners as per hex cellmodel

template<>
inline Foam::point Foam::boundBox::hexCorner<0>() const
{
    // == octCorner<0>()
    return min_;
}

template<>
inline Foam::point Foam::boundBox::hexCorner<1>() const
{
    // == octCorner<1>()
    return point(max_.x(), min_.y(), min_.z());
}

template<>
inline Foam::point Foam::boundBox::hexCorner<2>() const
{
    // == octCorner<3>()
    return point(max_.x(), max_.y(), min_.z());
}

template<>
inline Foam::point Foam::boundBox::hexCorner<3>() const
{
    // == octCorner<2>()
    return point(min_.x(), max_.y(), min_.z());
}

template<>
inline Foam::point Foam::boundBox::hexCorner<4>() const
{
    // == octCorner<4>()
    return point(min_.x(), min_.y(), max_.z());
}

template<>
inline Foam::point Foam::boundBox::hexCorner<5>() const
{
    // == octCorner<5>()
    return point(max_.x(), min_.y(), max_.z());
}

template<>
inline Foam::point Foam::boundBox::hexCorner<6>() const
{
    // == octCorner<7>()
    return max_;
}

template<>
inline Foam::point Foam::boundBox::hexCorner<7>() const
{
    // == octCorner<6>()
    return point(min_.x(), max_.y(), max_.z());
}

}  // End namespace Foam


// Non-specialized version is compile-time disabled
template<Foam::direction CornerNumber>
inline Foam::point Foam::boundBox::hexCorner() const
{
    static_assert(CornerNumber < 8, "Corner index [0..7]");
    return point();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::boundBox::boundBox()
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{}


inline Foam::boundBox::boundBox(const Foam::zero_one)
:
    min_(point::zero),
    max_(point::one)
{}


inline Foam::boundBox::boundBox(const point& p)
:
    min_(p),
    max_(p)
{}


inline Foam::boundBox::boundBox(const point& min, const point& max)
:
    min_(min),
    max_(max)
{}


inline Foam::boundBox::boundBox(const Pair<point>& bb)
:
    min_(bb.first()),
    max_(bb.second())
{}


inline Foam::boundBox::boundBox(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::boundBox::empty() const
{
    // Is empty/invalid if any component has (max < min), ie, !(min <= max)
    return
    (
        (max_.x() < min_.x())
     || (max_.y() < min_.y())
     || (max_.z() < min_.z())
    );
}


inline bool Foam::boundBox::good() const
{
    return !empty();
}


inline const Foam::point& Foam::boundBox::min() const noexcept
{
    return min_;
}


inline const Foam::point& Foam::boundBox::max() const noexcept
{
    return max_;
}


inline Foam::point& Foam::boundBox::min() noexcept
{
    return min_;
}


inline Foam::point& Foam::boundBox::max() noexcept
{
    return max_;
}


inline Foam::point Foam::boundBox::centre() const
{
    return 0.5 * (min_ + max_);
}


inline Foam::vector Foam::boundBox::span() const
{
    return (max_ - min_);
}


inline Foam::scalar Foam::boundBox::mag() const
{
    return min_.dist(max_);
}


inline Foam::scalar Foam::boundBox::magSqr() const
{
    return min_.distSqr(max_);
}


inline Foam::scalar Foam::boundBox::volume() const
{
    return cmptProduct(span());
}


inline Foam::scalar Foam::boundBox::minDim() const
{
    return cmptMin(span());
}


inline Foam::scalar Foam::boundBox::maxDim() const
{
    return cmptMax(span());
}


inline Foam::scalar Foam::boundBox::avgDim() const
{
    return cmptAv(span());
}


inline Foam::direction Foam::boundBox::minDir() const
{
    direction cmpt = 0;

    scalar best = ROOTVGREAT;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist < best && dist > 0)
        {
            best = dist;
            cmpt = dir;
        }
    }

    return cmpt;
}


inline Foam::direction Foam::boundBox::maxDir() const
{
    direction cmpt = 0;

    scalar best = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist > best)
        {
            best = dist;
            cmpt = dir;
        }
    }

    return cmpt;
}


inline int Foam::boundBox::nDim() const
{
    int ngood = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar dist = (max_[dir] - min_[dir]);
        if (dist < 0)
        {
            return -1;
        }
        else if (dist > 0)
        {
            ++ngood;
        }
    }

    return ngood;
}


inline void Foam::boundBox::reset()
{
    min_ = invertedBox.min();
    max_ = invertedBox.max();
}


inline void Foam::boundBox::reset(const Foam::zero_one)
{
    min_ = point::zero;
    max_ = point::one;
}


inline void Foam::boundBox::reset(const point& pt)
{
    min_ = pt;
    max_ = pt;
}


inline void Foam::boundBox::reset(const point& min, const point& max)
{
    min_ = min;
    max_ = max;
}


inline void Foam::boundBox::add(const boundBox& bb)
{
    min_ = ::Foam::min(min_, bb.min_);
    max_ = ::Foam::max(max_, bb.max_);
}


inline void Foam::boundBox::add(const point& pt)
{
    min_ = ::Foam::min(min_, pt);
    max_ = ::Foam::max(max_, pt);
}


inline void Foam::boundBox::add(const point& pt0, const point& pt1)
{
    add(pt0);
    add(pt1);
}


inline void Foam::boundBox::add(const Pair<point>& points)
{
    add(points.first());
    add(points.second());
}


inline void Foam::boundBox::add(const UList<point>& points)
{
    for (const point& p : points)
    {
        add(p);
    }
}


inline void Foam::boundBox::add(const tmp<pointField>& tpoints)
{
    add(tpoints());
    tpoints.clear();
}


inline void Foam::boundBox::grow(const scalar delta)
{
    min_.x() -= delta; min_.y() -= delta; min_.z() -= delta;
    max_.x() += delta; max_.y() += delta; max_.z() += delta;
}


inline void Foam::boundBox::grow(const vector& delta)
{
    min_ -= delta;
    max_ += delta;
}


inline void Foam::boundBox::inflate(const scalar factor)
{
    grow(factor*mag());
}


inline bool Foam::boundBox::box_box_overlaps
(
    const point& minA, const point& maxA,   // boxA
    const point& minB, const point& maxB    // boxB
)
{
    return
    (
        minA.x() <= maxB.x() && minB.x() <= maxA.x()
     && minA.y() <= maxB.y() && minB.y() <= maxA.y()
     && minA.z() <= maxB.z() && minB.z() <= maxA.z()
    );
}


inline bool Foam::boundBox::box_sphere_overlaps
(
    const point& corner0,
    const point& corner1,
    const point& centre,
    const scalar radiusSqr
)
{
    // Find out where centre is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar d0 = corner0[dir] - centre[dir];
        const scalar d1 = corner1[dir] - centre[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // centre inside both extrema. This component does not add any
            // distance.
        }
        else
        {
            distSqr += ::Foam::min(Foam::magSqr(d0), Foam::magSqr(d1));

            if (distSqr > radiusSqr)
            {
                return false;
            }
        }
    }

    return true;
}


inline bool Foam::boundBox::overlaps(const boundBox& bb) const
{
    return box_box_overlaps(min_, max_, bb.min(), bb.max());
}


inline bool Foam::boundBox::overlaps
(
    const point& centre,
    const scalar radiusSqr
) const
{
    return box_sphere_overlaps(min_, max_, centre, radiusSqr);
}


inline bool Foam::boundBox::contains(const point& pt) const
{
    return
    (
        min_.x() <= pt.x() && pt.x() <= max_.x()
     && min_.y() <= pt.y() && pt.y() <= max_.y()
     && min_.z() <= pt.z() && pt.z() <= max_.z()
    );
}


inline bool Foam::boundBox::contains(const boundBox& bb) const
{
    return contains(bb.min()) && contains(bb.max());
}


inline bool Foam::boundBox::containsInside(const point& pt) const
{
    return
    (
        min_.x() < pt.x() && pt.x() < max_.x()
     && min_.y() < pt.y() && pt.y() < max_.y()
     && min_.z() < pt.z() && pt.z() < max_.z()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
