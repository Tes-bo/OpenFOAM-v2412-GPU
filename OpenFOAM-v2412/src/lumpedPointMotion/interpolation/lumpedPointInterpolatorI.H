/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

inline Foam::lumpedPointInterpolator::lumpedPointInterpolator()
:
    lumpedPointInterpolator(-1)
{}


inline Foam::lumpedPointInterpolator::lumpedPointInterpolator(const label id)
:
    nearest_(id),
    next1_(id),
    next2_(id),
    weight1_(Zero),
    weight2_(Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::lumpedPointInterpolator::good() const
{
    return nearest_ != -1;
}


inline Foam::label Foam::lumpedPointInterpolator::nearest() const
{
    return nearest_;
}


inline Foam::label Foam::lumpedPointInterpolator::next() const
{
    return next1_;
}


inline Foam::label Foam::lumpedPointInterpolator::next1() const
{
    return next1_;
}


inline Foam::label Foam::lumpedPointInterpolator::next2() const
{
    return next2_;
}


inline Foam::scalar Foam::lumpedPointInterpolator::weight0() const
{
    return scalar(1) - weight1_ - weight2_;
}


inline Foam::scalar Foam::lumpedPointInterpolator::weight1() const
{
    return weight1_;
}


inline Foam::scalar Foam::lumpedPointInterpolator::weight2() const
{
    return weight2_;
}


inline void Foam::lumpedPointInterpolator::nearest(const label id)
{
    nearest_ = id;
    next(id, Zero);
}


inline void Foam::lumpedPointInterpolator::next
(
    const label id,
    scalar weight
)
{
    if (weight < scalar(0))
    {
        weight = 0;
    }
    else if (weight > scalar(1))
    {
        weight = 1;
    }

    next1_ = id;
    next2_ = id;
    weight1_ = weight;
    weight2_ = Zero;
}


inline void Foam::lumpedPointInterpolator::next
(
    const label id,
    scalar weight,
    const label position
)
{
    if (weight < scalar(0))
    {
        weight = 0;
    }
    else if (weight > scalar(1))
    {
        weight = 1;
    }

    if (!position)
    {
        next1_ = id;
        weight1_ = weight;
    }
    else
    {
        next2_ = id;
        weight2_ = weight;
    }
}


inline void Foam::lumpedPointInterpolator::set
(
    const triFace& ids,
    const barycentric2D& weights
)
{
    nearest_ = ids[0];
    next1_ = ids[1];
    next2_ = ids[2];

    weight1_ = weights[1];
    weight2_ = weights[2];
}


template<class T>
inline T Foam::lumpedPointInterpolator::interpolate(const UList<T>& input) const
{
    if (nearest_ == -1)
    {
        return Zero;
    }
    else if (next1_ == -1 || next1_ == nearest_)
    {
        return input[nearest_];
    }

    return
    (
        input[nearest_] * (1-weight1_-weight2_)
      + input[next1_] * (weight1_)
      + input[next2_] * (weight2_)
    );
}


// ************************************************************************* //
