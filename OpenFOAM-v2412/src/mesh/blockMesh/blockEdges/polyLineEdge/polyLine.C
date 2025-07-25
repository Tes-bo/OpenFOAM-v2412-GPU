/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "polyLine.H"
#include "SubList.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::polyLine::concat
(
    const point& p0,
    const pointField& intermediate,
    const point& p1
)
{
    auto tresult = tmp<pointField>::New(intermediate.size() + 2);
    auto& result = tresult.ref();

    // Intermediate points (knots)
    SubList<point>(result, intermediate.size(), 1) = intermediate;

    // Start/end points (knots)
    result.first() = p0;
    result.last() = p1;

    return tresult;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::polyLine::calcParam()
{
    lineLength_ = 0;
    param_.resize(points_.size());

    if (param_.size())
    {
        param_[0] = 0;

        for (label i=1; i < param_.size(); i++)
        {
            param_[i] = param_[i-1] + mag(points_[i] - points_[i-1]);
        }

        // Normalize on the interval 0-1
        lineLength_ = param_.last();
        for (label i=1; i < param_.size() - 1; i++)
        {
            param_[i] /= lineLength_;
        }
        param_.last() = 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLine::polyLine(const pointField& p, const bool)
:
    points_(p),
    lineLength_(0),
    param_()
{
    calcParam();
}


Foam::polyLine::polyLine
(
    const point& start,
    const pointField& intermediate,
    const point& end,
    const bool
)
:
    points_(polyLine::concat(start, intermediate, end)),
    lineLength_(0),
    param_()
{
    calcParam();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::polyLine::points() const noexcept
{
    return points_;
}


Foam::label Foam::polyLine::nSegments() const noexcept
{
    return points_.size()-1;
}


Foam::label Foam::polyLine::localParameter(scalar& lambda) const
{
    // Check endpoints
    if (lambda < SMALL)
    {
        lambda = 0;
        return 0;
    }
    else if (lambda > 1 - SMALL)
    {
        lambda = 1;
        return nSegments();
    }

    // Search table of cumulative distances to find which line-segment
    // we are on.
    // Check the upper bound.
    // Too small to bother with a binary search...

    label segment = 1;
    while (param_[segment] < lambda)
    {
        ++segment;
    }
    --segment;   // We want the corresponding lower bound

    // The local parameter [0-1] on this line segment
    lambda =
        (lambda - param_[segment])/(param_[segment+1] - param_[segment]);

    return segment;
}


Foam::point Foam::polyLine::position(const scalar mu) const
{
    // Check end-points
    if (mu < SMALL)
    {
        return points_.first();
    }
    else if (mu > 1 - SMALL)
    {
        return points_.last();
    }

    scalar lambda = mu;
    label segment = localParameter(lambda);
    return position(segment, lambda);
}


Foam::point Foam::polyLine::position
(
    const label segment,
    const scalar mu
) const
{
    // Out-of-bounds
    if (segment < 0)
    {
        return points_.first();
    }
    else if (segment > nSegments())
    {
        return points_.last();
    }

    const point& p0 = points_[segment];
    const point& p1 = points_[segment+1];

    // Special cases - no calculation needed
    if (mu <= 0.0)
    {
        return p0;
    }
    else if (mu >= 1.0)
    {
        return p1;
    }

    // Linear interpolation
    return points_[segment] + mu*(p1 - p0);
}


Foam::scalar Foam::polyLine::length() const noexcept
{
    return lineLength_;
}


// ************************************************************************* //
