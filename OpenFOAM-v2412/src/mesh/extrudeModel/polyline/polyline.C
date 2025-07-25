/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Ivor Clifford/Paul Scherrer Institut
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

#include "polyline.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyline, 0);

addToRunTimeSelectionTable(extrudeModel, polyline, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyline::polyline(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    geometry_(0),
    vertices_(coeffDict_.lookup("vertices")),
    segments_
    (
        coeffDict_.lookup("edges"),
        blockEdge::iNew(coeffDict_, geometry_, vertices_)
    ),
    x_(segments_.size() + 1),
    y_(segments_.size() + 1),
    relTol_(coeffDict_.getOrDefault<scalar>("toleranceCheck", SMALL))
{
    // Check continuity and smoothness of the supplied polyline
    for (label i=1; i < segments_.size(); ++i)
    {
        // Check continuity
        const vector x0 = segments_[i-1].position(1);
        const vector x1 = segments_[i].position(0);

        if (mag(x1-x0) > SMALL)
        {
            FatalErrorInFunction()
                << "Supplied polyline is not continuous." << endl
                << Foam::abort(FatalError);
        }

        // Check smoothness
        const vector v0 =
            normalised
            (
                segments_[i-1].position(1)
              - segments_[i-1].position(1-DELTA)
            );

        const vector v1 =
            normalised
            (
                segments_[i].position(DELTA)
              - segments_[i].position(0)
            );

        if ((v1 & v0) < (1 - relTol_))
        {
            FatalErrorInFunction()
                << "Supplied polyline is not smooth." << endl
                << Foam::abort(FatalError);
        }
    }

    // Calculate cumulative length along polyline
    x_[0] = 0;
    y_[0] = 0;
    scalar totalLength = 0;
    forAll(segments_, i)
    {
        totalLength += segments_[i].length();
        x_[i+1] = totalLength;
        y_[i+1] = i+1;
    }

    // Normalise cumulative length (0 <= x <= 1)
    x_ /= totalLength;

    // Position vector and direction at start of polyline
    positionAndDirection(0, p0_, n0_);

    if (debug)
    {
        Info
            << tab << "Polyline start: " << p0_ << nl
            << tab << "Polyline normal at start: " << n0_ << nl
            << tab << "Polyline end: "
            << segments_.last().position(1) << nl
            << tab << "Total length: " << totalLength << endl;
    }
}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point polyline::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    // Offset between supplied point and origin of polyline
    vector dp = (surfacePoint - p0_);

    // If this is the first layer, check whether the start of the
    // polyline seems to lie on the surface
    if (layer == 0)
    {
        if (mag((dp/mag(dp)) & n0_) > relTol_)
        {
            WarningInFunction()
                << "The starting point of the polyline does not appear "
                << "to lie of the supplied surface. Apparent absolute "
                << "misalignment is " << (dp & n0_) << endl;
        }
    }

    // Position and direction vector at end of layer
    vector p;
    vector n;
    positionAndDirection(sumThickness(layer), p, n);

    // Angle between normal vector and normal at origin
    scalar cosTheta = (n & n0_);

    // Rotate point to align with current normal vector
    if (cosTheta < (1-SMALL))
    {
        const vector axis = normalised(n0_ ^ n);

        dp = quaternion(axis, cosTheta, true).transform(dp);
    }

    return p + dp;
}


void polyline::positionAndDirection
(
    const scalar lambda,
    vector& p,
    vector& n
) const
{
    // Find associated segment and position for supplied lambda
    scalar y = interpolateXY(lambda, x_, y_);
    int i = floor(y);
    scalar s = y - i;
    if (i > segments_.size()-1)
    {
        i = segments_.size()-1;
        s = 1.0;
    }

    // Position vector
    p = segments_[i].position(s);

    // Normal vector at current position
    // Estimated normal vector using numerical differencing since
    // blockEdge doesn't include a normal function

    n = normalised
    (
        segments_[i].position(min(s + DELTA, 1))
      - segments_[i].position(max(s - DELTA, 0))
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
