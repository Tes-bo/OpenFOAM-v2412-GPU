/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "CatmullRomSpline.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::CatmullRomSpline::derivative
(
    const label segment,
    const scalar mu
) const
{
    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];

    // determine the end points
    point e0;
    point e1;

    if (segment == 0)
    {
        // end: simple reflection
        e0 = 2*p0 - p1;
    }
    else
    {
        e0 = points()[segment-1];
    }

    if (segment+1 == nSegments())
    {
        // end: simple reflection
        e1 = 2*p1 - p0;
    }
    else
    {
        e1 = points()[segment+2];
    }
    const point derivativePoint
    (
        0.5 *
        (
            (-e0 + p1)
          + mu *
            (
                2 * (2*e0 - 5*p0 + 4*p1 - e1)
              + mu * 3 * (-e0 + 3*p0 - 3*p1 + e1)
            )
        )
    );
    return mag(derivativePoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CatmullRomSpline::CatmullRomSpline
(
    const pointField& knots,
    const bool closed
)
:
    polyLine(knots, closed)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::CatmullRomSpline::position(const scalar mu) const
{
    // endpoints
    if (mu < SMALL)
    {
        return points().first();
    }
    else if (mu > 1 - SMALL)
    {
        return points().last();
    }

    scalar lambda = mu;
    label segment = localParameter(lambda);
    return position(segment, lambda);
}


Foam::point Foam::CatmullRomSpline::position
(
    const label segment,
    const scalar mu
) const
{
    // out-of-bounds
    if (segment < 0)
    {
        return points().first();
    }
    else if (segment > nSegments())
    {
        return points().last();
    }

    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];

    // special cases - no calculation needed
    if (mu <= 0.0)
    {
        return p0;
    }
    else if (mu >= 1.0)
    {
        return p1;
    }


    // determine the end points
    point e0;
    point e1;

    if (segment == 0)
    {
        // end: simple reflection
        e0 = 2*p0 - p1;
    }
    else
    {
        e0 = points()[segment-1];
    }

    if (segment+1 == nSegments())
    {
        // end: simple reflection
        e1 = 2*p1 - p0;
    }
    else
    {
        e1 = points()[segment+2];
    }


    return
        0.5 *
        (
            (2*p0)
          + mu *
            (
                (-e0 + p1)
              + mu *
                (
                    (2*e0 - 5*p0 + 4*p1 - e1)
                  + mu*(-e0 + 3*p0 - 3*p1 + e1)
                )
            )
        );
}


Foam::scalar Foam::CatmullRomSpline::length() const
{
    const solveScalar xi[5]=
    {
        -0.9061798459386639927976,
        -0.5384693101056830910363,
        0,
        0.5384693101056830910363,
        0.9061798459386639927976
    };
    const solveScalar wi[5]=
    {
        0.2369268850561890875143,
        0.4786286704993664680413,
        0.5688888888888888888889,
        0.4786286704993664680413,
        0.2369268850561890875143
    };
    scalar sum=0;
    for (label segment=0;segment<nSegments();segment++)
    {
        for (int i=0;i<5;i++)
        {
            sum+=wi[i]*derivative(segment,(xi[i]+1.0)/2.0)/2.0;
        }
    }

    return sum;
}


// ************************************************************************* //
