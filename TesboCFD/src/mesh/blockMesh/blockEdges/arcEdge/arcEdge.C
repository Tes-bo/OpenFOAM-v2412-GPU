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

#include "arcEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(arcEdge, 0);
    addToRunTimeSelectionTable(blockEdge, arcEdge, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockEdges::arcEdge::calcFromMidPoint
(
    const point& p1,
    const point& p3,
    const point& p2
)
{
    const vector a = p2 - p1;
    const vector b = p3 - p1;

    // Find centre of arcEdge
    const scalar asqr = a & a;
    const scalar bsqr = b & b;
    const scalar adotb = a & b;

    const scalar denom = asqr*bsqr - adotb*adotb;

    if (mag(denom) < ROOTVSMALL)
    {
        FatalErrorInFunction
            << denom
            << abort(FatalError);
    }

    const scalar fact = 0.5*(bsqr - adotb)/denom;

    const point centre = p1 + 0.5*a + fact*((a ^ b) ^ a);

    // Position vectors from centre
    const vector r1(p1 - centre);
    const vector r2(p2 - centre);
    const vector r3(p3 - centre);

    const scalar mag1(mag(r1));
    const scalar mag3(mag(r3));

    vector arcAxis(r1 ^ r3);

    // The radius from r1 and from r3 will be identical
    radius_ = mag(r3);


    // Determine the angle
    angle_ = acos((r1 & r3)/(mag1*mag3));

    // Check if the vectors define an exterior or an interior arcEdge
    if (((r1 ^ r2) & (r1 ^ r3)) < 0.0)
    {
        angle_ = constant::mathematical::twoPi - angle_;
    }

    if (angle_ <= constant::mathematical::pi)
    {
        if (mag(arcAxis)/(mag1*mag3) < 0.001)
        {
            arcAxis = r1 ^ r2;
        }
    }
    else
    {
        arcAxis = -arcAxis;
    }

    // Corresponding local cylindrical coordinate system
    cs_ = coordSystem::cylindrical(centre, arcAxis, r1);
}


void Foam::blockEdges::arcEdge::calcFromCentre
(
    const point& p1,
    const point& p3,
    const point& centre,
    bool adjustCentre,
    scalar rMultiplier
)
{
    // Position vectors from centre
    const vector r1(p1 - centre);
    const vector r3(p3 - centre);

    const scalar mag1(mag(r1));
    const scalar mag3(mag(r3));

    const vector chord(p3 - p1);

    const vector arcAxis(r1 ^ r3);

    // The average radius
    radius_ = 0.5*(mag1 + mag3);

    // The included angle
    angle_ = acos((r1 & r3)/(mag1*mag3));

    // TODO? check for 180 degrees (co-linear points)?

    bool needsAdjust = false;

    if (adjustCentre)
    {
        needsAdjust = !equal(mag1, mag3);

        if (!equal(rMultiplier, 1))
        {
            // The min radius is constrained by the chord,
            // otherwise bad things will happen.

            needsAdjust = true;
            radius_ *= rMultiplier;
            radius_ = max(radius_, (1.001*0.5*mag(chord)));
        }
    }

    if (needsAdjust)
    {
        // The centre is not equidistant to p1 and p3.
        // Use the chord and the arcAxis to determine the vector to
        // the midpoint of the chord and adjust the centre along this
        // line.

        const point newCentre =
        (
            (0.5 * (p3 + p1))                   // mid-chord point
          + sqrt(sqr(radius_) - 0.25 * magSqr(chord))
          * normalised(arcAxis ^ chord)         // mid-chord -> centre
        );

        //// Info<< nl << "Adjust centre. r1=" << mag1 << " r3=" << mag3
        ////     << " radius=" << radius_ << nl
        ////     << "angle=" << radToDeg(angle_) << ' '
        ////     << coordSystem::cylindrical(centre, arcAxis, r1) << nl;

        // Recalculate - do attempt to readjust
        calcFromCentre(p1, p3, newCentre, false);
    }
    else
    {
        // Corresponding local cylindrical coordinate system
        cs_ = coordSystem::cylindrical(centre, arcAxis, r1);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::arcEdge::arcEdge
(
    const pointField& points,
    const point& origin,
    const edge& fromTo
)
:
    blockEdge(points, fromTo),
    radius_(0),
    angle_(0),
    cs_()
{
    calcFromCentre(firstPoint(), lastPoint(), origin);
}


Foam::blockEdges::arcEdge::arcEdge
(
    const pointField& points,
    const edge& fromTo,
    const point& midPoint
)
:
    blockEdge(points, fromTo),
    radius_(0),
    angle_(0),
    cs_()
{
    calcFromMidPoint(firstPoint(), lastPoint(), midPoint);
}


Foam::blockEdges::arcEdge::arcEdge
(
    const pointField& points,
    const point& origin,
    const label from,
    const label to
)
:
    arcEdge(points, origin, edge(from,to))
{}


Foam::blockEdges::arcEdge::arcEdge
(
    const pointField& points,
    const label from,
    const label to,
    const point& midPoint
)
:
    arcEdge(points, edge(from,to), midPoint)
{}


Foam::blockEdges::arcEdge::arcEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces&,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    radius_(0),
    angle_(0),
    cs_()
{
    point p;

    token tok(is);
    if (tok.isWord())
    {
        // Can be
        //   - origin (0 0 0)
        //   - origin 1.2 (0 0 0)

        scalar rMultiplier = 1;

        is >> tok;
        if (tok.isNumber())
        {
            rMultiplier = tok.number();
        }
        else
        {
            is.putBack(tok);
        }

        is >> p;  // The origin (centre)

        calcFromCentre(firstPoint(), lastPoint(), p, true, rMultiplier);
    }
    else
    {
        is.putBack(tok);

        is >> p;  // A mid-point

        calcFromMidPoint(firstPoint(), lastPoint(), p);
    }

    if (debug)
    {
        Info<< "arc " << start_ << ' ' << end_ << ' '
            << position(0.5) << " origin " << cs_.origin() << " // ";
        cs_.rotation().write(Info);
        Info<< nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::arcEdge::position(const scalar lambda) const
{
    #ifdef FULLDEBUG
    if (lambda < -SMALL || lambda > 1 + SMALL)
    {
        InfoInFunction
            << "Limit parameter to [0-1] range: " << lambda << nl;
    }
    #endif

    if (lambda < SMALL)
    {
        return firstPoint();
    }
    else if (lambda >= 1 - SMALL)
    {
        return lastPoint();
    }

    return cs_.globalPosition(vector(radius_, (lambda*angle_), 0));
}


Foam::scalar Foam::blockEdges::arcEdge::length() const noexcept
{
    return (radius_*angle_);
}


// ************************************************************************* //
