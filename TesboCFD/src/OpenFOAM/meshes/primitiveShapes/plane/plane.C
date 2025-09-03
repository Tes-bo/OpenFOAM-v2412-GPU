/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "dictionary.H"
#include "plane.H"
#include "tensor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plane::makeUnitNormal
(
    const char * const caller,
    const bool notTest
)
{
    const scalar magNormal(Foam::mag(normal_));

    if (magNormal < VSMALL)
    {
        FatalErrorInFunction
            << "Plane normal has zero length.\nCalled from " << caller
            << abort(FatalError);
    }

    if (notTest)
    {
        normal_ /= magNormal;
    }
}


void Foam::plane::calcFromCoeffs
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d,
    const char * const caller
)
{
    if (mag(a) > VSMALL)
    {
        origin_ = vector((-d/a), 0, 0);
    }
    else if (mag(b) > VSMALL)
    {
        origin_ = vector(0, (-d/b), 0);
    }
    else if (mag(c) > VSMALL)
    {
        origin_ = vector(0, 0, (-d/c));
    }
    else
    {
        FatalErrorInFunction
            << "At least one plane coefficient must have a value"
            << abort(FatalError);
    }

    normal_ = vector(a, b, c);
    makeUnitNormal(caller);
}


void Foam::plane::calcFromEmbeddedPoints
(
    const point& point1,
    const point& point2,
    const point& point3,
    const char * const caller
)
{
    origin_ = (point1 + point2 + point3)/3;
    const vector line12 = point1 - point2;
    const vector line23 = point2 - point3;

    if
    (
        mag(line12) < VSMALL
     || mag(line23) < VSMALL
     || mag(point3-point1) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Bad points:" << point1 << ' ' << point2 << ' ' << point3
            << abort(FatalError);
    }

    normal_ = line12 ^ line23;

    makeUnitNormal(caller);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plane::plane(const vector& normalVector)
:
    normal_(normalVector),
    origin_(Zero)
{
    makeUnitNormal(FUNCTION_NAME);
}


Foam::plane::plane
(
    const point& originPoint,
    const vector& normalVector,
    const bool doNormalise
)
:
    normal_(normalVector),
    origin_(originPoint)
{
    makeUnitNormal(FUNCTION_NAME, doNormalise);
}


Foam::plane::plane(const UList<scalar>& coeffs)
{
    calcFromCoeffs
    (
        coeffs[0],
        coeffs[1],
        coeffs[2],
        coeffs[3],
        FUNCTION_NAME
    );
}


Foam::plane::plane(const FixedList<scalar,4>& coeffs)
{
    calcFromCoeffs
    (
        coeffs[0],
        coeffs[1],
        coeffs[2],
        coeffs[3],
        FUNCTION_NAME
    );
}


Foam::plane::plane(const point& a, const point& b, const point& c)
{
    calcFromEmbeddedPoints(a, b, c, FUNCTION_NAME);
}


Foam::plane::plane(const dictionary& dict)
:
    normal_(Zero),
    origin_(Zero)
{
    word planeType;
    dict.readIfPresent("planeType", planeType);

    if (planeType.empty())
    {
        const dictionary& coeffs = dict.optionalSubDict("pointAndNormalDict");

        origin_ = coeffs.get<point>("point");
        normal_ = coeffs.get<point>("normal");

        makeUnitNormal("point/normal");
    }
    else if (planeType == "pointAndNormal")
    {
        const dictionary& coeffs = dict.subDict("pointAndNormalDict");

        origin_ = coeffs.getCompat<point>("point", {{"basePoint", 1612}});
        normal_ = coeffs.getCompat<point>("normal", {{"normalVector", 1612}});

        makeUnitNormal("point/normal");
    }
    else if (planeType == "planeEquation")
    {
        const dictionary& subDict = dict.subDict("planeEquationDict");

        calcFromCoeffs
        (
            subDict.get<scalar>("a"),
            subDict.get<scalar>("b"),
            subDict.get<scalar>("c"),
            subDict.get<scalar>("d"),
            "planeEquation"  // caller name for makeUnitNormal
        );
    }
    else if (planeType == "embeddedPoints")
    {
        const dictionary& subDict = dict.subDict("embeddedPointsDict");

        calcFromEmbeddedPoints
        (
            subDict.get<point>("point1"),
            subDict.get<point>("point2"),
            subDict.get<point>("point3"),
            "embeddedPoints"  // caller name for makeUnitNormal
        );
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Invalid plane type: " << planeType << nl
            << "Valid options: (planeEquation embeddedPoints pointAndNormal)"
            << abort(FatalIOError);
    }
}


Foam::plane::plane(Istream& is)
:
    normal_(is),
    origin_(is)
{
    makeUnitNormal(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::scalar, 4> Foam::plane::planeCoeffs() const
{
    FixedList<scalar, 4> coeffs;

    const scalar magX = mag(normal_.x());
    const scalar magY = mag(normal_.y());
    const scalar magZ = mag(normal_.z());

    if (magX > magY)
    {
        if (magX > magZ)
        {
            coeffs[0] = 1;
            coeffs[1] = normal_.y()/normal_.x();
            coeffs[2] = normal_.z()/normal_.x();
        }
        else
        {
            coeffs[0] = normal_.x()/normal_.z();
            coeffs[1] = normal_.y()/normal_.z();
            coeffs[2] = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            coeffs[0] = normal_.x()/normal_.y();
            coeffs[1] = 1;
            coeffs[2] = normal_.z()/normal_.y();
        }
        else
        {
            coeffs[0] = normal_.x()/normal_.z();
            coeffs[1] = normal_.y()/normal_.z();
            coeffs[2] = 1;
        }
    }

    coeffs[3] =
    (
      - coeffs[0] * origin_.x()
      - coeffs[1] * origin_.y()
      - coeffs[2] * origin_.z()
    );

    return coeffs;
}


Foam::scalar Foam::plane::normalIntersect
(
    const point& pnt0,
    const vector& dir
) const
{
    const scalar denom = stabilise((dir & normal_), VSMALL);

    return ((origin_ - pnt0) & normal_)/denom;
}


Foam::plane::ray Foam::plane::planeIntersect(const plane& plane2) const
{
    // Mathworld plane-plane intersection. Assume there is a point on the
    // intersection line with z=0 and solve the two plane equations
    // for that (now 2x2 equation in x and y)
    // Better: use either z=0 or x=0 or y=0.

    const vector& n1 = this->normal();
    const vector& n2 = plane2.normal();

    const point& p1 = this->origin();
    const point& p2 = plane2.origin();

    const scalar n1p1 = n1 & p1;
    const scalar n2p2 = n2 & p2;

    const vector dir = n1 ^ n2;

    // Determine zeroed out direction (can be x,y or z) by looking at which
    // has the largest component in dir.
    const scalar magX = mag(dir.x());
    const scalar magY = mag(dir.y());
    const scalar magZ = mag(dir.z());

    direction iZero, i1, i2;

    if (magX > magY)
    {
        if (magX > magZ)
        {
            iZero = 0;
            i1 = 1;
            i2 = 2;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }
    else
    {
        if (magY > magZ)
        {
            iZero = 1;
            i1 = 2;
            i2 = 0;
        }
        else
        {
            iZero = 2;
            i1 = 0;
            i2 = 1;
        }
    }


    vector pt;

    pt[iZero] = 0;
    pt[i1] = (n2[i2]*n1p1 - n1[i2]*n2p2) / (n1[i1]*n2[i2] - n2[i1]*n1[i2]);
    pt[i2] = (n2[i1]*n1p1 - n1[i1]*n2p2) / (n1[i2]*n2[i1] - n1[i1]*n2[i2]);

    return ray(pt, dir);
}


Foam::point Foam::plane::planePlaneIntersect
(
    const plane& plane2,
    const plane& plane3
) const
{
    FixedList<scalar, 4> coeffs1(planeCoeffs());
    FixedList<scalar, 4> coeffs2(plane2.planeCoeffs());
    FixedList<scalar, 4> coeffs3(plane3.planeCoeffs());

    tensor a
    (
        coeffs1[0],coeffs1[1],coeffs1[2],
        coeffs2[0],coeffs2[1],coeffs2[2],
        coeffs3[0],coeffs3[1],coeffs3[2]
    );

    vector b(coeffs1[3],coeffs2[3],coeffs3[3]);

    return (inv(a) & (-b));
}


Foam::point Foam::plane::somePointInPlane(const scalar dist) const
{
    // ax + by + cz + d = 0
    const FixedList<scalar, 4> coeff(planeCoeffs());

    // Perturb the base-point
    point p = origin_ + point::uniform(dist);

    if (Foam::mag(coeff[2]) < SMALL)
    {
        if (Foam::mag(coeff[1]) < SMALL)
        {
            p[0] = -1.0*(coeff[1]*p[1] + coeff[2]*p[2] + coeff[3])/coeff[0];
        }
        else
        {
            p[1] = -1.0*(coeff[0]*p[0] + coeff[2]*p[2] + coeff[3])/coeff[1];
        }
    }
    else
    {
        p[2] = -1.0*(coeff[0]*p[0] + coeff[1]*p[1] + coeff[3])/coeff[2];
    }

    return p;
}


Foam::point Foam::plane::mirror(const point& p) const
{
    const vector mirroredPtDir = p - nearestPoint(p);

    if ((normal() & mirroredPtDir) > 0)
    {
        return p - 2.0*distance(p)*normal();
    }
    else
    {
        return p + 2.0*distance(p)*normal();
    }
}


void Foam::plane::writeDict(Ostream& os) const
{
    os.writeEntry("planeType", "pointAndNormal");

    os.beginBlock("pointAndNormalDict");

    os.writeEntry("point",  origin_);
    os.writeEntry("normal", normal_);

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const plane& pln)
{
    os << pln.normal() << token::SPACE << pln.origin();
    return os;
}


// ************************************************************************* //
