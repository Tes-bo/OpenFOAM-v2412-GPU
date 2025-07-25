/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

Description
    The search for nearest point on an ellipse or ellipsoid follows the
    description given by Geometric Tools (David Eberly), which also
    include some pseudo code. The content is CC-BY 4.0

https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf

    In the search algorithm, symmetry is exploited and the searching is
    confined to the first (+x,+y,+z) octant, and the radii are ordered
    from largest to smallest.

\*---------------------------------------------------------------------------*/

#include "searchableSphere.H"
#include "addToRunTimeSelectionTable.H"
#include <array>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSphere, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSphere,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSphere,
        dict,
        sphere
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// General handling
namespace Foam
{

// Dictionary entry with single scalar or vector quantity
inline static vector getRadius(const word& name, const dictionary& dict)
{
    if (token(dict.lookup(name)).isNumber())
    {
        return vector::uniform(dict.get<scalar>(name));
    }

    return dict.get<vector>(name);
}


// Test point for negative components, return the sign-changes
inline static unsigned getOctant(const point& p)
{
    unsigned octant = 0;

    if (p.x() < 0) { octant |= 1; }
    if (p.y() < 0) { octant |= 2; }
    if (p.z() < 0) { octant |= 4; }

    return octant;
}


// Apply sign-changes to point
inline static void applyOctant(point& p, unsigned octant)
{
    if (octant & 1) { p.x() = -p.x(); }
    if (octant & 2) { p.y() = -p.y(); }
    if (octant & 4) { p.z() = -p.z(); }
}


// Vector magnitudes

inline static scalar vectorMagSqr
(
    const scalar x,
    const scalar y
)
{
    return (sqr(x) + sqr(y));
}

inline static scalar vectorMagSqr
(
    const scalar x,
    const scalar y,
    const scalar z
)
{
    return (sqr(x) + sqr(y) + sqr(z));
}

inline static scalar vectorMag
(
    const scalar x,
    const scalar y
)
{
    return hypot(x, y);
}

inline static scalar vectorMag
(
    const scalar x,
    const scalar y,
    const scalar z
)
{
    return ::sqrt(vectorMagSqr(x, y, z));
}


} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Searching
namespace Foam
{

// Max iterations for root finding
static constexpr int maxIters = 100;

// Relative ellipse size within the root finding (1)
static constexpr scalar tolCloseness = 1e-3;


// Find root for distance to ellipse
static scalar findRootEllipseDistance
(
    const scalar r0,    //!< Ratio of major/minor
    const scalar z0,    //!< Search point y0, scaled by e0 (major)
    const scalar z1,    //!< Search point y1, scaled by e1 (minor)
    scalar g            //!< Evaluated ellipse, implicit form
)
{
    const scalar n0 = r0*z0;

    scalar s0 = z1 - 1;
    scalar s1 = (g < 0 ? 0 : vectorMag(n0, z1) - 1);
    scalar s = 0;

    int nIters = 0;
    while (nIters++ < maxIters)
    {
        s = (s0 + s1) / 2;
        if (equal(s, s0) || equal(s, s1))
        {
            break;
        }

        g = sqr(n0/(s+r0)) + sqr(z1/(s+1)) - 1;

        if (mag(g) < tolCloseness)
        {
            break;
        }
        else if (g > 0)
        {
            s0 = s;
        }
        else  // g < 0
        {
            s1 = s;
        }
    }

    #ifdef FULLDEBUG
    InfoInFunction
        << "Located root in " << nIters << " iterations" << endl;
    #endif

    return s;
}


// Find root for distance to ellipsoid
static scalar findRootEllipsoidDistance
(
    const scalar r0,    //!< Ratio of major/minor
    const scalar r1,    //!< Ratio of mezzo/minor
    const scalar z0,    //!< Search point y0, scaled by e0 (major)
    const scalar z1,    //!< Search point y1, scaled by e1 (mezzo)
    const scalar z2,    //!< Search point y2, scaled by e2 (minor)
    scalar g            //!< Evaluated ellipsoid, implicit form
)
{
    const scalar n0 = r0*z0;
    const scalar n1 = r1*z1;

    scalar s0 = z2 - 1;
    scalar s1 = (g < 0 ? 0 : vectorMag(n0, n1, z2) - 1);
    scalar s = 0;

    int nIters = 0;
    while (nIters++ < maxIters)
    {
        s = (s0 + s1) / 2;
        if (equal(s, s0) || equal(s, s1))
        {
            break;
        }

        g = vectorMagSqr(n0/(s+r0), n1/(s+r1), z2/(s+1)) - 1;

        if (mag(g) < tolCloseness)
        {
            break;
        }
        else if (g > 0)
        {
            s0 = s;
        }
        else  // g < 0
        {
            s1 = s;
        }
    }

    #ifdef FULLDEBUG
    InfoInFunction
        << "root at " << s << " found in " << nIters
        << " iterations" << endl;
    #endif

    return s;
}


// Distance (squared) to an ellipse (2D)
static scalar distanceToEllipse
(
    // [in] Ellipse characteristics. e0 >= e1
    const scalar e0, const scalar e1,
    // [in] search point. y0 >= 0, y1 >= 0
    const scalar y0, const scalar y1,
    // [out] nearest point on ellipse
    scalar& x0, scalar& x1
)
{
    if (equal(y1, 0))
    {
        // On the y1 = 0 axis

        const scalar numer0 = e0*y0;
        const scalar denom0 = sqr(e0) - sqr(e1);

        if (numer0 < denom0)
        {
            const scalar xde0 = numer0/denom0;  // Is always < 1

            x0 = e0*xde0;
            x1 = e1*sqrt(1 - sqr(xde0));

            return vectorMagSqr((x0-y0), x1);
        }

        // Fallthrough
        x0 = e0;
        x1 = 0;

        return sqr(y0-e0);
    }
    else if (equal(y0, 0))
    {
        // On the y0 = 0 axis, in the y1 > 0 half
        x0 = 0;
        x1 = e1;
        return sqr(y1-e1);
    }
    else
    {
        // In the y0, y1 > 0 quadrant

        const scalar z0 = y0 / e0;
        const scalar z1 = y1 / e1;
        scalar eval = sqr(z0) + sqr(z1);

        scalar g = eval - 1;

        if (mag(g) < tolCloseness)
        {
            x0 = y0;
            x1 = y1;

            if (!equal(eval, 1))
            {
                // Very close, scale accordingly.
                eval = sqrt(eval);
                x0 /= eval;
                x1 /= eval;
            }

            return sqr(x0-y0) + sqr(x1-y1);
        }


        // General search.
        // Uses root find to get tbar of F(t) on (-e1^2,+inf)

        // Ratio major/minor
        const scalar r0 = sqr(e0 / e1);

        const scalar sbar =
            findRootEllipseDistance(r0, z0, z1, g);

        x0 = r0 * y0 / (sbar + r0);
        x1 = y1 / (sbar + 1);

        // Re-evaluate
        eval = sqr(x0/e0) + sqr(x1/e1);

        if (!equal(eval, 1))
        {
            // Very close, scale accordingly.
            //
            // This is not exact - the point is projected at a
            // slight angle, but we are only correcting for
            // rounding in the first place.

            eval = sqrt(eval);
            x0 /= eval;
            x1 /= eval;
        }

        return sqr(x0-y0) + sqr(x1-y1);
    }

    // Code never reaches here
    FatalErrorInFunction
        << "Programming/logic error" << nl
        << exit(FatalError);

    return 0;
}


// Distance (squared) to an ellipsoid (3D)
static scalar distanceToEllipsoid
(
    // [in] Ellipsoid characteristics. e0 >= e1 >= e2
    const scalar e0, const scalar e1, const scalar e2,
    // [in] search point. y0 >= 0, y1 >= 0, y2 >= 0
    const scalar y0, const scalar y1, const scalar y2,
    // [out] nearest point on ellipsoid
    scalar& x0, scalar& x1, scalar& x2
)
{
    if (equal(y2, 0))
    {
        // On the y2 = 0 plane. Can use 2D ellipse finding

        const scalar numer0 = e0*y0;
        const scalar numer1 = e1*y1;
        const scalar denom0 = sqr(e0) - sqr(e2);
        const scalar denom1 = sqr(e1) - sqr(e2);

        if (numer0 < denom0 && numer1 < denom1)
        {
            const scalar xde0 = numer0/denom0;  // Is always < 1
            const scalar xde1 = numer1/denom1;  // Is always < 1

            const scalar disc = (1 - sqr(xde0) - sqr(xde1));

            if (disc > 0)
            {
                x0 = e0*xde0;
                x1 = e1*xde1;
                x2 = e2*sqrt(disc);

                return vectorMagSqr((x0-y0), (x1-y1), x2);
            }
        }

        // Fallthrough - use 2D form
        x2 = 0;
        return distanceToEllipse(e0,e1, y0,y1, x0,x1);
    }
    else if (equal(y1, 0))
    {
        // On the y1 = 0 plane, in the y2 > 0 half
        x1 = 0;
        if (equal(y0, 0))
        {
            x0 = 0;
            x2 = e2;
            return sqr(y2-e2);
        }
        else  // y0 > 0
        {
            return distanceToEllipse(e0,e2, y0,y2, x0,x2);
        }
    }
    else if (equal(y0, 0))
    {
        // On the y1 = 0 plane, in the y1, y2 > 0 quadrant
        x0 = 0;
        return distanceToEllipse(e1,e2, y1,y2, x1,x2);
    }
    else
    {
        // In the y0, y1, y2 > 0 octant

        const scalar z0 = y0/e0;
        const scalar z1 = y1/e1;
        const scalar z2 = y2/e2;
        scalar eval = vectorMagSqr(z0, z1, z2);

        scalar g = eval - 1;

        if (mag(g) < tolCloseness)
        {
            x0 = y0;
            x1 = y1;
            x2 = y2;

            if (equal(eval, 1))
            {
                // Exactly on the ellipsoid - we are done
                return 0;
            }

            // Very close, scale accordingly.
            eval = sqrt(eval);
            x0 /= eval;
            x1 /= eval;
            x2 /= eval;

            return vectorMagSqr((x0-y0), (x1-y1), (x2-y2));
        }


        // General search.
        // Compute the unique root tbar of F(t) on (-e2^2,+inf)

        const scalar r0 = sqr(e0/e2);
        const scalar r1 = sqr(e1/e2);

        const scalar sbar =
            findRootEllipsoidDistance(r0,r1, z0,z1,z2, g);

        x0 = r0*y0/(sbar+r0);
        x1 = r1*y1/(sbar+r1);
        x2 = y2/(sbar+1);

        // Reevaluate
        eval = vectorMagSqr((x0/e0), (x1/e1), (x2/e2));

        if (!equal(eval, 1))
        {
            // Not exactly on ellipsoid?
            //
            // Scale accordingly. This is not exact - the point
            // is projected at a slight angle, but we are only
            // correcting for rounding in the first place.

            eval = sqrt(eval);
            x0 /= eval;
            x1 /= eval;
            x2 /= eval;
        }

        return vectorMagSqr((x0-y0), (x1-y1), (x2-y2));
    }

    // Code never reaches here
    FatalErrorInFunction
        << "Programming/logic error" << nl
        << exit(FatalError);

    return 0;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::searchableSphere::componentOrder
Foam::searchableSphere::getOrdering(const vector& radii)
{
    #ifdef FULLDEBUG
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (radii[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Radii must be positive, non-zero: " << radii << endl
                << exit(FatalError);
        }
    }
    #endif

    std::array<uint8_t, 3> idx{0, 1, 2};

    // Reverse sort by magnitude (largest first...)
    // Radii are positive (checked above, or just always true)
    std::stable_sort
    (
        idx.begin(),
        idx.end(),
        [&](uint8_t a, uint8_t b){ return radii[a] > radii[b]; }
    );

    componentOrder order{idx[0], idx[1], idx[2], shapeType::GENERAL};

    if (equal(radii[order.major], radii[order.minor]))
    {
        order.shape = shapeType::SPHERE;
    }
    else if (equal(radii[order.major], radii[order.mezzo]))
    {
        order.shape = shapeType::OBLATE;
    }
    else if (equal(radii[order.mezzo], radii[order.minor]))
    {
        order.shape = shapeType::PROLATE;
    }

    return order;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableSphere::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    // Handle special cases first

    if (order_.shape == shapeType::SPHERE)
    {
        // Point relative to origin, simultaneously the normal on the sphere
        const vector n(sample - origin_);
        const scalar magN = mag(n);

        // It is a sphere, all radii are identical

        if (nearestDistSqr >= sqr(magN - radii_[0]))
        {
            info.hitPoint
            (
                origin_
              + (magN < ROOTVSMALL ? vector(radii_[0],0,0) : (radii_[0]*n/magN))
            );
            info.setIndex(0);
        }

        return info;
    }


    //
    // Non-sphere
    //

    // Local point relative to the origin
    vector relPt(sample - origin_);

    // Detect -ve octants
    const unsigned octant = getOctant(relPt);

    // Flip everything into positive octant.
    // That is what the algorithm expects.
    applyOctant(relPt, octant);


    // TODO - quick reject for things that are too far away

    point& near = info.point();
    scalar distSqr{0};

    if (order_.shape == shapeType::OBLATE)
    {
        // Oblate (major = mezzo > minor) - use 2D algorithm
        // Distance from the minor axis to relPt
        const scalar axisDist = hypot(relPt[order_.major], relPt[order_.mezzo]);

        // Distance from the minor axis to near
        scalar nearAxis;

        distSqr = distanceToEllipse
        (
            radii_[order_.major], radii_[order_.minor],
            axisDist,             relPt[order_.minor],
            nearAxis,             near[order_.minor]
        );

        // Now nearAxis is the ratio, by which their components have changed
        nearAxis /= (axisDist + VSMALL);

        near[order_.major] = relPt[order_.major] * nearAxis;
        near[order_.mezzo] = relPt[order_.mezzo] * nearAxis;
        // near[order_.minor] = already calculated
    }
    else if (order_.shape == shapeType::PROLATE)
    {
        // Prolate (major > mezzo = minor) - use 2D algorithm
        // Distance from the major axis to relPt
        const scalar axisDist = hypot(relPt[order_.mezzo], relPt[order_.minor]);

        // Distance from the major axis to near
        scalar nearAxis;

        distSqr = distanceToEllipse
        (
            radii_[order_.major], radii_[order_.minor],
            relPt[order_.major],  axisDist,
            near[order_.major],   nearAxis
        );

        // Now nearAxis is the ratio, by which their components have changed
        nearAxis /= (axisDist + VSMALL);

        // near[order_.major] = already calculated
        near[order_.mezzo] = relPt[order_.mezzo] * nearAxis;
        near[order_.minor] = relPt[order_.minor] * nearAxis;
    }
    else  // General case
    {
        distSqr = distanceToEllipsoid
        (
            radii_[order_.major], radii_[order_.mezzo], radii_[order_.minor],
            relPt[order_.major],  relPt[order_.mezzo],  relPt[order_.minor],
            near[order_.major],   near[order_.mezzo],   near[order_.minor]
        );
    }

    // Flip everything back to original octant
    applyOctant(near, octant);

    // From local to global
    near += origin_;


    // Accept/reject based on distance
    if (distSqr <= nearestDistSqr)
    {
        info.setHit();
    }

    return info;
}


// From Graphics Gems - intersection of sphere with ray
void Foam::searchableSphere::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    if (order_.shape == shapeType::SPHERE)
    {
        vector dir(end-start);
        const scalar magSqrDir = magSqr(dir);

        if (magSqrDir > ROOTVSMALL)
        {
            dir /= Foam::sqrt(magSqrDir);

            const vector relStart(start - origin_);

            const scalar v = -(relStart & dir);

            const scalar disc = sqr(radius()) - (magSqr(relStart) - sqr(v));

            if (disc >= 0)
            {
                const scalar d = Foam::sqrt(disc);

                const scalar nearParam = v - d;
                const scalar farParam = v + d;

                if (nearParam >= 0 && sqr(nearParam) <= magSqrDir)
                {
                    near.hitPoint(start + nearParam*dir, 0);
                }

                if (farParam >= 0 && sqr(farParam) <= magSqrDir)
                {
                    far.hitPoint(start + farParam*dir, 0);
                }
            }
        }

        return;
    }


    // General case

    // Similar to intersection of sphere with ray (Graphics Gems),
    // but we scale x/y/z components according to radii
    // to have a unit spheroid for the interactions.
    // When finished, we unscale to get the real points

    // Note - can also be used for the spherical case

    const point relStart = scalePoint(start);

    vector dir(scalePoint(end) - relStart);
    const scalar magSqrDir = magSqr(dir);

    if (magSqrDir > ROOTVSMALL)
    {
        dir /= Foam::sqrt(magSqrDir);

        const scalar v = -(relStart & dir);

        const scalar disc = scalar(1) - (magSqr(relStart) - sqr(v));

        if (disc >= 0)
        {
            const scalar d = Foam::sqrt(disc);

            const scalar nearParam = v - d;
            const scalar farParam = v + d;

            if (nearParam >= 0 && sqr(nearParam) <= magSqrDir)
            {
                near.hitPoint(unscalePoint(relStart + nearParam*dir), 0);
            }
            if (farParam >= 0 && sqr(farParam) <= magSqrDir)
            {
                far.hitPoint(unscalePoint(relStart + farParam*dir), 0);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const point& origin,
    const scalar radius
)
:
    searchableSphere(io, origin, vector::uniform(radius))
{}


Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const point& origin,
    const vector& radii
)
:
    searchableSurface(io),
    origin_(origin),
    radii_(radii),
    order_(getOrdering(radii_))  // NB: use () not {} for copy initialization
{
    bounds().min() = (centre() - radii_);
    bounds().max() = (centre() + radii_);
}


Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSphere
    (
        io,
        dict.getCompat<vector>("origin", {{"centre", -1806}}),
        getRadius("radius", dict)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::searchableSphere::surfacePoint
(
    const scalar theta,
    const scalar phi
) const
{
    return point
    (
        origin_.x() + radii_.x() * cos(theta)*sin(phi),
        origin_.y() + radii_.y() * sin(theta)*sin(phi),
        origin_.z() + radii_.z() * cos(phi)
    );
}


Foam::vector Foam::searchableSphere::surfaceNormal
(
    const scalar theta,
    const scalar phi
) const
{
    // Normal is (x0/r0^2, x1/r1^2, x2/r2^2)

    return vector
    (
        cos(theta)*sin(phi) / radii_.x(),
        sin(theta)*sin(phi) / radii_.y(),
        cos(phi) / radii_.z()
    ).normalise();
}


bool Foam::searchableSphere::overlaps(const boundBox& bb) const
{
    if (order_.shape == shapeType::SPHERE)
    {
        return bb.overlaps(origin_, sqr(radius()));
    }

    if (!bb.good())
    {
        return false;
    }

    // Code largely as per
    // boundBox::overlaps(const point& centre, const scalar radiusSqr)
    // but normalized for a unit size

    // Find out where centre is in relation to bb.
    // Find nearest point on bb.

    // Note: no major advantage in treating sphere specially

    scalar distSqr = 0;
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar d0 = bb.min()[dir] - origin_[dir];
        const scalar d1 = bb.max()[dir] - origin_[dir];

        if ((d0 > 0) == (d1 > 0))
        {
            // Both min/max are on the same side of the origin
            // ie, box does not span spheroid in this direction

            if (Foam::mag(d0) < Foam::mag(d1))
            {
                distSqr += Foam::sqr(d0/radii_[dir]);
            }
            else
            {
                distSqr += Foam::sqr(d1/radii_[dir]);
            }

            if (distSqr > 1)
            {
                return false;
            }
        }
    }

    return true;
}


const Foam::wordList& Foam::searchableSphere::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_.first() = "region0";
    }
    return regions_;
}


void Foam::searchableSphere::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize(1);
    radiusSqr.resize(1);

    centres[0] = origin_;
    radiusSqr[0] = Foam::sqr(radius());

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


void Foam::searchableSphere::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.resize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableSphere::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection.
        // If none intersected take second one.

        findLineAll(start[i], end[i], info[i], b);

        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSphere::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection.
        // Discard far intersection

        findLineAll(start[i], end[i], info[i], b);

        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSphere::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.resize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;

        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].resize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].resize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].resize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableSphere::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize(info.size());
    region = 0;
}


void Foam::searchableSphere::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.resize(info.size());

    forAll(info, i)
    {
        if (info[i].hit())
        {
            if (order_.shape == shapeType::SPHERE)
            {
                // Special case (sphere)
                normal[i] = normalised(info[i].point() - origin_);
            }
            else
            {
                // General case
                // Normal is (x0/r0^2, x1/r1^2, x2/r2^2)

                normal[i] = scalePoint(info[i].point());

                normal[i].x() /= radii_.x();
                normal[i].y() /= radii_.y();
                normal[i].z() /= radii_.z();
                normal[i].normalise();
            }
        }
        else
        {
            // Set to what?
            normal[i] = Zero;
        }
    }
}


void Foam::searchableSphere::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.resize(points.size());

    if (order_.shape == shapeType::SPHERE)
    {
        // Special case. Minor advantage in treating specially

        const scalar rad2 = sqr(radius());

        forAll(points, pointi)
        {
            const point& p = points[pointi];

            volType[pointi] =
            (
                (magSqr(p - origin_) <= rad2)
              ? volumeType::INSIDE : volumeType::OUTSIDE
            );
        }

        return;
    }

    // General case - could also do component-wise (manually)
    // Evaluate:  (x/r0)^2 + (y/r1)^2 + (z/r2)^2 - 1 = 0
    // [sphere]:  x^2 + y^2 + z^2 - R^2 = 0

    forAll(points, pointi)
    {
        const point p = scalePoint(points[pointi]);

        volType[pointi] =
        (
            (magSqr(p) <= 1)
          ? volumeType::INSIDE : volumeType::OUTSIDE
        );
    }
}


// ************************************************************************* //
