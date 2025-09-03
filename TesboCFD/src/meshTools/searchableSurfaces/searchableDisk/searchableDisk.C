/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "searchableDisk.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableDisk, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableDisk,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableDisk,
        dict,
        disk
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableDisk::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    vector v(sample - origin());

    // Decompose sample-origin into normal and parallel component
    const scalar parallel = (v & normal());

    // Remove the parallel component and normalise
    v -= parallel * normal();

    const scalar magV = mag(v);

    v.normalise();

    // Clip to inner/outer radius
    info.setPoint(origin() + radialLimits_.clamp(magV)*v);

    if (info.point().distSqr(sample) < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


void Foam::searchableDisk::findLine
(
    const point& start,
    const point& end,
    pointIndexHit& info
) const
{
    info = pointIndexHit(false, Zero, -1);

    vector v(start - origin());

    // Decompose sample-origin into normal and parallel component
    const scalar parallel = (v & normal());

    if (Foam::sign(parallel) == Foam::sign(plane::signedDistance(end)))
    {
        return;
    }

    // Remove the parallel component and normalise
    v -= parallel * normal();

    const scalar magV = mag(v);

    v.normalise();

    // Set (hit or miss) to intersection of ray and plane of disk
    info.setPoint(origin() + magV*v);

    if (radialLimits_.contains(magV))
    {
        info.setHit();
        info.setIndex(0);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableDisk::searchableDisk
(
    const IOobject& io,
    const point& originPoint,
    const vector& normalVector,
    const scalar outerRadius,
    const scalar innerRadius
)
:
    searchableSurface(io),
    plane(originPoint, normalVector),
    radialLimits_(innerRadius, outerRadius)
{
    // Rough approximation of bounding box

    // See searchableCylinder
    vector span
    (
        sqrt(sqr(normal().y()) + sqr(normal().z())),
        sqrt(sqr(normal().x()) + sqr(normal().z())),
        sqrt(sqr(normal().x()) + sqr(normal().y()))
    );
    span *= outerRadius;

    bounds().min() = origin() - span;
    bounds().max() = origin() + span;
}


Foam::searchableDisk::searchableDisk
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableDisk
    (
        io,
        dict.get<point>("origin"),
        dict.get<vector>("normal"),
        dict.get<scalar>("radius"),
        dict.getOrDefault<scalar>("innerRadius", 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableDisk::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_.first() = "region0";
    }
    return regions_;
}


void Foam::searchableDisk::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize(1);
    radiusSqr.resize(1);

    centres[0] = origin();
    radiusSqr[0] = sqr(radialLimits_.max());

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


void Foam::searchableDisk::findNearest
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


void Foam::searchableDisk::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.resize(start.size());

    forAll(start, i)
    {
        findLine(start[i], end[i], info[i]);
    }
}


void Foam::searchableDisk::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine(start, end, info);
}


void Foam::searchableDisk::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit inter;
        findLine(start[i], end[i], inter);

        if (inter.hit())
        {
            info[i].setSize(1);
            info[i][0] = inter;
        }
        else
        {
            info[i].clear();
        }
    }
}


void Foam::searchableDisk::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.resize(info.size());
    region = 0;
}


void Foam::searchableDisk::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normals
) const
{
    normals.resize(info.size());
    normals = normal();
}


void Foam::searchableDisk::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for disk."
        << exit(FatalError);
}


// ************************************************************************* //
