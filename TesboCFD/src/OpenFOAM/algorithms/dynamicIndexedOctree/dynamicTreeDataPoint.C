/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "dynamicTreeDataPoint.H"
#include "dynamicIndexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(dynamicTreeDataPoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicTreeDataPoint::dynamicTreeDataPoint
(
    const DynamicList<point>& points
)
:
    points_(points)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox
Foam::dynamicTreeDataPoint::bounds(const labelUList& indices) const
{
    return treeBoundBox(points_, indices);
}


Foam::volumeType Foam::dynamicTreeDataPoint::getVolumeType
(
    const dynamicIndexedOctree<dynamicTreeDataPoint>& oc,
    const point& sample
) const
{
    return volumeType::UNKNOWN;
}


bool Foam::dynamicTreeDataPoint::overlaps
(
    const label index,
    const treeBoundBox& searchBox
) const
{
    return searchBox.contains(centre(index));
}


bool Foam::dynamicTreeDataPoint::overlaps
(
    const label index,
    const point& centre,
    const scalar radiusSqr
) const
{
    return (centre.distSqr(this->centre(index)) <= radiusSqr);
}


void Foam::dynamicTreeDataPoint::findNearest
(
    const labelUList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    for (const label index : indices)
    {
        const point& pt = centre(index);

        const scalar distSqr = sample.distSqr(pt);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = pt;
        }
    }
}


void Foam::dynamicTreeDataPoint::findNearest
(
    const labelUList& indices,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& minIndex,
    point& linePoint,
    point& nearestPoint
) const
{
    const treeBoundBox lnBb(ln.box());

    // Best so far
    scalar nearestDistSqr = linePoint.distSqr(nearestPoint);

    for (const label index : indices)
    {
        const point& pt = centre(index);

        if (tightest.contains(pt))
        {
            // Nearest point on line
            pointHit pHit = ln.nearestDist(pt);
            const scalar distSqr = sqr(pHit.distance());

            if (distSqr < nearestDistSqr)
            {
                nearestDistSqr = distSqr;
                minIndex = index;
                linePoint = pHit.point();
                nearestPoint = pt;

                tightest = lnBb;
                tightest.grow(pHit.distance());
            }
        }
    }
}


// ************************************************************************* //
