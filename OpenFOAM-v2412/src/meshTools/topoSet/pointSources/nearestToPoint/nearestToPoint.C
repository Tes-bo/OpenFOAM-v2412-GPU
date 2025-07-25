/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "nearestToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearestToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, nearestToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, nearestToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, nearestToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, nearestToPoint, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        nearestToPoint,
        word,
        nearest
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        nearestToPoint,
        istream,
        nearest
    );
}


Foam::topoSetSource::addToUsageTable Foam::nearestToPoint::usage_
(
    nearestToPoint::typeName,
    "\n    Usage: nearestToPoint (pt0 .. ptn)\n\n"
    "    Select the nearest point for each of the points pt0 ..ptn\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearestToPoint::combine(topoSet& set, const bool add) const
{
    // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(points_.size());

    // Do linear search since usually just a few points.
    const pointField& pts = mesh_.points();

    forAll(points_, pointi)
    {
        if (pts.size())
        {
            label minPointi = 0;
            scalar minDistSqr = magSqr(pts[minPointi] - points_[pointi]);

            for (label i = 1; i < pts.size(); i++)
            {
                scalar distSqr = magSqr(pts[i] - points_[pointi]);

                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    minPointi = i;
                }
            }

            const point& minPt = pts[minPointi];
            nearest[pointi].first() = pointIndexHit(true, minPt, minPointi);
            nearest[pointi].second() = Tuple2<scalar, label>
            (
                magSqr(minPt-points_[pointi]),
                Pstream::myProcNo()
            );
        }
    }


    Pstream::listCombineReduce(nearest, mappedPatchBase::nearestEqOp());

    for (const auto& near : nearest)
    {
        if (near.second().second() == Pstream::myProcNo())
        {
            addOrDelete(set, near.first().index(), add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearestToPoint::nearestToPoint
(
    const polyMesh& mesh,
    const pointField& points
)
:
    topoSetPointSource(mesh),
    points_(points)
{}


Foam::nearestToPoint::nearestToPoint
(
    const polyMesh& mesh,
    pointField&& points
)
:
    topoSetPointSource(mesh),
    points_(std::move(points))
{}


Foam::nearestToPoint::nearestToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    points_(dict.get<pointField>("points"))
{}


Foam::nearestToPoint::nearestToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    points_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearestToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding points nearest to " << points_ << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing points nearest to " << points_ << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
