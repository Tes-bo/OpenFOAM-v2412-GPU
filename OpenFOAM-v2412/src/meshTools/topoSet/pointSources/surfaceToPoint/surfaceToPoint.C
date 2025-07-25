/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "surfaceToPoint.H"
#include "polyMesh.H"
#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "cpuTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, surfaceToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, surfaceToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, surfaceToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, surfaceToPoint, istream);
}


Foam::topoSetSource::addToUsageTable Foam::surfaceToPoint::usage_
(
    surfaceToPoint::typeName,
    "\n    Usage: surfaceToPoint <surface> <near> <inside> <outside>\n\n"
    "    <surface> name of triSurface\n"
    "    <near> scalar; include points with coordinate <= near to surface\n"
    "    <inside> boolean; whether to include points on opposite side of"
    " surface normal\n"
    "    <outside> boolean; whether to include points on this side of"
    " surface normal\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceToPoint::combine(topoSet& set, const bool add) const
{
    cpuTime timer;

    triSurface surf(surfName_, surfType_, scale_);

    if (verbose_)
    {
        Info<< "    Read surface from " << surfName_
            << " in = "<< timer.cpuTimeIncrement() << " s" << nl << endl;
    }

    // Construct search engine on surface
    triSurfaceSearch querySurf(surf);

    if (includeInside_ || includeOutside_)
    {
        boolList pointInside(querySurf.calcInside(mesh_.points()));

        forAll(pointInside, pointi)
        {
            if (pointInside[pointi] ? includeInside_ : includeOutside_)
            {
                addOrDelete(set, pointi, add);
            }
        }
    }

    if (nearDist_ > 0)
    {
        const pointField& meshPoints = mesh_.points();

        // Box dimensions to search in octree.
        const vector span(nearDist_, nearDist_, nearDist_);

        forAll(meshPoints, pointi)
        {
            const point& pt = meshPoints[pointi];

            pointIndexHit inter = querySurf.nearest(pt, span);

            if (inter.hit() && inter.point().dist(pt) < nearDist_)
            {
                addOrDelete(set, pointi, add);
            }
        }
    }
}


void Foam::surfaceToPoint::checkSettings() const
{
    if (nearDist_ < 0 && !includeInside_ && !includeOutside_)
    {
        FatalErrorInFunction
            << "Illegal point selection specification."
            << " Result would be either all or no points." << endl
            << "Please set one of includeInside or includeOutside"
            << " to true, set nearDistance to a value > 0"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceToPoint::surfaceToPoint
(
    const polyMesh& mesh,
    const fileName& surfName,
    const scalar nearDist,
    const bool includeInside,
    const bool includeOutside
)
:
    topoSetPointSource(mesh),
    surfName_(surfName),
    surfType_(),
    scale_(-1),
    nearDist_(nearDist),
    includeInside_(includeInside),
    includeOutside_(includeOutside)
{
    checkSettings();
}


Foam::surfaceToPoint::surfaceToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    surfName_(dict.get<fileName>("file").expand()),
    surfType_(dict.getOrDefault<word>("fileType", word::null)),
    scale_(dict.getOrDefault<scalar>("scale", -1)),
    nearDist_(dict.get<scalar>("nearDistance")),
    includeInside_(dict.get<bool>("includeInside")),
    includeOutside_(dict.get<bool>("includeOutside"))
{
    checkSettings();
}


Foam::surfaceToPoint::surfaceToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    surfName_(checkIs(is)),
    surfType_(),
    scale_(-1),
    nearDist_(readScalar(checkIs(is))),
    includeInside_(readBool(checkIs(is))),
    includeOutside_(readBool(checkIs(is)))
{
    checkSettings();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding points in relation to surface "
                << surfName_ << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing points in relation to surface "
                << surfName_ << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
