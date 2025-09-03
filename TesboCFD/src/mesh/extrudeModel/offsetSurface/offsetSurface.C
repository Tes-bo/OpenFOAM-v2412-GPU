/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "offsetSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "barycentric2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(offsetSurface, 0);

addToRunTimeSelectionTable(extrudeModel, offsetSurface, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

offsetSurface::offsetSurface(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    project_(coeffDict_.getOrDefault("project", false))
{
    // Read surface
    fileName baseName(coeffDict_.get<fileName>("baseSurface").expand());
    baseSurfPtr_.reset(new triSurface(baseName));

    // Read offsetted surface
    fileName offsetName(coeffDict_.get<fileName>("offsetSurface").expand());
    offsetSurfPtr_.reset(new triSurface(offsetName));


    const triSurface& b = *baseSurfPtr_;
    const triSurface& o = *offsetSurfPtr_;

    if
    (
        b.size() != o.size()
     || b.nPoints() != o.nPoints()
     || b.nEdges() != o.nEdges()
    )
    {
        FatalIOErrorInFunction(dict)
            << "offsetSurface:\n    " << offsetName
            << " has different topology than the baseSurface:\n    "
            << baseName << endl
            << exit(FatalIOError);
    }

    // Construct search engine
    baseSearchPtr_.reset(new triSurfaceSearch(b));

    // Construct search engine
    offsetSearchPtr_.reset(new triSurfaceSearch(o));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

offsetSurface::~offsetSurface()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point offsetSurface::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    if (layer == 0)
    {
        return surfacePoint;
    }
    else
    {
        pointField samples(1, surfacePoint);
        scalarField nearestDistSqr(1, GREAT);
        List<pointIndexHit> info;
        baseSearchPtr_().findNearest(samples, nearestDistSqr, info);

        label triI = info[0].index();


        const triSurface& base = baseSurfPtr_();
        const triPointRef baseTri(base[triI].tri(base.points()));
        const barycentric2D bary = baseTri.pointToBarycentric(surfacePoint);

        const triSurface& offset = offsetSurfPtr_();
        const triPointRef offsetTri(offset[triI].tri(offset.points()));

        const point offsetPoint
        (
            bary[0]*offsetTri.a()
          + bary[1]*offsetTri.b()
          + bary[2]*offsetTri.c()
        );

        point interpolatedPoint
        (
            surfacePoint + sumThickness(layer)*(offsetPoint-surfacePoint)
        );


        // Either return interpolatedPoint or re-project onto surface (since
        // snapping might not have do so exactly)

        if (project_)
        {
            // Re-project onto surface
            offsetSearchPtr_().findNearest
            (
                pointField(1, interpolatedPoint),
                scalarField(1, GREAT),
                info
            );
            return info[0].hitPoint();
        }
        else
        {
            return interpolatedPoint;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
