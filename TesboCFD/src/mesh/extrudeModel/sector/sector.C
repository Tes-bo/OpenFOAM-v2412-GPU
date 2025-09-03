/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "sector.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sector, 0);

addToRunTimeSelectionTable(extrudeModel, sector, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sector::sector(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    refPoint_(coeffDict_.getCompat<point>("point", {{"axisPt", -1812}})),
    axis_(coeffDict_.get<vector>("axis").normalise()),
    angle_(degToRad(coeffDict_.get<scalar>("angle")))
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point sector::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    scalar sliceAngle;

    // For the case of a single layer extrusion assume a
    // symmetric sector about the reference plane is required
    if (nLayers_ == 1)
    {
        if (layer == 0)
        {
            sliceAngle = -angle_/2.0;
        }
        else
        {
            sliceAngle = angle_/2.0;
        }
    }
    else
    {
        sliceAngle = angle_*sumThickness(layer);
    }

    // Find projection onto axis (or rather decompose surfacePoint
    // into vector along edge (proj), vector normal to edge in plane
    // of surface point and surface normal.
    point d = surfacePoint - refPoint_;

    d -= (axis_ & d)*axis_;

    scalar dMag = mag(d);

    point edgePt = surfacePoint - d;

    // Rotate point around sliceAngle.
    point rotatedPoint = edgePt;

    if (dMag > VSMALL)
    {
        vector n = (d/dMag) ^ axis_;

        rotatedPoint +=
          + cos(sliceAngle)*d
          - sin(sliceAngle)*mag(d)*n; // Use either n or surfaceNormal
    }

    return rotatedPoint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
