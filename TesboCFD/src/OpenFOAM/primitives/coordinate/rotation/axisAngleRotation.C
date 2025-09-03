/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "axisAngleRotation.H"
#include "dictionary.H"
#include "quaternion.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

    defineTypeName(axisAngle);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        axisAngle,
        dictionary
    );

} // End namespace coordinateRotation
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordinateRotations::axisAngle::checkSpec()
{
    if (mag(angle_) < VSMALL || mag(axis_) < SMALL)
    {
        clear(); // identity rotation
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::axisAngle::rotation
(
    const vector& axis,
    const scalar angle,
    bool degrees
)
{
    if (mag(angle) < VSMALL || mag(axis) < SMALL)
    {
        return sphericalTensor::I;  // identity rotation
    }

    return quaternion(axis, (degrees ? degToRad(angle) : angle)).R();
}


Foam::tensor Foam::coordinateRotations::axisAngle::rotation
(
    const vector::components axis,
    const scalar angle,
    bool degrees
)
{
    vector rotAxis(Zero);
    rotAxis[axis] = 1;

    return axisAngle::rotation(rotAxis, angle, degrees);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::axisAngle::axisAngle()
:
    coordinateRotation(),
    axis_ (0,0,1),  // e3 = global Z
    angle_(Zero),
    degrees_(true)
{}


Foam::coordinateRotations::axisAngle::axisAngle(const axisAngle& crot)
:
    coordinateRotation(),
    axis_(crot.axis_),
    angle_(crot.angle_),
    degrees_(crot.degrees_)
{
    checkSpec();
}


Foam::coordinateRotations::axisAngle::axisAngle
(
    const vector& axis,
    scalar angle,
    bool degrees
)
:
    coordinateRotation(),
    axis_(axis),
    angle_(angle),
    degrees_(degrees)
{
    checkSpec();
}


Foam::coordinateRotations::axisAngle::axisAngle
(
    const vector::components axis,
    scalar angle,
    bool degrees
)
:
    coordinateRotation(),
    axis_(Zero),
    angle_(angle),
    degrees_(degrees)
{
    axis_[axis] = 1;
}


Foam::coordinateRotations::axisAngle::axisAngle(const dictionary& dict)
:
    axisAngle
    (
        dict.get<vector>("axis"),
        dict.get<scalar>("angle"),
        dict.getOrDefault("degrees", true)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::axisAngle::clear()
{
    axis_  = vector(0,0,1);  // e3 = global Z
    angle_ = Zero;
}


Foam::tensor Foam::coordinateRotations::axisAngle::R() const
{
    return rotation(axis_, angle_, degrees_);
}


void Foam::coordinateRotations::axisAngle::write(Ostream& os) const
{
    os  << "rotation axis: " << axis_
        << " angle(" << (degrees_ ? "deg" : "rad") << "): " << angle_;
}


void Foam::coordinateRotations::axisAngle::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());
    os.writeEntry("axis",  axis_);
    os.writeEntry("angle", angle_);
    if (!degrees_)
    {
        os.writeEntry("degrees", "false");
    }

    os.endBlock();
}


// ************************************************************************* //
