/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "EulerCoordinateRotation.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

    defineTypeName(euler);

    // Standard short name
    addNamedToRunTimeSelectionTable
    (
        coordinateRotation,
        euler,
        dictionary,
        euler
    );

    // Long name - Compat 1806
    addAliasToRunTimeSelectionTable
    (
        coordinateRotation,
        euler,
        dictionary,
        euler,
        EulerRotation,
        1806
    );

} // End namespace coordinateRotations
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::euler::rotation
(
    const eulerOrder order,
    const vector& angles,
    bool degrees
)
{
    scalar angle1(angles.x());  // Rotation #1
    scalar angle2(angles.y());  // Rotation #2
    scalar angle3(angles.z());  // Rotation #3

    if (degrees)
    {
        angle1 *= degToRad();
        angle2 *= degToRad();
        angle3 *= degToRad();
    }

    const scalar c1(cos(angle1)); const scalar s1(sin(angle1));
    const scalar c2(cos(angle2)); const scalar s2(sin(angle2));
    const scalar c3(cos(angle3)); const scalar s3(sin(angle3));

    // https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix

    switch (order)
    {
        // Proper Euler angles

        case eulerOrder::XZX:  // X1-Z2-X3 rotation
        {
            return tensor
            (
                (   c2  ), (      -c3*s2      ), (      s2*s3        ),
                ( c1*s2 ), ( c1*c2*c3 - s1*s3 ), ( -c3*s1 - c1*c2*s3 ),
                ( s1*s2 ), ( c1*s3 + c2*c3*s1 ), (  c1*c3 - c2*s1*s3 )
            );
            break;
        }

        case eulerOrder::XYX:  // X1-Y2-X3 rotation
        {
            return tensor
            (
                (    c2  ), (      s2*s3       ), (      c3*s2 ),
                (  s1*s2 ), ( c1*c3 - c2*s1*s3 ), ( -c1*s3 - c2*c3*s1 ),
                ( -c1*s2 ), ( c3*s1 + c1*c2*s3 ), (  c1*c2*c3 - s1*s3 )
            );
            break;
        }

        case eulerOrder::YXY:  // Y1-X2-Y3 rotation
        {
            return tensor
            (
                ( c1*c3 - c2*s1*s3 ), ( s1*s2 ), ( c1*s3 + c2*c3*s1 ),
                (     s2*s3        ), (   c2  ), (    -c3*s2        ),
                ( -c3*s1 -c1*c2*s3 ), ( c1*s2 ), ( c1*c2*c3 - s1*s3 )
            );
            break;
        }

        case eulerOrder::YZY:  // Y1-Z2-Y3 rotation
        {
            return tensor
            (
                ( c1*c2*c3 - s1*s3 ), ( -c1*s2 ), ( c3*s1 + c1*c2*s3 ),
                (        c3*s2     ), (   c2   ), (     s2*s3        ),
                (-c1*s3 - c2*c3*s1 ), (  s1*s2 ), ( c1*c3 - c2*s1*s3 )
            );
            break;
        }

        case eulerOrder::ZYZ:  // Z1-Y2-Z3 rotation
        {
            return tensor
            (
                ( c1*c2*c3 - s1*s3 ), ( -c3*s1 - c1*c2*s3 ), ( c1*s2 ),
                ( c1*s3 + c2*c3*s1 ), (  c1*c3 - c2*s1*s3 ), ( s1*s2 ),
                (     -c3*s2       ), (      s2*s3 ),        (   c2  )
            );
            break;
        }

        case eulerOrder::ZXZ:  // Z1-X2-Z3 rotation
        {
            return tensor
            (
                ( c1*c3 - c2*s1*s3 ), ( -c1*s3 - c2*c3*s1 ), (  s1*s2 ),
                ( c3*s1 + c1*c2*s3 ), ( c1*c2*c3 - s1*s3  ), ( -c1*s2 ),
                (     s2*s3        ), (      c3*s2        ), (    c2  )
            );
            break;
        }


        // Tait-Bryan angles

        case eulerOrder::XZY:  // X1-Z2-Y3 rotation
        {
            return tensor
            (
                (      c2*c3       ), (  -s2  ), (       c2*s3      ),
                ( s1*s3 + c1*c3*s2 ), ( c1*c2 ), ( c1*s2*s3 - c3*s1 ),
                ( c3*s1*s2 - c1*s3 ), ( c2*s1 ), ( c1*c3 + s1*s2*s3 )
            );
            break;
        }

        case eulerOrder::XYZ:  // X1-Y2-Z3 rotation
        {
            return tensor
            (
                (      c2*c3       ), (    -c2*s3        ), (    s2  ),
                ( c1*s3 + c3*s1*s2 ), ( c1*c3 - s1*s2*s3 ), ( -c2*s1 ),
                ( s1*s3 - c1*c3*s2 ), ( c3*s1 + c1*s2*s3 ), (  c1*c2 )
            );
            break;
        }

        case eulerOrder::YXZ:  // Y1-X2-Z3 rotation
        {
            return tensor
            (
                ( c1*c3 + s1*s2*s3 ), ( c3*s1*s2 - c1*s3 ), ( c2*s1 ),
                (     c2*s3        ), (        c2*c3     ), (  -s2  ),
                ( c1*s2*s3 - c3*s1 ), ( c1*c3*s2 + s1*s3 ), ( c1*c2 )
            );
            break;
        }

        case eulerOrder::YZX:  // Y1-Z2-X3 rotation
        {
            return tensor
            (
                (  c1*c2 ), ( s1*s3 - c1*c3*s2 ), ( c3*s1 + c1*s2*s3 ),
                (  s2    ), ( c2*c3            ), ( -c2*s3           ),
                ( -c2*s1 ), ( c1*s3 + c3*s1*s2 ), ( c1*c3 - s1*s2*s3 )
            );
            break;
        }

        case eulerOrder::ZYX:  // Z1-Y2-X3 rotation
        {
            return tensor
            (
                ( c1*c2 ), ( c1*s2*s3 - c3*s1 ), ( s1*s3 + c1*c3*s2 ),
                ( c2*s1 ), ( c1*c3 + s1*s2*s3 ), ( c3*s1*s2 - c1*s3 ),
                (  -s2  ), (      c2*s3       ), (        c2*c3     )
            );
            break;
        }

        case eulerOrder::ZXY:  // Z1-X2-Y3 rotation
        {
            return tensor
            (
                ( c1*c3 - s1*s2*s3 ), ( -c2*s1 ), ( c1*s3 + c3*s1*s2 ),
                ( c3*s1 + c1*s2*s3 ), (  c1*c2 ), ( s1*s3 - c1*c3*s2 ),
                (    -c2*s3 ),        (   s2   ), (     c2*c3        )
            );
            break;
        }

        default:
            FatalErrorInFunction
                << "Unknown euler rotation order "
                << int(order) << abort(FatalError);
    }

    return sphericalTensor::I;  // identity rotation
}


Foam::tensor Foam::coordinateRotations::euler::rotation
(
    const vector& angles,
    bool degrees
)
{
    return rotation(eulerOrder::ZXZ, angles, degrees);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::euler::euler()
:
    coordinateRotation(),
    angles_(Zero),
    degrees_(true),
    order_(eulerOrder::ZXZ)
{}


Foam::coordinateRotations::euler::euler(const euler& crot)
:
    coordinateRotation(),
    angles_(crot.angles_),
    degrees_(crot.degrees_),
    order_(crot.order_)
{}


Foam::coordinateRotations::euler::euler
(
    const vector& angles,
    bool degrees
)
:
    coordinateRotation(),
    angles_(angles),
    degrees_(degrees),
    order_(eulerOrder::ZXZ)
{}


Foam::coordinateRotations::euler::euler
(
    scalar angle1,
    scalar angle2,
    scalar angle3,
    bool degrees
)
:
    coordinateRotation(),
    angles_(angle1, angle2, angle3),
    degrees_(degrees),
    order_(eulerOrder::ZXZ)
{}


Foam::coordinateRotations::euler::euler(const dictionary& dict)
:
    coordinateRotation(),
    angles_(dict.get<vector>("angles")),
    degrees_(dict.getOrDefault("degrees", true)),
    order_
    (
        quaternion::eulerOrderNames.getOrDefault
        (
            "order",
            dict,
            quaternion::eulerOrder::ZXZ
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::euler::clear()
{
    angles_ = Zero;
    degrees_ = true;
}


Foam::tensor Foam::coordinateRotations::euler::R() const
{
    return euler::rotation(order_, angles_, degrees_);
}


void Foam::coordinateRotations::euler::write(Ostream& os) const
{
    os  << "euler-angles(" << (degrees_ ? "deg" : "rad") << "): " << angles_;
}


void Foam::coordinateRotations::euler::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());
    os.writeEntry("angles", angles_);
    if (!degrees_)
    {
        os.writeEntry("degrees", "false");
    }

    // writeEntryIfDifferent, but with enumerated name
    if (order_ != eulerOrder::ZXZ)
    {
        os.writeEntry("order", quaternion::eulerOrderNames[order_]);
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
