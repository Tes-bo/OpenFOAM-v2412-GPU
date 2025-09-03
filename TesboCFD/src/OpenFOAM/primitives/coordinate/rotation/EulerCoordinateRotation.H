/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

Class
    Foam::coordinateRotations::euler

Description
    A coordinateRotation defined in the z-x-z (intrinsic) Euler convention.

    The 3 rotations are defined in the Euler intrinsic convention
    (around Z, around X' and around Z'').
    The order of the parameter arguments matches this rotation order.

    For reference and illustration, see
    https://en.wikipedia.org/wiki/Euler_angles

    \verbatim
    rotation
    {
        type    euler;
        angles  (0 0 180);
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd | Default
        type        | Type name: euler                      | yes  |
        angles      | Rotation angles (usually z-x-z order) | yes  |
        degrees     | Angles are in degrees                 | no   | true
        order       | Rotation order                        | no   | zxz
    \endtable

Note
    The rotation order is usually z-x-z, but can also be something like
    "rollPitchYaw" etc.
    Also accepts "EulerRotation" (OpenFOAM-v1806) for the type.

SourceFiles
    EulerCoordinateRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_euler_H
#define Foam_coordinateRotations_euler_H

#include "coordinateRotation.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                 Class coordinateRotations::euler Declaration
\*---------------------------------------------------------------------------*/

class euler
:
    public coordinateRotation
{
public:

    // Public Types

    //- Euler-angle rotation order
    using eulerOrder = quaternion::eulerOrder;


private:

    // Private Data

        //- The rotation angles
        vector angles_;

        //- Angles measured in degrees
        bool degrees_;

        //- The Euler-angle rotation order (default: zxz)
        eulerOrder order_;


public:

    //- Runtime type information
    TypeNameNoDebug("euler");


    // Constructors

        //- Default construct - an identity transform
        euler();

        //- Copy construct
        euler(const euler& crot);

        //- Construct from Euler intrinsic rotation angles (z-x-z)
        euler(const vector& angles, bool degrees);

        //- Construct from Euler intrinsic rotation angles (z-x-z)
        euler(scalar angle1, scalar angle2, scalar angle3, bool degrees);

        //- Construct from dictionary
        explicit euler(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~euler() = default;


    // Static Member Functions

        //- Rotation tensor calculated for the intrinsic Euler
        //- angles in z-x-z order
        static tensor rotation(const vector& angles, bool degrees=false);

        //- Rotation tensor calculated for given order and angles
        static tensor rotation
        (
            const eulerOrder order,
            const vector& angles,
            bool degrees=false
        );


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- The rotation tensor calculated for the specified Euler angles.
        virtual tensor R() const;

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations

//- Compatibility typedef 1806
typedef coordinateRotations::euler EulerCoordinateRotation;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
