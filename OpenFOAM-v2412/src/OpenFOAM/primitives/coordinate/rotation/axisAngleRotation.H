/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::coordinateRotations::axisAngle

Description
    A coordinateRotation specified by a rotation axis and a rotation angle
    about that axis.

    \verbatim
    rotation
    {
        type        axisAngle;
        axis        (1 0 0);
        angle       90;
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd | Default
        type        | Type name: axisAngle                  | yes  |
        axis        | Axis of rotation (vector)             | yes  |
        angle       | Rotation angle                        | yes  |
        degrees     | The angle is in degrees               | no   | true
    \endtable

Note
    The rotation axis is normalized internally.

SourceFiles
    axisAngleRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_axisAngle_H
#define Foam_coordinateRotations_axisAngle_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
               Class coordinateRotations::axisAngle Declaration
\*---------------------------------------------------------------------------*/

class axisAngle
:
    public coordinateRotation
{
    // Private Data

        //- The rotation axis
        vector axis_;

        //- The rotation angle
        scalar angle_;

        //- Angle measured in degrees
        bool degrees_;


    // Private Member Functions

        //- Check specification for an identity rotation
        void checkSpec();


public:

    //- Runtime type information
    TypeNameNoDebug("axisAngle");


    // Constructors

        //- Default construct. Axis = Z, angle = 0.
        axisAngle();

        //- Copy construct
        axisAngle(const axisAngle& crot);

        //- Construct from axis and angle
        axisAngle(const vector& axis, scalar angle, bool degrees);

        //- Construct from x/y/z axis enumeration and angle
        axisAngle(const vector::components axis, scalar angle, bool degrees);

        //- Construct from dictionary
        explicit axisAngle(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~axisAngle() = default;


    // Static Member Functions

        //- The rotation tensor for given axis/angle
        static tensor rotation
        (
            const vector& axis,
            const scalar angle,
            bool degrees=false
        );

        //- Rotation tensor calculated for given axis and angle
        static tensor rotation
        (
            const vector::components axis,
            const scalar angle,
            bool degrees=false
        );


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- Calculate and return the rotation tensor
        //- calculated from axis and angle
        virtual tensor R() const;

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
