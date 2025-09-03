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
    Foam::coordinateRotations::starcd

Description
    A coordinateRotation defined by the STAR-CD convention.

    The 3 rotations are defined in the STAR-CD convention
    (around Z, around X' and around Y'').
    The order of the parameter arguments matches this rotation order.

    - the rotation angles are in degrees, unless otherwise explicitly specified:

    \verbatim
    rotation
    {
        type    starcd;
        angles  (0 0 180);
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd | Default
        type        | Type name: starcd                     | yes  |
        angles      | The z-x-y rotation angles             | yes  |
        degrees     | Angles are in degrees                 | no   | true
    \endtable

Note
    Also accepts "STARCDRotation" (OpenFOAM-v1806) for the type.

SourceFiles
    STARCDCoordinateRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_starcd_H
#define Foam_coordinateRotations_starcd_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                 Class coordinateRotations::starcd Declaration
\*---------------------------------------------------------------------------*/

class starcd
:
    public coordinateRotation
{
    // Private Data

        //- The rotation angles
        vector angles_;

        //- Angles measured in degrees
        bool degrees_;


public:

    //- Runtime type information
    TypeNameNoDebug("starcd");


    // Constructors

        //- Default construct - an identity transform
        starcd();

        //- Copy construct
        starcd(const starcd& crot);

        //- Construct from rotation vector
        starcd(const vector& rotZrotXrotY, bool degrees);

        //- Construct from components of rotation vector
        starcd(scalar rotZ, scalar rotX, scalar rotY, bool degrees);

        //- Construct from dictionary
        explicit starcd(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~starcd() = default;


    // Static Member Functions

        //- Rotation tensor calculated for the specified STARCD angles
        //- interpreted as rotate-Z, rotate-X, rotate-Y
        static tensor rotation(const vector& angles, bool degrees);


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- The rotation tensor calculated for the specified STARCD angles.
        virtual tensor R() const;

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations

//- Compatibility typedef 1806
typedef coordinateRotations::starcd STARCDCoordinateRotation;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
