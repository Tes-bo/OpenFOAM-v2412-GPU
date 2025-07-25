/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::coordinateRotations::cylindrical

Description
    A special purpose coordinateRotation that is generally for use in
    combination with a cylindricalCS when only the orientation of the
    local Z-axis is relevant.

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd | Default
        type        | Type name: cylindrical                | yes  |
        axis        | The z-axis                            | yes  |
        e3          | Alias for 'axis'                      | no   |
    \endtable

SourceFiles
    cylindricalRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_cylindrical_H
#define Foam_coordinateRotations_cylindrical_H

#include "axesRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
              Class coordinateRotations::cylindrical Declaration
\*---------------------------------------------------------------------------*/

class cylindrical
:
    public coordinateRotations::axes
{
public:

    //- Runtime type information
    TypeNameNoDebug("cylindrical");


    // Constructors

        //- Default construct - an identity transform
        cylindrical() = default;

        //- Copy construct
        cylindrical(const cylindrical& crot);

        //- Construct from single axis (as e3) using best-guess for the
        //- second axis.
        explicit cylindrical(const vector& axis);

        //- Construct from dictionary
        explicit cylindrical(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~cylindrical() = default;


    // Static Member Functions

        //- The rotation tensor calculated from axes and order.
        //  The input axes will be normalised.
        static tensor rotation(const vector& axis);


    // Member Functions

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
