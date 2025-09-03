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
    Foam::coordinateRotations::identity

Description
    An identity coordinateRotation.

    \verbatim
    rotation
    {
        type    none;
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Reqd | Default
        type        | Type name: none                       | yes  |
    \endtable

SourceFiles
    identityRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_identity_H
#define Foam_coordinateRotations_identity_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                Class coordinateRotations::identity Declaration
\*---------------------------------------------------------------------------*/

class identity
:
    public coordinateRotation
{
public:

    //- Runtime type information
    TypeNameNoDebug("none");


    // Generated Methods

        //- Default construct
        identity() = default;

        //- Copy construct
        identity(const identity&) = default;

        //- Copy assignment
        identity& operator=(const identity&) = default;


    // Constructors

        //- Construct from dictionary
        explicit identity(const dictionary& unused);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~identity() = default;


    // Member Functions

        //- Reset specification (no-op)
        virtual void clear();

        //- Return an identity rotation tensor
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
