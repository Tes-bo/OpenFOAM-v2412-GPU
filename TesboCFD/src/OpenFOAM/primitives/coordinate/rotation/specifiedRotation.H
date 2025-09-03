/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::coordinateRotations::specified

Description
    An user-specified coordinateRotation, primarily to be used internally
    within coding when the rotation matrix is already known but
    needs to be wrapped as a coordinateRotation for use in a coordinate
    system.

Note
    Currently no runtime selection mechanism, since specifying a
    rotation matrix via a dictionary can be fragile due to rounding
    factors and uncertainty if a forward or reverse rotation matrix
    is intended.

SourceFiles
    specifiedRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotations_specified_H
#define Foam_coordinateRotations_specified_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
               Class coordinateRotations::specified Declaration
\*---------------------------------------------------------------------------*/

class specified
:
    public coordinateRotation
{
    // Private Data

        //- The rotation matrix
        tensor Rmatrix_;


public:

    //- Runtime type information
    TypeNameNoDebug("specified");


    // Constructors

        //- Default construct - an identity matrix
        specified();

        //- Construct from transformation matrix
        explicit specified(const tensor& rot);

        //- Construct from dictionary
        explicit specified(const dictionary& unused);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return coordinateRotation::Clone(*this);
        }


    //- Destructor
    virtual ~specified() = default;


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- Return the rotation tensor
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
