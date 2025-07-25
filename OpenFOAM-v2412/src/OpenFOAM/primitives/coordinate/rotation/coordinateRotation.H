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

Namespace
    Foam::coordinateRotations

Description
    Namespace for coordinate system rotations.

Class
    Foam::coordinateRotation

Description
    User specification of a coordinate rotation.

    \verbatim
    rotation
    {
        type    axes
        e1      (1 0 0);
        e2      (0 1 0);
    }
    \endverbatim

    Types of coordinateRotations:
      -# \link coordinateRotations::identity none \endlink
      -# \link coordinateRotations::axes axes \endlink
      -# \link coordinateRotations::axisAngle axisAngle \endlink
      -# \link coordinateRotations::cylindrical cylindrical \endlink
      -# \link coordinateRotations::euler euler \endlink
      -# \link coordinateRotations::starcd starcd \endlink

SourceFiles
    coordinateRotation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateRotation_H
#define Foam_coordinateRotation_H

#include "autoPtr.H"
#include "vector.H"
#include "tensor.H"
#include "dictionary.H"
#include "IOobjectOption.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class coordinateRotation
{
protected:

    // Protected Member Functions

        //- Determine best-guess for an orthogonal axis
        static vector findOrthogonal(const vector& axis);


public:

    //- Runtime type information
    TypeNameNoDebug("coordinateRotation");

    // Declare run-time constructor selection table from dictionary
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateRotation,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        // Uses all default constructors

        //- Return a clone
        virtual autoPtr<coordinateRotation> clone() const = 0;


    // Factory Methods

        //- Clone a coordinate rotation
        template<class Derived>
        static autoPtr<coordinateRotation> Clone(const Derived& crot)
        {
            return autoPtr<coordinateRotation>(new Derived(crot));
        }

        //- Select construct the specified coordinate rotation type
        //
        //  An empty modelType will be treated as "axes" (eg, e1/e3)
        static autoPtr<coordinateRotation> New
        (
            const word& modelType,
            const dictionary& dict
        );

        //- Select construct from dictionary (requires the "type" entry)
        static autoPtr<coordinateRotation> New(const dictionary& dict);


    //- Destructor
    virtual ~coordinateRotation() = default;


    // Member Functions

        //- Reset specification
        virtual void clear() = 0;

        //- Calculate and return the rotation tensor
        virtual tensor R() const = 0;


    // Write

        //- Write information
        virtual void write(Ostream& os) const = 0;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
