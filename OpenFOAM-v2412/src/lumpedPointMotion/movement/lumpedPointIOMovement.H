/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::lumpedPointIOMovement

Description
    IO-registered version of lumpedPointMovement.

SourceFiles
    lumpedPointMovement.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointIOMovement_H
#define lumpedPointIOMovement_H

#include "lumpedPointMovement.H"
#include "regIOobject.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class lumpedPointIOMovement Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointIOMovement
:
    public lumpedPointMovement,
    public regIOobject
{
    // Private Member Functions

        //- No copy construct
        lumpedPointIOMovement(const lumpedPointIOMovement&) = delete;

        //- No copy assignment
        void operator=(const lumpedPointIOMovement&) = delete;


public:

    //- Runtime type information
    TypeNameNoDebug("lumpedPointMovement");


    // Static Member Functions

        //- Find the movement object or nullptr if not found
        static lumpedPointIOMovement* getMovementObject
        (
            const objectRegistry& obr
        );

        //- Create a movement object in the registry by
        //- reading system dictionary
        static autoPtr<lumpedPointIOMovement> New
        (
            const objectRegistry& obr,
            label ownerId = -1
        );


    // Constructors

        //- Construct from IOobject, optionally with some owner information
        explicit lumpedPointIOMovement(const IOobject& io, label ownerId = -1);


    //- Destructor
    virtual ~lumpedPointIOMovement() = default;


    // Member Functions

        //- The readData member function used by regIOobject
        bool readData(Istream& is);

        //- The writeData member function required by regIOobject
        bool writeData(Ostream& os) const;
};


// IOstream Operators
Ostream& operator<<(Ostream& os, const lumpedPointIOMovement& obj);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
