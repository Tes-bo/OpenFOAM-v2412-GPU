/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fieldDictionary

Description
    Read field as dictionary (without mesh).

\*---------------------------------------------------------------------------*/

#ifndef Foam_fieldDictionary_H
#define Foam_fieldDictionary_H

#include "dictionary.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldDictionary Declaration
\*---------------------------------------------------------------------------*/

class fieldDictionary
:
    public regIOobject,
    public dictionary
{
    // Private Data

        //- The specified field type
        const word fieldType_;


public:

    // Constructors

        //- Construct from IOobject and specified (overloaded) type name
        fieldDictionary(const IOobject& io, const word& fieldType)
        :
            regIOobject(io),
            dictionary(readStream(fieldType)),
            fieldType_(fieldType)
        {
            close();
        }


    // Member Functions

        //- Redefine type name to be of the instantiated type
        virtual const word& type() const
        {
            return fieldType_;
        }

        //- The writeData function required by regIOobject write operation
        virtual bool writeData(Ostream& os) const
        {
            dictionary::write(os, false);
            return os.good();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
