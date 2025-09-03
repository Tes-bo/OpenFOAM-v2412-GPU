/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::removeRegisteredObject

Group
    grpUtilitiesFunctionObjects

Description
    Removes registered objects if present in the database.

Usage
    Example of function object specification:
    \verbatim
    removeRegisteredObject1
    {
        type        removeRegisteredObject;
        libs        (utilityFunctionObjects);
        ...
        objects     (obj1 obj2);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: removeRegisteredObject | yes |
        objects      | objects to remove       | yes         |
    \endtable

See also
    Foam::functionObject

SourceFiles
    removeRegisteredObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_removeRegisteredObject_H
#define functionObjects_removeRegisteredObject_H

#include "regionFunctionObject.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class removeRegisteredObject Declaration
\*---------------------------------------------------------------------------*/

class removeRegisteredObject
:
    public regionFunctionObject
{
    // Private Data

        //- Names of objects to control
        wordList objectNames_;


    // Private member functions

        //- No copy construct
        removeRegisteredObject(const removeRegisteredObject&) = delete;

        //- No copy assignment
        void operator=(const removeRegisteredObject&) = delete;


public:

    //- Runtime type information
    TypeName("removeRegisteredObject");


    // Constructors

        //- Construct from Time and dictionary
        removeRegisteredObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~removeRegisteredObject() = default;


    // Member Functions

        //- Read the removeRegisteredObject data
        virtual bool read(const dictionary& dict);

        //- Remove the registered objects
        virtual bool execute();

        //- Do nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
