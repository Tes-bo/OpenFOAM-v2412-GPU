/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::syncObjects

Group
    grpUtilitiesFunctionObjects

Description
    Copies content of local objectRegistry to all processors.

    Currently only supports IOFields. It is used to distribute IOFields
    from different processors/regions/worlds
    when doing loose coupling. To be used in combination with 'mapped'
    boundary conditions.

Usage
    Example of function object specification:
    \verbatim
    syncObjects1
    {
        type        syncObjects;
        libs        ("libutilityFunctionObjects.so");
        ...

        // Where is data located relative to runTime. Given as a filename
        // with every '/' indicating a sub-objectRegistry w.r.t. runTime.
        // Local data is under <root>/send/processorXXX. After execution
        // (potentially remote) data will be under the corresponding
        // <root>/receive/processorYYY objectRegistry.
        //root        "level0/level1/level2";
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required     | Default value
        type         | type name: syncObjects  | yes          |
        root         | relative location of data  | no          | ""
    \endtable

See also
    Foam::functionObject

SourceFiles
    syncObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_syncObjects_H
#define functionObjects_syncObjects_H

#include "token.H"
#include "functionObject.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::syncObjects Declaration
\*---------------------------------------------------------------------------*/

class syncObjects
:
    public functionObject
{
private:

    // Private data

        //- Top-level registry
        const objectRegistry& obr_;

        //- objectRegistry location relative to top-level
        fileName root_;


    // Private Member Functions

        //- No copy construct
        syncObjects(const syncObjects&) = delete;

        //- No copy assignment
        void operator=(const syncObjects&) = delete;


protected:

        //- Do all: synchronise all IOFields and objectRegistry
        void sync();


public:

    //- Runtime type information
    TypeName("syncObjects");


    // Constructors

        //- Construct from Time and dictionary
        syncObjects
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~syncObjects() = default;


    // Member Functions

        const objectRegistry& obr() const
        {
            return obr_;
        }

        const fileName& root() const
        {
            return root_;
        }

        //- Read the syncObjects data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
