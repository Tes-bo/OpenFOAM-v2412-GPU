/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::writeObjects

Group
    grpUtilitiesFunctionObjects

Description
    Allows specification of different writing frequency of objects registered
    to the database.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    writeObjects1
    {
        // Mandatory entries
        type          writeObjects;
        libs          (utilityFunctionObjects);

        // Optional entries
        writeOption   <word>;

        // Conditional entries

            // Option-1
            field     <word>;

            // Option-2
            fields    (<wordRes>);

            // Option-3
            objects   (<wordRes>);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type | Reqd    | Deflt
      type       | Type name: writeObjects   | word | yes     | -
      libs       | Library name: utilityFunctionObjects | word | yes  | -
      writeOption | Select objects with the specified write mode | no | anyWrite
      field      | Name of field to write    | word | no      | -
      fields     | Names of fields to write  | wordRes | no   | -
      objects    | Names of objects to write | wordRes | no   | -
    \endtable

    Options for the \c writeOption entry:
    \vartable
      autoWrite | Objects set to write at output time
      noWrite   | Objects set to not write by default
      anyWrite  | Any option of the previous two
      log       | Only report registered objects without writing objects
    \endvartable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink

SourceFiles
    writeObjects.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_writeObjects_H
#define Foam_functionObjects_writeObjects_H

#include "functionObject.H"
#include "wordRes.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::writeObjects Declaration
\*---------------------------------------------------------------------------*/

class writeObjects
:
    public functionObject
{
public:

    // Public Data Types

        //- Re-enumeration defining the write options,
        //- Naming based on the IOobjectOption::writeOption
        enum writeOption
        {
            NO_WRITE,
            AUTO_WRITE,
            ANY_WRITE,
            LOG
        };

        //- Names for writeOption
        static const Enum<writeOption> writeOptionNames_;

private:

    // Private Data

        //- Reference to registry
        const objectRegistry& obr_;

        //- To only write objects of defined writeOption
        writeOption writeOption_;

        //- Names of objects to control
        wordRes objectNames_;


    // Private Member Functions

        //- No copy construct
        writeObjects(const writeObjects&) = delete;

        //- No copy assignment
        void operator=(const writeObjects&) = delete;


public:

    //- Runtime type information
    TypeName("writeObjects");


    // Constructors

        //- Construct from Time and dictionary
        writeObjects
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeObjects() = default;


    // Member Functions

        //- Read the writeObjects data
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
