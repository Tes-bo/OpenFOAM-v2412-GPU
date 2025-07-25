/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::abort

Group
    grpUtilitiesFunctionObjects

Description
    Watches for presence of the named trigger file in the case directory
    and signals a simulation stop (or other) event if found.

    The presence of the trigger file is only checked on the master process.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

    Example of function object specification:
    \verbatim
    abort
    {
        type    abort;
        libs    (utilityFunctionObjects);

        file    "<case>/GOODBYE";
        action  writeNow
    }
    \endverbatim

    \heading Function object usage
    \table
        Property  | Description                         | Required | Default
        type      | Type name: abort                    | yes |
        file      | The trigger filename                | no  | \<case\>/name
        action    | The default action to trigger       | no  | nextWrite
    \endtable

    When the trigger file is found, it is checked for the following
    content which corresponds to actions.

    - \c action=noWriteNow
       : triggers Foam::Time::saNoWriteNow (stop without writing data)
    - \c action=writeNow
       : triggers Foam::Time::saWriteNow (stop and write data)
    - \c action=nextWrite
       : triggers Foam::Time::saNextWrite (stop after next normal data write)
    - \c action=endTime
       : triggers Foam::Time::saEndTime (continue simulation to the end)
    - Anything else (empty file, no action=, ...)
       : use the default action
    .

Note
    The trigger file is considered "sticky". This means that once detected
    and processed, the trigger is duly noted and the file will not be
    rechecked. It is not possible or desirable to 'untrigger' an action.

SourceFiles
    abort.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_abort_H
#define functionObjects_abort_H

#include "timeFunctionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class abort Declaration
\*---------------------------------------------------------------------------*/

class abort
:
    public functionObjects::timeFunctionObject
{
    // Private Data

        //- The fully-qualified name of the trigger file
        fileName file_;

        //- The default action (defined in dictionary)
        Time::stopAtControls defaultAction_;

        //- Only trigger the action once
        bool triggered_;


    // Private Member Functions

        //- No copy construct
        abort(const abort&) = delete;

        //- No copy assignment
        void operator=(const abort&) = delete;


public:

    //- Runtime type information
    TypeName("abort");


    // Constructors

        //- Construct from Time and dictionary
        abort(const word& name, const Time& runTime, const dictionary& dict);


    //- Destructor
    virtual ~abort() = default;


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Check existence of the file and take action
        virtual bool execute();

        //- No-op
        virtual bool write();

        //- Remove the trigger file after the final time-loop.
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
