/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::systemCall

Group
    grpUtilitiesFunctionObjects

Description
    Executes system calls, entered in the form of string lists.

    Calls can be made at the following points in the calculation:
    - every time step
    - every output time
    - end of the calculation

Usage
    Example of the function object specification:
    \verbatim
    systemCall1
    {
        type        systemCall;
        libs        (utilityFunctionObjects);
        ...
        executeCalls
        (
            "echo execute"
        );
        writeCalls
        (
            "echo === writing data ==="
        );
        endCalls
        (
            "echo === echoing .bashrc ==="
            "cat ~/.bashrc"
            "echo \*\*\* done \*\*\*"
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description              | Required    | Default value
        type         | type name: systemCall    | yes         |
        executeCalls | list of calls on execute | yes         |
        writeCalls   | list of calls on write   | yes         |
        endCalls     | list of calls on end     | yes         |
        master       | execute on master only   | no          | false
    \endtable

Note
    Since this function object executes system calls, there is a potential
    security risk.  In order to use the \c systemCall function object, the
    \c allowSystemOperations must be set to '1'; otherwise, system calls will
    not be allowed.

    Additionally, since the system commands are normally sent via the shell,
    special shell character may require backslash escaping.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    systemCall.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_systemCall_H
#define functionObjects_systemCall_H

#include "functionObject.H"
#include "stringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class systemCall Declaration
\*---------------------------------------------------------------------------*/

class systemCall
:
    public functionObject
{
protected:

    // Private data

        //- List of calls to execute - every step
        stringList executeCalls_;

        //- List of calls to execute - write steps
        stringList writeCalls_;

        //- List of calls to execute when exiting the time-loop
        stringList endCalls_;

        //- Perform system calls on the master only
        bool masterOnly_;


    // Protected Member Functions

        //- Dispatch specified calls
        label dispatch(const stringList& calls);

        //- No copy construct
        systemCall(const systemCall&) = delete;

        //- No copy assignment
        void operator=(const systemCall&) = delete;


public:

    //- Runtime type information
    TypeName("systemCall");


    // Constructors

        //- Construct from Time and dictionary
        systemCall
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~systemCall() = default;


    // Member Functions

        //- Read the system calls
        virtual bool read(const dictionary& dict);

        //- Execute the "executeCalls" at each time-step
        virtual bool execute();

        //- Write, execute the "writeCalls"
        virtual bool write();

        //- Execute the "endCalls" at the final time-loop
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
