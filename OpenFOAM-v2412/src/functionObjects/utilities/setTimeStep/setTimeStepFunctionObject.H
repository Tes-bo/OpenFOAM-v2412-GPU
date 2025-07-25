/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::setTimeStepFunctionObject

Group
    grpUtilitiesFunctionObjects

Description
    This function object overrides the calculation time step.

    Can only be used with solvers with adjustTimeStep control (e.g.
    pimpleFoam).  It makes no attempt to co-operate with other time step
    'controllers', e.g. maxCo, other functionObjects. Supports 'enabled'
    flag but none of the other options 'timeStart', 'timeEnd', 'writeControl'
    etc.

Usage
    Example of function object specification to manipulate the time step:
    \verbatim
    setTimeStep1
    {
        type        setTimeStep;
        libs        (utilityFunctionObjects);

        deltaT      table ((0 5e-4)(0.01 1e-3));

        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: setTimeStep  | yes         |
        enabled      | On/off switch           | no          | yes
        deltaT       | Time step value         | yes         |
        timeStart    | Start time              | no          | 0
        timeEnd      | End time                | no          | GREAT
        deltaTCoeff  | Time step change limit  | no          | GREAT
    \endtable

SourceFiles
    setTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setTimeStepFunctionObject_H
#define functionObjects_setTimeStepFunctionObject_H

#include "timeFunctionObject.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class setTimeStepFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setTimeStepFunctionObject
:
    public functionObjects::timeFunctionObject
{
    // Private Member Data

        //- Time step function/table
        autoPtr<Function1<scalar>> timeStepPtr_;


    // Private Member Functions

        //- No copy construct
        setTimeStepFunctionObject(const setTimeStepFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const setTimeStepFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("setTimeStep");


    // Constructors

        //- Construct from components
        setTimeStepFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor
    virtual ~setTimeStepFunctionObject() = default;


    // Member Functions

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Execute does nothing
        virtual bool execute();

        //- Write does nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
