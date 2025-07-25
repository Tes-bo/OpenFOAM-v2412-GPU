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
    Foam::functionObjects::setTimeStepFaRegionsFunctionObject

Group
    grpUtilitiesFunctionObjects

Description
    This function object controls the time step for classes of the type
    \c regionFaModel. It reads \c regionFaMaxCo entry from \c controlDict
    and evaluate the time step based on the finite area Courant Number.

    Can only be used with solvers using \c adjustTimeStep control (e.g.
    \c pimpleFoam). It makes no attempt to co-operate with other time step
    'controllers', e.g. \c maxCo, other functionObjects. Supports \c enabled
    flag but none of the other options \c timeStart, \c timeEnd, \c writeControl
    etc.

Usage
    Example of function object specification to manipulate the time step:
    \verbatim
    setTimeStep1
    {
        // Mandatory entries
        type        setTimeStepFaRegion;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                    | Type | Reqd | Deflt
      type         | Type name: setTimeStepFaRegion | word | yes  | -
      enabled      | On/off switch                  | bool | no   | yes
    \endtable

    The inherited entries are elaborated in:
     - \link timeFunctionObject.H \endlink
     - \link regionFaModel.H \endlink

SourceFiles
    setTimeStepFaRegionsFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setTimeStepFaRegionsFunctionObject_H
#define functionObjects_setTimeStepFaRegionsFunctionObject_H

#include "timeFunctionObject.H"
#include "regionFaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::regionModels;

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class setTimeStepFaRegionsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setTimeStepFaRegionsFunctionObject
:
    public functionObjects::timeFunctionObject
{
    // Private Member Functions

        //- No copy construct
        setTimeStepFaRegionsFunctionObject
        (
            const setTimeStepFaRegionsFunctionObject&
        ) = delete;

        //- No copy assignment
        void operator=(const setTimeStepFaRegionsFunctionObject&) = delete;

        //- Return minimum deltaT from fa regions
        scalar regionDeltaT() const;


public:

    //- Runtime type information
    TypeName("setTimeStepFaRegion");


    // Constructors

        //- Construct from components
        setTimeStepFaRegionsFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor
    virtual ~setTimeStepFaRegionsFunctionObject() = default;


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
