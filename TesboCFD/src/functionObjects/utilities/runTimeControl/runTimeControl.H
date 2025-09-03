/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::runTimeControl

Group
    grpUtilitiesFunctionObjects

Description
    Controls when the calculation is terminated based on satisfying
    user-specified conditions.

    Optionally specify a number of write steps before the calculation is
    terminated.  Here, a write is performed each time that all conditions are
    satisfied.

SourceFiles
    runTimeControl.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControl_H
#define functionObjects_runTimeControl_H

#include "fvMeshFunctionObject.H"
#include "Map.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

// Forward declaration of classes
class runTimeCondition;

/*---------------------------------------------------------------------------*\
                       Class runTimeControl Declaration
\*---------------------------------------------------------------------------*/

class runTimeControl
:
    public fvMeshFunctionObject
{
public:

    // Public enumerations

        enum class satisfiedAction
        {
            ABORT,      //!< "abort" - write and emit a FatalError
            END,        //!< "end" - write and exit cleanly
            SET_TRIGGER //!< "setTrigger" - trigger another condition
        };

        static Enum<satisfiedAction> satisfiedActionNames;


private:

    // Private data

        //- List of conditions to satisfy
        PtrList<runTimeCondition> conditions_;

        //- Map to define group IDs
        Map<label> groupMap_;

        //- Number of write steps before exiting
        label nWriteStep_;

        //- Current number of steps written
        label writeStepI_;

        //- Action to take when conditions are satisfied
        satisfiedAction satisfiedAction_;

        //- Trigger index if satisfiedAction is setTrigger
        label triggerIndex_;

        //- Active flag
        //  Used in the trigger case to bypass any evaluations after the
        //  trigger has been set
        bool active_;

        //- Can be restarted flag
        //  Used in the trigger case after the trigger has been set to allow
        //  this object to be restarted/reset the active flag
        bool canRestart_;


    // Private Member Functions

        //- No copy construct
        runTimeControl(const runTimeControl&) = delete;

        //- No copy assignment
        void operator=(const runTimeControl&) = delete;


public:

    //- Runtime type information
    TypeName("runTimeControl");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        runTimeControl
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~runTimeControl() = default;


    // Member Functions

        //- Read the runTimeControl data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the runTimeControl and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
