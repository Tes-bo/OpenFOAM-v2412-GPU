/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjects::runTimeControls::runTimeCondition

Description
    Base class for run time conditions

SourceFiles
    runTimeCondition.C
    runTimeConditionNew.C
    runTimeCondition.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControls_runTimeCondition_H
#define functionObjects_runTimeControls_runTimeCondition_H

#include "stateFunctionObject.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class runTimeCondition Declaration
\*---------------------------------------------------------------------------*/

class runTimeCondition
{

protected:

    // Protected data

        //- Condition name
        word name_;

        //- Reference to the object registry
        const objectRegistry& obr_;

        //- State
        stateFunctionObject& state_;

        //- On/off switch
        bool active_;

        //- Reference to the condition dictionary
        dictionary& conditionDict_;

        //- Group index - if applied, all conditions in a group must be
        //  satisfied before condition is met
        label groupID_;


    // Protected Member Functions

        //- Set the condition dictionary (create if necessary)
        dictionary& setConditionDict();

        //- Return const access to the conditions dictionary
        const dictionary& conditionDict() const;

        //- Return non-const access to the conditions dictionary
        dictionary& conditionDict();


public:

    //- Runtime type information
    TypeName("runTimeCondition");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        runTimeCondition,
        dictionary,
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            stateFunctionObject& state
        ),
        (name, obr, dict, state)
    );


    //- Constructor
    runTimeCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~runTimeCondition() = default;

    //- Selector
    static autoPtr<runTimeCondition> New
    (
        const word& conditionName,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    // Public Data

        //- Switch to send output to Info
        Switch log;


    // Public Member Functions

        //- Return the condition name
        virtual const word& name() const;

        //- Return the active flag
        virtual bool active() const;

        //- Return the group index
        virtual label groupID() const;

        //- Apply the condition
        virtual bool apply() = 0;

        //- Write
        virtual void write() = 0;

        //- Reset
        virtual void reset() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
