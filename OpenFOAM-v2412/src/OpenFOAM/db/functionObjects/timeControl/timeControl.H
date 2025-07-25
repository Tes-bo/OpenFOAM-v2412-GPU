/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::timeControl

Description
    General time dependent execution controller.
    The execution parameters are given by the "Control" and (optionally)
    the "Interval", with the default being to execute every time-step.

    For example, an "execute" control every tenth write time:
    \verbatim
        executeControl  writeTime;
        executeInterval 10;
    \endverbatim

    See Foam::functionObject for a list of known selection types.

SourceFiles
    timeControl.C

\*---------------------------------------------------------------------------*/

#ifndef timeControl_H
#define timeControl_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeControl Declaration
\*---------------------------------------------------------------------------*/

class timeControl
{
public:

    //- The time control options
    enum timeControls
    {
        ocNone = 0,           //!< No execution
        ocAlways,             //!< Always execute
        ocTimeStep,           //!< Execution coupled to time-step (default)
        ocWriteTime,          //!< Execution coupled to write-time
        ocRunTime,            //!< Use run-time for execution
        ocAdjustableRunTime,  //!< Currently identical to "runTime"
        ocClockTime,          //!< Use clock time for execution
        ocCpuTime,            //!< Use CPU time for execution
        ocOnEnd               //!< Execute on end of run
    };


private:

        //- The selection names for timeControls enumerations
        static const Enum<timeControls> controlNames_;


    // Private Data

        //- Time object
        const Time& time_;

        //- The prefix for the control name (eg, "write", "execute", ...)
        const word prefix_;

        //- Type of time control
        timeControls timeControl_;

        //- Execution interval steps (timeStep | writeTime).
        //  A value <= 1 means execute at every time-step or writeTime.
        label intInterval_;

        //- Execution interval (clockTime | runTime | cpuTime | adjustable)
        scalar interval_;

        //- Index of previous execution
        label executionIndex_;


    // Private Member Functions

        //- No copy construct
        timeControl(const timeControl&) = delete;

        //- No copy assignment
        void operator=(const timeControl&) = delete;


public:

    // Constructors

        //- Construct a control object that executes each time-step.
        //  For places where a time control object is required, but should
        //  not actually intervene.
        explicit timeControl(const Time& runTime, const word& prefix = "");

        //- Construct from Time object, dictionary and the prefix for the
        //- control name.
        timeControl
        (
            const Time& runTime,
            const dictionary& dict,
            const word& prefix
        );


    //- Destructor
    ~timeControl() = default;


    // Member Functions

        //- Identify if a timeControl object is present in the dictionary
        //  Matches prefix + "Control"
        static bool entriesPresent(const dictionary& dict, const word& prefix);

        //- Read from dictionary
        void read(const dictionary& dict);

        //- Return the Time
        inline const Time& time() const;

        //- Return the name (prefix)
        inline const word& name() const;

        //- Return the named control enumeration as its 'type'
        inline const word& type() const;

        //- Reset control to 'always' - ie, no intervention
        void clear();

        //- Return the control enumeration
        inline timeControls control() const;

        //- Return true if the control will always execute - ie, no intervention
        inline bool always() const;

        //- Flag to indicate whether to execute
        bool execute();

        //- Return interval
        inline scalar interval() const;

        //- Return the index of the previous execution
        inline label executionIndex() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
