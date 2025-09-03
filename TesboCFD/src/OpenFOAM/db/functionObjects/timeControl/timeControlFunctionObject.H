/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::timeControl

Description
    Wrapper around functionObjects to add time control.

    Adds
    - timeStart, timeEnd activation
    - deltaT modification to hit writeInterval (note: if adjustableRunTime
      also in controlDict this has to be an integer multiple)
    - smooth deltaT variation through optional deltaTCoeff parameter

Note
    Since the timeIndex is used directly from Foam::Time, it is unaffected
    by user-time conversions. For example, Foam::engineTime might cause \a
    writeInterval to be degrees crank angle, but the functionObject
    execution \a interval would still be in timestep.

    The function object can be limited to operate within a time range using
    the timeStart and timEnd options.  All objects are read (and the
    OutputFilter allocated) on construction.  However, if a timeEnd is
    supplied, the object will call the 'end' function of the filter
    at the timeEnd time and destroy the filter.
    Any other callback (execute(), write(), timeSet() etc) will only operate
    if within the timeStart, timeEnd time range. Note that the time range
    uses 0.5 * deltaT as a comparison tolerance to account for precision errors.

SourceFiles
    timeControlFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_timeControl_H
#define functionObjects_timeControl_H

#include "timeFunctionObject.H"
#include "dictionary.H"
#include "timeControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::timeControl Declaration
\*---------------------------------------------------------------------------*/

class timeControl
:
    public functionObjects::timeFunctionObject
{
public:

    // Public Enumerations

        //- Control mode
        enum class controlMode
        {
            TIME,
            TRIGGER,
            TIME_OR_TRIGGER,
            TIME_AND_TRIGGER
        };

        static const Enum<controlMode> controlModeNames_;


private:

    // Private Data

        //- Input dictionary
        dictionary dict_;


        // Optional user inputs

            //- Control mode (combination of time/trigger behaviour)
            controlMode controlMode_;

            //- Activation time - defaults to -VGREAT (ie, active)
            scalar timeStart_;

            //- De-activation time - defaults to VGREAT
            scalar timeEnd_;

            //- Activation trigger index - defaults to labelMax (ie, inactive)
            label triggerStart_;

            //- De-activation trigger index - defaults to labelMax
            label triggerEnd_;

            //- Max time step change
            scalar deltaTCoeff_;

            //- Number of steps before the dump-time during which deltaT
            //  may be changed (valid for adjustableRunTime)
            label nStepsToStartTimeChange_;


        //- Execute controls
        Foam::timeControl executeControl_;

        //- Write controls
        Foam::timeControl writeControl_;

        //- The functionObject to execute
        autoPtr<functionObject> foPtr_;

        //- Time index of the last execute call
        label executeTimeIndex_;


        // For time-step change limiting (deltaTCoeff supplied) only:

            //- Store the old deltaT value
            scalar deltaT0_;

            //- Store the series expansion coefficient value
            scalar seriesDTCoeff_;



    // Private Member Functions

        //- Read relevant dictionary entries
        void readControls();

        //- Returns true if within time bounds
        bool active() const;

        //- Return expansion ratio (deltaT change) that gives overall time
        static scalar calcExpansion
        (
            const scalar startRatio,
            const scalar y,
            const label n
        );

        //- Calculate deltaT change such that the next write interval is
        //  obeyed. Updates requiredDeltaTCoeff
        void calcDeltaTCoeff
        (
            scalar& requiredDeltaTCoeff,
            const scalar wantedDT,
            const label nSteps,
            const scalar presentTime,
            const scalar timeToNextWrite,
            const bool rampDirectionUp
        );

        //- No copy construct
        timeControl(const timeControl&) = delete;

        //- No copy assignment
        void operator=(const timeControl&) = delete;


public:

    //- Runtime type information
    TypeName("timeControl");


    // Constructors

        //- Construct from components
        timeControl
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return the input dictionary
            inline const dictionary& dict() const;

            //- Return the execute control object
            inline const Foam::timeControl& executeControl() const;

            //- Return the write control object
            inline const Foam::timeControl& writeControl() const;

            //- Return the functionObject filter
            inline const functionObject& filter() const;


        // Function object control

            //- Helper function to identify if a timeControl object is present
            //  in the dictionary
            static bool entriesPresent(const dictionary& dict);

            //- Called at each ++ or += of the time-loop.
            //  postProcess overrides the usual executeControl behaviour and
            //  forces execution (used in post-processing mode)
            virtual bool execute();

            //- Execute using the specified subIndex.
            virtual bool execute(const label subIndex);

            //- Called at each ++ or += of the time-loop.
            //  postProcess overrides the usual writeControl behaviour and
            //  forces writing (used in post-processing mode)
            virtual bool write();

            //- Called when Time::run() determines that the time-loop exits
            virtual bool end();

            //- Called at the end of Time::adjustDeltaT() if adjustTime is true
            virtual bool adjustTimeStep();

            //- Did any file get changed during execution?
            virtual bool filesModified() const;

            //- Read and set the function object if its data have changed
            virtual bool read(const dictionary&);

            //- Update for changes of mesh
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Update for changes of mesh
            virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeControlFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
