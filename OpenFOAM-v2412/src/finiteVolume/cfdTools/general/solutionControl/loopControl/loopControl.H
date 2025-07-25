/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::loopControl

Description
    A class for managing arbitrary loops with the ability to invoke
    function object execution.

Usage
    Examples of function object specification:
    \verbatim
    SIMPLE
    {
        energyCoupling
        {
            iterations  100;
            onLoop      ();
            onConverged ( externalCoupled  "loopThings.*" );

            convergence
            {
                "h"     1e-3;
            }
        }
    }
    \endverbatim

    Where the loop entries comprise:
    \table
        Property    | Description           | Required      | Default
        enabled     | activate/deactivate loop | no            | true
        iteration   | times to loop         | no            | 0
        timeStart   | begin time for loop activation  | no  | -VGREAT
        timeEnd     | end time of loop activation     | no  | VGREAT
        interval    | sub-interval to execute onLoop  | no  | 0
        onLoop      | function object names to call at executeInterval | no
        onConverged | function object names to call when converged | no
        onEnd       | function object names to call when loop ends | no
        convergence | dictionary of convergence values to check | no
    \endtable

    The function object names listed by \c onLoop, \c onConverged, \c onEnd
    must implement an \c execute(int) method.
    If the time controls \c timeStart or \c timeEnd are used for the loop,
    these values are only inspected upon creation, not during execution.

SeeAlso
    fvSolution

SourceFiles
    loopControl.C

\*---------------------------------------------------------------------------*/

#ifndef loopControl_H
#define loopControl_H

#include "subLoopTime.H"
#include "dictionary.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class loopControl Declaration
\*---------------------------------------------------------------------------*/

class loopControl
:
    public subLoopTime
{
    // Private Member Functions

        //- Reset
        void clear();

        //- Read settings from dictionary
        void read(const dictionary& dict);

        //- Execute specified function names
        bool checkConverged() const;

        //- No copy construct
        loopControl(const loopControl&) = delete;

        //- No copy assignment
        void operator=(const loopControl&) = delete;


protected:

    // Protected data

        //- Name of the loop control (the lookup dictionary name).
        word name_;

        //- The interval to execute onLoop function-objects
        label interval_;

        //- Dictionary for checking convergence (all regions)
        dictionary convergenceDict_;

        //- Function object names to fire during the loop (at executeInterval)
        wordRes onLoop_;

        //- Function object names to fire on convergence
        wordRes onConverged_;

        //- Function object names to fire when the loop exits without
        //- convergence
        wordRes onEnd_;

        //- Convergence tests passed
        bool converged_;

public:

    // Constructors

        //- Construct from time with fixed number of cycles
        //  \param runTime  the top-level time
        //  \param nCycles  the number of times to loop
        //  \param loopName  the name of the loop
        loopControl
        (
            Time& runTime,
            const label nCycles,
            const word& dictName = "loop"
        );

        //- Construct from fvSolution dictionary based on time and the name
        //- of the controlling algorithm
        //  \param runTime  the top-level time
        //  \param algorithmName the name of the fvSolution dictionary,
        //      typically PIMPLE or SIMPLE
        //  \param dictName  the name of the control dictionary
        loopControl
        (
            Time& runTime,
            const word& algorithmName,
            const word& dictName = "loop"
        );

        //- Construct from fvSolution dictionary based on time and the name
        //- of the controlling algorithm
        //  \param runTime  the top-level time
        //  \param algorithmDict the fvSolution algorithm dictionary,
        //      typically PIMPLE or SIMPLE
        //  \param dictName  the name of the control dictionary
        loopControl
        (
            Time& runTime,
            const dictionary& algorithmDict,
            const word& dictName = "loop"
        );


    //- Destructor
    ~loopControl();


    // Member Functions

        //- Name of the loop control
        inline const word& name() const
        {
            return name_;
        }

        //- The interval to execute onLoop function-objects
        inline label interval() const
        {
            return interval_;
        }

        //- True if looping is active, increments the index and executes
        //- the onLoop and onConverged functions.
        //  Example usage,
        //  \code
        //      while (control.loop())
        //      {
        //          solve;
        //      }
        //  \endcode
        bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write name and state (on/off, index/total) to Ostream
Ostream& operator<<(Ostream& os, const loopControl& ctrl);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
