/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::SIMPLEControlSingleRun

Description
    SIMPLE control class for single runs (i.e. not optimisation).
    Time acts as in simpleFoam, with all solver control read through
    optimisationDict

\*---------------------------------------------------------------------------*/

#ifndef SIMPLEControlSingleRun_H
#define SIMPLEControlSingleRun_H

#include "SIMPLEControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SIMPLEControlSingleRun Declaration
\*---------------------------------------------------------------------------*/

class SIMPLEControlSingleRun
:
    public SIMPLEControl
{
protected:

    // Protected Data

        //- Start time
        scalar startTime_;

        //- End time
        scalar endTime_;


    // Protected Member Functions

        //- Read controls from optimisationDict
        virtual bool read();

        //- Set end time if number of iters has changed
        virtual void readIters();

        //- Check whether endTime has been overwritten by a run-time
        //- modification of controlDict
        void checkEndTime(bool& isRunning);


private:

        //- No copy construct
        SIMPLEControlSingleRun(const SIMPLEControlSingleRun&) = delete;

        //- No copy assignment
        void operator=(const SIMPLEControlSingleRun&) = delete;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("singleRun");


    // Constructors

        //- Construct from mesh
        SIMPLEControlSingleRun
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~SIMPLEControlSingleRun() = default;


    // Member Functions

        // Solution control

            //- Whether to call time.write() or not
            virtual bool write(const bool valid = true) const;

            //- Write fields, even if it is not a writeTime
            void writeNow();

        // Evolution

            //- Loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
