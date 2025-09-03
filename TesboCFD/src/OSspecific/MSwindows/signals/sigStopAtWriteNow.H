/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::sigStopAtWriteNow

Description
    Signal handler to write and stop the job.
    The interrupt is defined by OptimisationSwitches::stopAtWriteNowSignal

SourceFiles
    sigStopAtWriteNow.C

See also
    Foam::JobInfo

\*---------------------------------------------------------------------------*/

#ifndef Foam_sigStopAtWriteNow_H
#define Foam_sigStopAtWriteNow_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                      Class sigStopAtWriteNow Declaration
\*---------------------------------------------------------------------------*/

class sigStopAtWriteNow
{
    // Private Data

        //- Signal number to use
        static int signal_;


    // Private Member Functions

        //- Handler for caught signals
        static void sigHandler(int);


public:

    //- Allow setter access to signal_
    friend class addstopAtWriteNowSignalToOpt;


    // Constructors

        //- Default construct
        sigStopAtWriteNow() noexcept = default;

        //- Construct with Time reference
        explicit sigStopAtWriteNow(const Time& runTime, bool verbose=false);


    //- Destructor
    ~sigStopAtWriteNow();


    // Member Functions

        //- Is active?
        static bool active() noexcept { return (signal_ > 0); }

        //- Signal number being used
        static int signalNumber() noexcept { return signal_; }

        //- Set/reset signal handler
        static void set(bool verbose=false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
