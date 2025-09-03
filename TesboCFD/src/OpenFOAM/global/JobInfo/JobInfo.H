/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::JobInfo

Description
    Helper class for recording information about run/finished jobs,
    acts like global singleton.

    Writes the following files:
      - $FOAM_JOB_DIR/runningJobs
      - $FOAM_JOB_DIR/finishedJobs

    If FOAM_JOB_DIR is unset, defaults to ~/.OpenFOAM/jobControl

Note
    JobInfo is treated as a largely failsafe operation.
    If the directories cannot be written to,
    only a warning (not an error) is emitted.

SourceFiles
    JobInfo.C

\*---------------------------------------------------------------------------*/

#ifndef JobInfo_H
#define JobInfo_H

#include "dictionary.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JobInfo Declaration
\*---------------------------------------------------------------------------*/

class JobInfo
:
    public dictionary
{
    // Private Data

        //- Full path to "running" job information file (master only)
        fileName running_;

        //- Full path to "finished" job information file (master only)
        fileName finished_;

        //- Timing information
        cpuTime cpuTime_;


    // Private Member Functions

        //- Relocate job file from "running" to "finished" directory.
        //  Invalidates job contents, clears 'constructed'.
        void jobEnding();

        //- Add elapsed times, termination type, remove file
        //- from "running" and create in "finished" directory.
        //  Invalidates job contents, clears 'constructed'.
        void jobEnding(const word& terminationType);

        //- No copy construct
        JobInfo(const JobInfo&) = delete;

        //- No copy assignment
        void operator=(const JobInfo&) = delete;


public:

    //- Global value for constructed job info
    static bool constructed;

    //- Global value for writeJobInfo enabled
    static bool writeJobInfo;


    // Constructors

        //- Default construct
        JobInfo();


    //- Destructor, relocates the job file from running to finished.
    ~JobInfo();


    // Static Member Functions

        //- Disallow JobInfo by forcing writeJobInfo (InfoSwitch) off.
        static void disable() noexcept;

        //- Simple shutdown (finalize) of JobInfo
        //  Relocates the job file from "running" to "finished".
        //  Invalidates job contents, clears 'constructed'.
        static void shutdown();

        //- Exit or abort shutdown (finalize) of JobInfo
        static void shutdown(bool isAbort);


    // Member Functions

        //- Write job info to its file in the "running" jobs directory
        void write() const;

        //- Job end with "normal" termination
        void stop();

        //- Job end with "exit" termination
        void exit();

        //- Job end with "abort" termination
        void abort();

        //- Relocate job file from "running" to "finished" directory.
        //  Invalidates job contents, clears 'constructed'.
        void signalEnd();


    // Housekeeping

        //- Same as stop
        void end() { stop(); }
};


// Job information bookkeeping (global)
extern JobInfo jobInfo;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
