/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::ParRunControl

Description
    Helper class for initializing parallel jobs from the command arguments,
    storing 'dry-run' state etc.
    Also handles cleanup of parallel (or serial) jobs.

Note
    In the meantime the class name may be slightly misleading.

\*---------------------------------------------------------------------------*/

#ifndef Foam_argListRunControl_H
#define Foam_argListRunControl_H

#include "UPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParRunControl Declaration
\*---------------------------------------------------------------------------*/

class ParRunControl
{
    // Private Data

        //- The dry-run level
        int dryRun_;

        //- The verbosity level
        int verbose_;

        //- True if this is (or will be) a parallel run
        bool parallel_;

        //- Uses distributed roots
        bool distributed_;

        //- MPI threads are desired
        bool needsThread_;


public:

    // Constructors

        //- Default construct
        ParRunControl() noexcept
        :
            dryRun_(0),
            verbose_(0),
            parallel_(false),
            distributed_(false),
            needsThread_(false)
        {}


    //- Destructor. Shutdown (finalize) MPI as required
    ~ParRunControl()
    {
        if (parallel_)
        {
            Info<< "Finalising parallel run" << endl;
        }
        UPstream::shutdown();
    }


    // Member Functions - General Control

        //- Return the dry-run level (default: 0)
        int dryRun() const noexcept { return dryRun_; }

        //- Increase the dry-run level
        void incrDryRun(int level = 1) noexcept { dryRun_ += level; }

        //- Change dry-run level, returns old value
        int dryRun(const int level) noexcept
        {
            int old(dryRun_);
            dryRun_ = level;
            return old;
        }

        //- Return the verbosity level (default: 0)
        int verbose() const noexcept { return verbose_; }

        //- Increase the verbosity level
        void incrVerbose(int level = 1) noexcept { verbose_ += level; }

        //- Change verbosity level, returns old value
        int verbose(const int level) noexcept
        {
            int old(verbose_);
            verbose_ = level;
            return old;
        }


    // Member Functions - Parallel Control

        //- True if this is (or will be) a parallel run
        bool parRun() const noexcept
        {
            return parallel_;
        }

        //- Set as parallel run on/off, return the previous value.
        //  Use with \b extreme caution if runPar() has already been
        //  called.
        bool parRun(const bool on) noexcept
        {
            bool old(parallel_);
            parallel_ = on;
            return old;
        }

        //- True if a parallel run and uses distributed roots.
        bool distributed() const noexcept
        {
            return (parallel_ && distributed_);
        }

        //- Set use of distributed roots, but only if actually parallel
        void distributed(bool on) noexcept
        {
            distributed_ = (parallel_ && on);
        }

        //- True if MPI threads are desired (default: false)
        bool threads() const noexcept
        {
            return needsThread_;
        }

        //- Set preference for use of MPI threads
        void threads(bool on) noexcept
        {
            needsThread_ = on;
        }

        //- Initialize UPstream for a parallel run
        void runPar(int& argc, char**& argv)
        {
            if (!UPstream::init(argc, argv, needsThread_))
            {
                Info<< "Failed to start parallel run" << endl;
                UPstream::exit(1);
            }
            parallel_ = true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
