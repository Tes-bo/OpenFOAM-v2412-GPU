/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::functionObjects::parProfiling

Group
    grpUtilitiesFunctionObjects

Description
    Simple (simplistic) mpi-profiling.

Usage
    Example of function object specification:
    \verbatim
    profiling
    {
        type  parProfiling;
        libs  (utilityFunctionObjects);

        // Report stats on exit only (instead of every time step)
        executeControl  onEnd;
        writeControl    none;
        detail          0;
    }
    \endverbatim

SourceFiles
    parProfiling.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_parProfiling_H
#define Foam_functionObjects_parProfiling_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class parProfiling Declaration
\*---------------------------------------------------------------------------*/

class parProfiling
:
    public functionObject
{
    // Private Data

        //- The reporting level
        //  0: summary, 1: per-proc times, 2: per-proc times/counts
        int reportLevel_;

public:

    // Generated Methods

        //- No copy construct
        parProfiling(const parProfiling&) = delete;

        //- No copy assignment
        void operator=(const parProfiling&) = delete;


    //- Runtime type information
    TypeName("parProfiling");


    // Constructors

        //- Construct from Time and dictionary. Enables profilingPstream
        parProfiling
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor. Disables profilingPstream
    virtual ~parProfiling();


    // Member Functions

        //- Report the current profiling information
        void report();

        //- Report
        virtual bool execute();

        //- Do nothing
        virtual bool write();

        //- Disables profilingPstream
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
