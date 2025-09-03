/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::sigSegv

Description
    Signal handler for SEGV interrupt.

    The standard interrupt handler is overridden to ensure that the
    runningJob file is removed.
    Can be used either directly through the static member functions or
    through the scope of the object (constructor sets trapping; destructor
    restores original).

See also
    Foam::JobInfo

SourceFiles
    sigSegv.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sigSegv_H
#define Foam_sigSegv_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigSegv Declaration
\*---------------------------------------------------------------------------*/

class sigSegv
{
    // Private Data

        //- Signal trapping enabled?
        static bool sigActive_;


    // Private Member Functions

        //- Handler for caught signals
        static void sigHandler(int);


public:

    // Constructors

        //- Default construct
        sigSegv();


    //- Destructor
    ~sigSegv();


    // Member functions

        //- Activate SIGSEGV signal handler
        static void set(bool verbose=false);

        //- Deactivate SIGSEGV signal handler
        static void unset(bool verbose=false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
