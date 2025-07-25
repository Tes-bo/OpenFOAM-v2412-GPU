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
    Foam::subLoopTime

Description
    A class for managing sub-loops referencing Time.

SourceFiles
    subLoopTime.C

SeeAlso
    subCycleTime

\*---------------------------------------------------------------------------*/

#ifndef subLoopTime_H
#define subLoopTime_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class subLoopTime Declaration
\*---------------------------------------------------------------------------*/

class subLoopTime
{
    // Private Member Functions

        //- No copy construct
        subLoopTime(const subLoopTime&) = delete;

        //- No copy assignment
        void operator=(const subLoopTime&) = delete;


protected:

    // Protected data

        Time& time_;

        //- The current sub-cycle index
        label index_;

        //- The total number of cycles in the sub-cycle
        label total_;


public:

    // Constructors

        //- Construct from time with fixed number of cycles
        //  \param runTime  the top-level time
        //  \param nCycles  the number of times to loop
        subLoopTime(Time& runTime, const label nCycles);


    //- Destructor
    ~subLoopTime();


    // Member Functions

        //- The current cycle index
        inline label index() const
        {
            return index_;
        }

        //- The total number of cycles
        inline label nCycles() const
        {
            return total_;
        }

        //- Return true if looping is active
        bool status() const;

        //- Force terminate the loop
        void stop();

        //- True if looping is active, increments the index.
        //  Example usage,
        //  \code
        //      while (control.loop())
        //      {
        //          solve;
        //      }
        //  \endcode
        bool loop();


    // Member Operators

        //- Prefix increment
        subLoopTime& operator++();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
