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
    Foam::SIMPLEControlOpt

Description
    SIMPLE control class for optimisation runs. Each time is sub-cycled and
    corresponds to one optimisation cycle

\*---------------------------------------------------------------------------*/

#ifndef SIMPLEControlOpt_H
#define SIMPLEControlOpt_H

#include "SIMPLEControl.H"
#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SIMPLEControlOpt Declaration
\*---------------------------------------------------------------------------*/

class SIMPLEControlOpt
:
    public SIMPLEControl
{
protected:

    // Protected Data

        //- Number of iterations for the first optimisation cycle
        label nInitialIters_;

        //- Subcycled time
        autoPtr<subCycleTime> subCycledTimePtr_;

        //- deltaT of the sub-cycled time
        scalar deltaTSubSycle_;


    // Protected Member Functions

        //- Read controls from optimisationDict
        virtual bool read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Maximum number of solver iterations
        const label& nIters() const;

        //- Reset deltaT in case controlDict has been re-written at run-time
        void resetDeltaT();


private:

        //- No copy construct
        SIMPLEControlOpt(const SIMPLEControlOpt&) = delete;

        //- No copy assignment
        void operator=(const SIMPLEControlOpt&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("steadyOptimisation");


    // Constructors

        //- Construct from mesh
        SIMPLEControlOpt
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~SIMPLEControlOpt() = default;


    // Member Functions

        // Solution control

            //- Whether to call time.write() or not
            virtual bool write(const bool valid = true) const;

        // Evolution

            //- Loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
