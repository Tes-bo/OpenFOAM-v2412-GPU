/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::SIMPLEControl

Description
    SIMPLE control class to supply convergence information/checks for
    the SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef SIMPLEControl_H
#define SIMPLEControl_H

#include "solverControl.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SIMPLEControl Declaration
\*---------------------------------------------------------------------------*/

class SIMPLEControl
:
    public solverControl,
    public simpleControl
{

protected:

    // Protected Data

        //- Optimisation type
        const word& managerType_;

        //- Number of SIMPLE iterations
        label nIters_;

        //- Pressure reference cell
        label  pRefCell_;
        //
        //- Pressure reference value
        scalar pRefValue_;


private:

        //- No copy construct
        SIMPLEControl(const SIMPLEControl&) = delete;

        //- No copy assignment
        void operator=(const SIMPLEControl&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("SIMPLEControl");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            SIMPLEControl,
            dictionary,
            (
                fvMesh& mesh,
                const word& managerType,
                const solver& solver
            ),
            (mesh, managerType, solver)
        );


    // Constructors

        //- Construct from mesh
        SIMPLEControl
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~SIMPLEControl() = default;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<SIMPLEControl> New
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    // Member Functions

        virtual bool read();

        virtual void readIters();

        void checkMeanSolution() const;

        // Access

            //- Return the solution dictionary
            inline virtual const dictionary dict() const;

            //- Return pressure reference cell
            inline label& pRefCell();

            //- Return pressure reference value
            inline scalar& pRefValue();

        // Solution control

            //- Whether to call time.write() or not
            virtual bool write(const bool valid = true) const = 0;

        // Evolution

            //- Loop
            virtual bool loop() = 0;

            //- Is the solver converged
            virtual bool converged();

            //- Perform this iteration?
            virtual bool performIter();

            //- Increment iteration counter
            virtual void incrementIter();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "SIMPLEControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
