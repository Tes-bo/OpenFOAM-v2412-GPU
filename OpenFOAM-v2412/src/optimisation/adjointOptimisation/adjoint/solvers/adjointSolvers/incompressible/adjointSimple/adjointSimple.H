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
    Foam::adjointSimple

Description
    Solution of the adjoint PDEs for incompressible, steady-state flows

    Reference:
    \verbatim
        For the development of the adjoint PDEs

            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255-299.
            http://doi.org/10.1007/s11831-014-9141-9
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef adjointSimple_H
#define adjointSimple_H

#include "incompressibleAdjointSolver.H"
#include "SIMPLEControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class adjointSimple Declaration
\*---------------------------------------------------------------------------*/

class adjointSimple
:
    public incompressibleAdjointSolver
{

private:

    // Private Member Functions

        //- No copy construct
        adjointSimple(const adjointSimple&) = delete;

        //- No copy assignment
        void operator=(const adjointSimple&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<SIMPLEControl> solverControl_;

        //- Reference to incompressibleAdjointVars
        //  Used for convenience and to avoid repetitive dynamic_casts
        //  Same as getAdjointVars()
        incompressibleAdjointVars& adjointVars_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;


    // Protected Member Functions

        //- Allocate incompressibleAdjointVars and return reference to be used
        //- for convenience in the rest of the class.
        incompressibleAdjointVars& allocateVars();

        //- Compute continuity errors
        void continuityErrors();

        //- Accumulate the sensitivities integrand before calculating them
        virtual void preCalculateSensitivities();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointSimple");


    // Constructors

        //- Construct from mesh and dictionary
        adjointSimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    //- Destructor
    virtual ~adjointSimple() = default;


    // Member Functions

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Steps to be executed before each main SIMPLE iteration
            virtual void preIter();

            //- The main SIMPLE iter
            virtual void mainIter();

            //- Steps to be executed before each main SIMPLE iteration
            virtual void postIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Functions to be called before loop
            virtual void preLoop();


            // Source terms to be added to the adjoint equations

                //- Source terms for the momentum equations
                virtual void addMomentumSource(fvVectorMatrix& matrix);

                //- Source terms for the continuity equation
                virtual void addPressureSource(fvScalarMatrix& matrix);


            //- Update primal based quantities
            //- related to the objective functions.
            //  Also writes the objective function values to files.
            //  Written here and not in the postLoop function of the primal
            //  to make sure we don't pollute the objective files with
            //  objectives of non-converged linearSearch iterations
            virtual void updatePrimalBasedQuantities();

            //- Add fvOptions for TopO
            virtual void addTopOFvOptions() const;

            //- Write average iteration
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
