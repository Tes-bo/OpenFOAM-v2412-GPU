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
    Foam::optimisationManager

Description
    Abstract base class for optimisation methods

SourceFiles
    optimisationManager.C

\*---------------------------------------------------------------------------*/

#ifndef optimisationManager_H
#define optimisationManager_H

#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "primalSolver.H"
#include "designVariables.H"
#include "designVariablesUpdate.H"
#include "adjointSolverManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class optimisationManager Declaration
\*---------------------------------------------------------------------------*/

class optimisationManager
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Reference to the time
        Time& time_;

        //- Design variables of the optimisation problem.
        //  Constructed before the adjoint solvers, since
        //  some objective functions depend on their values
        autoPtr<designVariables> designVars_;

        //- List of primal solvers to be included in the optimisation
        PtrList<primalSolver> primalSolvers_;

        //- List of adjoint solver managers to be included in the optimisation
        PtrList<adjointSolverManager> adjointSolverManagers_;

        //- Type of the optimisation manager (singleRun, (un)steadyOptimisation)
        const word managerType_;

        //- Helper class managing parts of the optimisation.
        //  Allocated only when an actual optimisation is conducted.
        autoPtr<designVariablesUpdate> dvUpdate_;

        //- Switch defining if the design variables should be updated or not.
        Switch shouldUpdateDesignVariables_;


    // Protected Member Functions

        //- Reset time
        //  Does nothing in base, useful in unsteady adjoint
        virtual void resetTime
        (
            const dimensionedScalar startTime,
            const label startTimeIndex,
            const scalar endTime
        );

        //- Update design variables
        virtual void moveDesignVariables();

        //- Update design variables. Multiplication with line search step
        //  happens inside the update(direction) function
        virtual void moveDesignVariables
        (
            const scalarField& direction
        );

        //- Update design variables using a line-search
        void lineSearchUpdate();

        //- Update design variables using a fixed step
        void fixedStepUpdate();

        //- Initialization. Construct primal and adjoint solvers
        virtual void initialize();


private:

    // Private Member Functions

        //- No copy construct
        optimisationManager(const optimisationManager&) = delete;

        //- No copy assignment
        void operator=(const optimisationManager&) = delete;


public:

    //- Runtime type information
    TypeName("optimisationManager");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            optimisationManager,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        optimisationManager(fvMesh& mesh);


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<optimisationManager> New(fvMesh& mesh);


    //- Destructor
    virtual ~optimisationManager() = default;


    // Member Functions

        // Access

            //- Get the primal solvers
            inline PtrList<primalSolver>& primalSolvers()
            {
                return primalSolvers_;
            }

            //- Get the adjoint solver managers
            inline PtrList<adjointSolverManager>& adjointSolverManagers()
            {
                return adjointSolverManagers_;
            }

            //- Get the design variables
            inline autoPtr<designVariables>& getDesignVariables()
            {
                return designVars_;
            }

            //- Get the mechanism supporting the update of the design variables
            inline autoPtr<designVariablesUpdate>& getOptimisationType()
            {
                return dvUpdate_;
            }


        // Evolution

            //- Changes in case of run-time update of optimisationDict
            virtual bool read();

            //- Prefix increment
            virtual optimisationManager& operator++() = 0;

            //- Postfix increment, this is identical to the prefix increment
            virtual optimisationManager& operator++(int) = 0;

            //- Return true if end of optimisation run.
            //  Also, updates the design variables if needed
            virtual bool checkEndOfLoopAndUpdate() = 0;

            //- Return true if end of optimisation run
            virtual bool end() = 0;

            //- Whether to update the design variables
            virtual bool update() = 0;

            //- Update design variables.
            //  Might employ a line search to find a correction satisfying the
            //  step convergence criteria
            virtual void updateDesignVariables();

            //- Solve all primal equations
            virtual void solvePrimalEquations();

            //- Solve all adjoint equations
            virtual void solveAdjointEquations();

            //- Compute all adjoint sensitivities
            virtual void computeSensitivities();

            //- Clear all adjoint sensitivities
            virtual void clearSensitivities();

            //- Solve all primal equations
            virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
