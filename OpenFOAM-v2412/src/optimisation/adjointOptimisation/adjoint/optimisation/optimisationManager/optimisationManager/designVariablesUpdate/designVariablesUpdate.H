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
    Foam::designVariablesUpdate

Description
    A class encapsulating functionality neccessary to perform an optimisation
    loop, such as updating the design variables,  checking the
    sufficient reduction/adhetion of objective and constraint values in each
    optimisation cycle by performing a line-search and checking the overall
    convergence of the optimisation loop, if the corresponding criteria are
    provided.

    Kept separate from optimisationManager to isolate functionality required
    only when the update of the design variables is performed.

SourceFiles
    designVariablesUpdate.C

\*---------------------------------------------------------------------------*/

#ifndef designVariablesUpdate_H
#define designVariablesUpdate_H

#include "adjointSolverManager.H"
#include "designVariables.H"
#include "updateMethod.H"
#include "lineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class designVariablesUpdate Declaration
\*---------------------------------------------------------------------------*/

class designVariablesUpdate
{
protected:

    // Protected data

        fvMesh& mesh_;
        const dictionary dict_;
        PtrList<adjointSolverManager>& adjointSolvManagers_;
        autoPtr<designVariables>& designVars_;

        //- Method to update the design variables based on the provided
        //- sensitivity derivatives
        autoPtr<updateMethod> updateMethod_;

        //- Line search mechanism to approximate the update step length
        autoPtr<lineSearch> lineSearch_;

        // Variables related to the computation of the CPU cost

            //- Output file
            OFstream CPUcostFile_;

            //- Primal solutions per optimisation cycle
            label nPrimalsPerCycle_;

            //- Adjoint solutions per optimisation cycle
            label nAdjointsPerCycle_;

            //- Primal evaluations performed so far
            label nPrimalSolutions_;

            //- Adjoint evaluations performed so far
            label nAdjointSolutions_;

            //- CPU cost (in seconds)
            scalar CPUcost_;

        //  Convergence criteria

            //- The maximum of the correction/designVariables values
            //- must be lower that this threshold to consider the run converged
            autoPtr<scalar> designVarsThreshold_;

            //- The relative update of the objective value w.r.t. to its last
            //- value should be smaller than this value to considered the run
            //- converged
            autoPtr<scalar> objectiveThreshold_;

            //- Is at least a single convergence criterion defined
            bool convergenceCriteriaDefined_;

            //- In case of a constrained optimisation, the sum of positive
            //- constraints should be lower than this value to consider the
            //- run converged (i.e. this tolerates some deviation from
            //- satisfying all constraints)
            scalar feasibilityThreshold_;


    // Protected Member Functions

        //- Get the number of adjoint solvers that correspond to constraints
        label nConstraints
        (
            PtrList<adjointSolverManager>& adjointSolverManagers
        ) const;

        //- Get total number of adjoint solvers
        label nAdjointSolvers() const;

        //- Write CPU cost header
        void writeCPUcostHeader();

        //- Write to cost file
        void writeToCostFile(bool zeroAdjointSolns = false);


private:

    // Private Member Functions

        //- No copy construct
        designVariablesUpdate(const designVariablesUpdate&) = delete;

        //- No copy assignment
        void operator=(const designVariablesUpdate&) = delete;


public:

    //- Runtime type information
    TypeName("designVariablesUpdate");


    // Constructors

        //- Construct from components
        designVariablesUpdate
        (
            fvMesh& mesh,
            const dictionary& dict,
            PtrList<adjointSolverManager>& adjointSolverManagers,
            autoPtr<designVariables>& designVars
        );

    //- Destructor
    virtual ~designVariablesUpdate() = default;


    // Member Functions

        //- Update design variables
        void update();

        //- Update design variables based on a given correction
        void update(const scalarField& correction);

        //- Compute update direction
        tmp<scalarField> computeDirection();

        //- Compute cumulative objective and constraint gradients
        void updateGradientsAndValues();

        //- Sum objective values from all adjointSolverManagers
        scalar computeObjectiveValue();

        //- Set the old objective value known by the updateMethod
        //  Used to check convergence of the optimisation loop
        void setOldObjectiveValue();

        //- Compute the merit function of the optimisation problem.
        //  Could be different than the objective function in case of
        //  constraint optimisation
        scalar computeMeritFunction();

        //- Derivative of the merit function
        scalar meritFunctionDirectionalDerivative();

        //- Update old correction. Needed for quasi-Newton Methods
        void updateOldCorrection(const scalarField&);

        //- Write useful quantities to files
        void write();

        //- Steps to be executed after the susccessfull update of the design
        //- varibles, i.e. the last step of line search or the simple update
        //- in the fixedStep approach
        void postUpdate(const scalarField& oldCorrection);

        //- Check if the optimisation loop has converged based on the provided
        //- criteria
        //  May terminate the program
        void checkConvergence();

        //- Get access to design variables
        inline autoPtr<designVariables>& getDesignVariables()
        {
            return designVars_;
        }

        //- Get a reference to the line search object
        inline autoPtr<lineSearch>& getLineSearch()
        {
            return lineSearch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
