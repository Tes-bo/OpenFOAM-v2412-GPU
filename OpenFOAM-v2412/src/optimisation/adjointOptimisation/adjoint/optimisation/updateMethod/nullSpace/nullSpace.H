/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 PCOpt/NTUA
    Copyright (C) 2023 FOSS GP
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
    Foam::nullSpace

Description
    Update design variables using a null space approach.
    Can handle inequality and bound constraints.

    Reference:
    \verbatim
        Feppon, F., Allaire, G., & Dapogny, C. (2020).
        Null space gradient flows for constrained optimization with
        applications to shape optimization.
        ESAIM: COCV, 26, 90.
        https://doi.org/10.1051/cocv/2020015
    \endverbatim

SourceFiles
    nullSpace.C

\*---------------------------------------------------------------------------*/

#ifndef nullSpace_H
#define nullSpace_H

#include "constrainedOptimisationMethod.H"
#include "updateMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class nullSpace Declaration
\*---------------------------------------------------------------------------*/

class nullSpace
:
    public constrainedOptimisationMethod,
    public updateMethod
{
protected:

    // Protected data

        //- Are bound constraints included?
        bool includeBoundConstraints_;

        //- Lagrange multipliers for flow-reated constraints
        scalarField mu_;

        //- Lagrange multipliers for the lower bound constraints
        scalarField l_;

        //- Lagrange multipliers for the upper bound constraints
        scalarField u_;


        // Dual variables

            //- Solve the dual problem?
            bool solveDualProblem_;

            //- Lagrange multipliers of the dual problem for flow-related
            //- constraints
            scalarField dualMu_;

            //- Lagrange multipliers of the dual problem for the lower bound
            //- constraints
            scalarField dualL_;

            //- Lagrange multipliers of the dual problem for the upper bound
            //- constraints
            scalarField dualU_;


        // Constraint subsets

            //- List of saturated or violated constraints
            //  List[0]: flow related constraints
            //  List[1]: lower bound constraints
            //  List[2]: upper bound constraints
            //  Sub-list meaning is the same for all labelListLists that follow
            labelListList iTilda_;

            //- List of saturated or violated constraints (up to epsConstr_)
            labelListList iTildaEps_;

            //- List of constraints that must remain active
            //  Determines the null space update
            labelListList iHat_;

            //- List of constraints the values of which need to be reduced
            //  Determines the range space update
            labelListList iRangeSpace_;


        // Fields holding the updates of the Lagrange multipliers of the primal
        // and dual problems

            scalarField deltaMu_;
            scalarField deltaL_;
            scalarField deltaU_;
            scalarField deltaDualMu_;
            scalarField deltaDualL_;
            scalarField deltaDualU_;


        //- Infinitesimal quantity
        //  Updated during the inner iterations of the dual problem
        scalar eps_;

        //- Maxmimum number of Newton iterations for the dual problem
        label maxNewtonIters_;

        //- Maxmimum number of line search iterations for each Newton iteration
        //- of the dual problem
        label maxLineSearchIters_;

        //- Maxmimum number of CG iterations for obtaining the null space and
        //- range space updates
        label maxCGIters_;

        //- Tolerance of the dual problem
        scalar dualTolerance_;


        // Constant parameters

            //- Value for considering a constraint as marginally active
            //  Used to avoid the frequent change of the active set of
            //  constraints
            scalar epsConstr_;

            //- Value to perturb the design variables by, if all of them
            //- lay on their bounds at the beginning of the optimisation
            scalar epsPerturb_;

            //- Multiplier of the null space update
            scalar aJ_;

            //- Multiplier of the range space update
            scalar aC_;

            //- Change aJ and aC adaptively
            bool adaptiveStep_;

            //- Last cycle on which to reset aJ
            label lastAcceleratedCycle_;

            //- Max change of the design variables, pertaining to the objective
            //- reduction.
            //  If adaptiveStep is set to true, aJ will be set in such a way
            //  to obtain a maxDVChange_ for the design variables, for all
            //  optimisation cycles up to lastAcceleratedCycle_
            autoPtr<scalar> maxDVChange_;

            //- Whether to impose maxDVChange strictly on all optimisation
            //- cycles or just up to the lastAcceleratedCycle
            bool strictMaxDVChange_;

            //- Target reduction of active constraints (range in [0-1])
            //  1: active constraints will become zero (first-order approx)
            //  0: no change in the constraints
            scalar targetConstraintReduction_;

            //- Downplay the importance of the bound constraint reduction
            //- by this ammount
            scalar bcMult_;


    // Protected Member Functions

        //- Update sizes of fields related to the constraints
        void initialise();

        //- Zero the updates computed in the previous optimisation cycle
        void zeroUpdates();


        // Functions related to the solution of the dual subproblem

            //- Solve the dual problem for computing the Lagrange multiplier
            //- using a Newton solver
            void solveDualProblem();

            //- Compute and return residuals based on the current solution
            tmp<scalarField> computeResiduals();

            //- Compute direction for the Newton problem
            void computeNewtonDirection();

            //- Perform line search and return max residual corresponding to
            //- the updated solution
            scalar lineSearch();

            //- Update the current solution using the known direction and the
            //- given step length
            void updateSolution(const scalar step);

            //- Adjust step to satisfy cireteria
            void adjustStep
            (
                scalar& step,
                const scalar value,
                const scalar update
            );


        // Functions related to the computation of the design variables update

            //- Update the constraint subsets
            void updateViolatedConstraintsSubsets();

            //- Update the constraint subsets
            void updateNullAndRangeSpaceSubsets();

            //- Update violated constraints indices (iTilda and iTildaEps)
            void updateViolatedIndices
            (
                const label i,
                const scalarField& constraints
            );

            //- Update constraint indices related to the null and range space
            //- part of the correction
            void updateCorrectionIndices
            (
                const label i,
                const scalarField& LagrangeMults,
                const scalarField& dual
            );

            //- Compute the update of the design variables as a combination of
            //- null space and range space parts
            void correction();

            //- A & v
            //  where A is the Jacobian of all constraints, identified by
            //  the labelList addressings, w.r.t. the active design variables
            tmp<scalarField> Av
            (
                const scalarField& v,
                const labelListList& subsets
            );

            //- A.T & v
            //  where A is the Jacobian of all constraints, identified by
            //  the labelList addressings, w.r.t. the active design variables
            tmp<scalarField> ATv
            (
                const scalarField& v,
                const labelListList& subsets
            );

            //- Collect all constraint values corresponding to the provided
            //- indices
            tmp<scalarField> activeConstraints
            (
                const labelListList& subsets
            );

            //- Computes the parts of ksiJ and ksiC that require a system
            //- solution
            //  Formulated in terms of the flow constraints, which are usually
            //  few
            tmp<scalarField> constraintRelatedUpdate
            (
                const scalarField& rhs,
                const labelListList& subsets
            );

            //- Print statistics on the number of flow- and bound-related
            //- constraints included in the subset
            void statistics
            (
                const labelListList& subset,
                const word& description
            );


private:

    // Private Member Functions

        //- No copy construct
        nullSpace(const nullSpace&) = delete;

        //- No copy assignment
        void operator=(const nullSpace&) = delete;


public:

    //- Runtime type information
    TypeName("nullSpace");


    // Constructors

        //- Construct from components
        nullSpace
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~nullSpace() = default;


    // Member Functions

        //- Compute design variables correction
        virtual void computeCorrection();

        //- Compute the merit function for line search
        virtual scalar computeMeritFunction();

        //- Directional derivative of the merit function, in the direction of
        //- the correction
        virtual scalar meritFunctionDirectionalDerivative();

        //- Write useful quantities to files
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
