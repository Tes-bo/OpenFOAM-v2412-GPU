/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 PCOpt/NTUA
    Copyright (C) 2020-2021 FOSS GP
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
    Foam::ISQP

Description
    An L-BFGS-based SQP algorithm for computing the update of the design
    variables in the presence of inequality constraints. The QP problem is
    solved using the interior point method (hence the I in the classe name,
    which is not standard in the terminology used in the literature). The
    (potentially dense) linear system formulated by the interior point method
    is solved using Conjugate Gradient and a choice from three preconditioners
    (diagonal, inverse Hessian, Sherman-Morrison), using matrix-vector products
    to avoid storing the LHS matrix.

    Bound constraints on the design variables will also be included, if set by
    the designVariables known by the updateMethod. If the QP problem is
    infeasible, the algorithm can still be used by setting includeExtraVars_
    to true, to allow a computation of an update, despite not being able to
    satisfy all the constraints of the QP problem.  Alternatively,
    targetConstraintReduction can be set to a number lower than 1 to help with
    the feasibility of the dual problem.


SourceFiles
    ISQP.C

\*---------------------------------------------------------------------------*/

#ifndef ISQP_H
#define ISQP_H

#include "LBFGS.H"
#include "SQPBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class ISQP Declaration
\*---------------------------------------------------------------------------*/

class ISQP
:
    public LBFGS,
    public SQPBase
{
public:

        //- Enumeration of preconditioners for the subproblem
        enum preconditioners
        {
            diag, invHessian, ShermanMorrison
        };

        //- Names of preconditioners for the subproblem
        static const Enum<preconditioners> preconditionerNames;


protected:

    // Protected data

        //- Should Lagrange multipliers be allocated
        bool doAllocateLagrangeMultipliers_;

        //- Are bound constraints included?
        bool includeBoundConstraints_;

        //- Are additional design variables included?
        //  These are introduced to find relaxed solutions, even if the
        //  original problem does not have any feasible points
        bool includeExtraVars_;

        //- The set of design variables being updated during the subproblem.
        //  Size is that of the active design variables.
        //  Correction will end up being the difference of this and the old
        //  design variables.
        scalarField p_;

        // Lagrange multipliers and slack variables

            //- Slack variables for the inequality constraints
            scalarField gs_;

            //- Lagrange multipliers for the lower bound constraints
            autoPtr<scalarField> lTilda_;

            //- Slack variables the lower bound constraints
            autoPtr<scalarField> ls_;

            //- Lagrange multipliers for the upper bound constraints
            autoPtr<scalarField> uTilda_;

            //- Slack variables the upper bound constraints
            autoPtr<scalarField> us_;

            //- Extra variables for finding solutions even in infeasible
            //- problems
            autoPtr<scalarField> extraVars_;

            //- Lagrange multipliers for positive extra variables
            autoPtr<scalarField> z_;

            //- Multiplier of the additional variables y in the Lagrangian, to
            //- make them 'expensive'
            autoPtr<Function1<scalar>> c_;


        // Fields holding updates of the design, Lagrange and slack variables

            scalarField deltaP_;
            scalarField deltaLamda_;
            scalarField deltaGs_;
            autoPtr<scalarField> deltaLTilda_;
            autoPtr<scalarField> deltaLs_;
            autoPtr<scalarField> deltaUTilda_;
            autoPtr<scalarField> deltaUs_;
            autoPtr<scalarField> deltaExtraVars_;
            autoPtr<scalarField> deltaZ_;


        //- Infinitesimal quantity
        //  Updated during the inner iterations of the subproblem
        scalar eps_;

        //- Final eps quantity to be reached during the solution of the
        //- subproblem
        scalar epsMin_;

        //- Maxmimum number of Newton iterations for the subproblem
        label maxNewtonIters_;

        //- Maxmimum number of line search iterations for each iteration of the
        //- subproblem
        label maxLineSearchIters_;

        //- Maxmimum number of iterations for the deltaX equation
        label maxDxIters_;

        //- Relative tolerance of the CG solver, solving for deltaP_
        scalar relTol_;

        //- Which preconditioner to use for the solution of the subproblem
        label preconType_;

        //- Percentage reduction for the constraints.
        //  Can be used to relax QP problems with no feasible points
        scalar cRed_;

        //- Disable damping
        //  Generally discouraged but could be useful in cases with no
        //  constraints or only bound ones
        bool disableDamping_;

        //- File including the l1 merit function
        autoPtr<OFstream> meritFunctionFile_;


    // Protected Member Functions

        //- Update sizes of fields related to the active design variables
        void updateSizes();

        //- Allocate multipliers for the bound constraints
        void allocateBoundMultipliers();

        //- Allocate Lagrange multipliers for the inequality constraints
        void allocateLagrangeMultipliers();

        //- Update the vectors accosiated with the Hessian matrix
        void updateYS();

        //- Allocate fields related to constraints
        void initialize();


        // Functions related to the solution of the primal-dual subproblem

            //- Solve subproblem using a Newton optimiser
            void solveSubproblem();

            //- Compute direction for the Newton problem
            void computeNewtonDirection();

            //- Zero the updates computed in the previous optimisation cycle
            void zeroUpdates();

            //- Solve the equation for deltaX, which is the expensive part of
            //- the Newtopn step.
            //  All other corrections can be computed based on this
            void solveDeltaPEqn();

            //- Compute the RHS for the deltaX equation
            //  Currently not in use
            tmp<scalarField> computeRHSForDeltaX(const scalarField& FDx);

            //- CG algorithm for the solution of the deltaP eqn
            void CGforDeltaP(const scalarField& FDx);

            //- Procudt of the LHS of the deltaP eqn with a vector
            tmp<scalarField> matrixVectorProduct
            (
                const scalarField& vector
            );

            //- Preconditioner-vector product for CG
            //  In case a diagonal preconditioner is used, it is stored in
            //  precon for all CG iterations
            tmp<scalarField> preconVectorProduct
            (
                const scalarField& vector,
                autoPtr<scalarField>& precon
            );

            //- Diagonal preconditioner of the deltaP eqn
            tmp<scalarField> diagPreconditioner();

            //- Use the Sherman-Morrison formula to compute the matrix-free
            //- product of the approximate inverse of the LHS with a vector
            tmp<scalarField> ShermanMorrisonPrecon(const scalarField& vector);

            //- Recursive (naive) implementation of the rank-1 update
            //  WIP, efficient for up to 2 flow-related constraints
            tmp<scalarField> ShermanMorrisonRank1Update
            (
                const PtrList<scalarField>& r1Updates,
                const scalarField& p,
                const refPtr<scalarField>& diag,
                const scalarField& mult,
                label n
            );

            //- Perform line search and return max residual corresponding to
            //- the updated solution
            scalar lineSearch();

            //- Adjust step to satisfy cireteria
            void adjustStep
            (
                scalar& step,
                const scalar value,
                const scalar update
            );

            //- Update the current solution using the known direction and the
            //- given step length
            void updateSolution(const scalar step);


        // Residuals of the various KKT conditions

            //- Compute and return residuals based on the current solution
            tmp<scalarField> computeResiduals();

            //- Residual of the gradient of the Lagrangian
            //  Size is that of the active design variables
            tmp<scalarField> resFL();

            //- Product of the inverse Hessian with the residual of the
            //- gradient of the Lagrangian.
            //  Avoid the formation of the Hessian matrix.
            //  Size is that of the active design variables.
            tmp<scalarField> invHFL();

            //- Residual of the inequality constraint slackness
            tmp<scalarField> resFGs();

            //- Residual of the complementarity slackness for the
            //- inequality constraints
            tmp<scalarField> resFlamda();

            //- Residual of the lower bounds slackness
            tmp<scalarField> resFls();

            //- Residual of the upper bounds slackness
            tmp<scalarField> resFus();

            //- Residual of the complementarity slackness for the
            //- lower bound constraints
            tmp<scalarField> resFlTilda();

            //- Residual of the complementarity slackness for the
            //- upper bound constraints
            tmp<scalarField> resFuTilda();

            //- Residual of the Lagrangian derivative wrt the extra variables
            tmp<scalarField> resFExtraVars();

            //- Residual of the complementarity slackness for the
            //- extra variables
            tmp<scalarField> resFz();


        //- Store old fields needed for the next iter
        void storeOldFields();

        //- Get the part the merit function that depends on the constraints
        virtual scalar meritFunctionConstraintPart() const;


private:

    // Private Member Functions

        //- No copy construct
        ISQP(const ISQP&) = delete;

        //- No copy assignment
        void operator=(const ISQP&) = delete;


public:

    //- Runtime type information
    TypeName("ISQP");


    // Constructors

        //- Construct from components
        ISQP
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~ISQP() = default;


    // Member Functions

        //- Compute design variables correction
        void computeCorrection();

        //- Compute merit function. Could be different than the objective
        //- in the presence of constraints
        virtual scalar computeMeritFunction();

        //- Derivative of the merit function. Could be different than the
        //- objective derivative in the presence of constraints
        virtual scalar meritFunctionDirectionalDerivative();

        //- Update old correction. Useful for quasi-Newton methods coupled with
        //- line search
        virtual void updateOldCorrection(const scalarField& oldCorrection);

        //- Write useful quantities to files
        virtual bool writeData(Ostream& os) const;

        //- Write merit function information
        virtual bool writeAuxiliaryData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
