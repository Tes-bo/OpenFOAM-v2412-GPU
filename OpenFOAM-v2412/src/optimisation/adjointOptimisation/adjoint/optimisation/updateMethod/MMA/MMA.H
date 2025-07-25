/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::MMA

Description
    Update design variables using the Method of Moving Assymptotes.
    Can handle inequality constraints.

    Reference:
    \verbatim
        Svanberg, K. (1987)
        The method of moving asymptotes—a new method for structural
        optimization
        International Journal for Numerical Methods in Engineering, 24(2),
        359-373.
        https://doi.org/10.1002/nme.1620240207
    \endverbatim

SourceFiles
    MMA.C

\*---------------------------------------------------------------------------*/

#ifndef MMA_H
#define MMA_H

#include "constrainedOptimisationMethod.H"
#include "updateMethod.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class MMA Declaration
\*---------------------------------------------------------------------------*/

class MMA
:
    public constrainedOptimisationMethod,
    public updateMethod
{
protected:

    // Protected data

        //- The previous design variables. Used to compute new bounds
        scalarField x0_;

        //- The twice previous design variables. Used to compute new bounds
        scalarField x00_;

        //- The set of design variables being updated during the subproblem.
        //  Correction will end up being the difference of this and the
        //  old values of the design variables
        scalarField xNew_;

        //- Old objective value
        //  Needed by GCMMA
        scalar oldObjectiveValue_;

        //- Old constraint values
        //  Needed by GCMMA
        scalarField oldCValues_;

        //- Objective/Constraint values and approximations in the new point
        //  Needed by GCMMA
        PtrList<scalarField> valsAndApproxs_;

        //- Second term in the approximation function
        scalar z_;

        //- Term multiplying z in the objective function
        scalar alpha0_;

        //- Terms multyplying z in the constraint functions
        scalarField alpha_;

        //- y in the constraints functions
        scalarField y_;

        //- y multipliers in the objective function
        autoPtr<Function1<scalar>> c_;

        //- y^2 multipliers in the objective function
        scalarField d_;

        //- Lower design variables asymptotes
        scalarField lower_;

        //- Upper design variables asymptotes
        scalarField upper_;

        //- Lower design variables constraints
        scalarField a_;

        //- Upper design variables constraints
        scalarField b_;

        //- Constants in the (p,q) functions.
        //  Size is cValues_.size() + 1 (the objective)
        scalarField rho_;

        //- Bound the rho value with an upper value?
        bool boundRho_;

        //- Correct the design variables
        bool correctDVs_;

        // Dual variables

            //- Constraint Lagrange multipliers
            scalarField lamda_;

            //- Lagrange multipliers for the lower limits constraints
            scalarField ksi_;

            //- Lagrange multipliers for the upper limits constraints
            scalarField Eta_;

            //- Lagrange multipliers for the y constraints
            scalarField mu_;

            //- Lagrange multiplier for the z constraint
            scalar zeta_;

            //- Slack variables for the inequality constraints
            scalarField s_;


        // Fields holding updates of the design, Lagrange and slack variables

            scalarField deltaLamda_;
            scalar deltaZ_;
            scalarField deltaX_;
            scalarField deltaY_;
            scalarField deltaS_;
            scalar deltaZeta_;
            scalarField deltaMu_;
            scalarField deltaEta_;
            scalarField deltaKsi_;


        //- Infinitesimal quantity
        //  Updated during the inner iterations of the subproblem
        scalar eps_;

        //- Maxmimum number of Newton iterations for the subproblem
        label maxNewtonIters_;

        //- Maxmimum number of line search iterations for each iteration of the
        //- subproblem
        label maxLineSearchIters_;

        //- Initialize every subproblem with x = 0.5*(a + b)
        //- and reset Lagrange multipliers
        //  The solution of the subproblem can become relatively expensive
        //  in some cases, hence the GCMMA subproblems are usually
        //  initialized with the solution of the previous inner iteration,
        //  to save up some time. Occassionally, this can cause the solution
        //  of the subproblem to be trapped in a non-feasible point.
        //  In such cases, the solution of each subproblem could be
        //  initialized from scratch
        bool initializeEverySubproblem_;


        // Constant parameters

            //- Lower assymprote reduction rate
            scalar asymDecr_;

            //- Upper assymprote increase rate
            scalar asymIncr_;

            //- Used to update the assymptotes in the first optimisation cycle
            scalar sInit_;

            //- Movement of the unatainable upper and lower bounds
            scalar move_;

            //- Constant in p, q functions
            scalar raa0_;


        // Parameters related to the update of rho in GCMMA

            //- Multiplier of the mean derivative used for the initialisation
            //- of rho
            scalar maxInitRhoMult_;

            //- Multiplier for the max. rho value during its update
            scalar maxRhoMult_;

            //- Change rho constants in each iteration?
            bool variableRho_;

            //- Change rho constants in each iteration?
            bool dynamicRhoInitialisation_;

            //- The rho values obtained by the dynamic rho initialisation
            //- might be too conservative. Multiply with this number to relax
            //- them
            scalar dynamicRhoMult_;


        // Normalisation based on the gradient

            //- Perform the normalisation
            bool normalise_;

            //- Perform the normalisation in each optimisation cycle
            bool continuousNormalisation_;

            //- Normalisation factor for the objective
            autoPtr<scalar> Jnorm_;

            //- Normalisation factors for the constraints
            tmp<scalarField> Cnorm_;

            //- Constaint weight after the normalisation
            scalar cw_;

            //- Constaint weight after the normalisation
            label lastNormalisationStep_;


    // Protected Member Functions

        //- Update sizes of fields related to the active design variables
        void updateSizes();

        //- Initialize rho constants in the (p, q) functions
        //  These control how aggressive or conservative the method will be
        void initializeRho();

        //- Update the bounds associated with the design variables
        void updateBounds();

        //- Allocate fields related to constraints
        void initialize();

        //- Store old objcective and constraint values
        //  Needed by GCMMA
        void storeOldValues();

        //- p and q functions, used to approximate the objective and contraints
        //  Computed based on the current set of design variables, not updated
        //  during the iterations of the subproblem
        tmp<scalarField> p
        (
            const scalarField& derivs,
            const scalar r,
            const scalarField& x
        );
        tmp<scalarField> p(const scalarField& derivs, const scalar r);
        tmp<scalarField> q
        (
            const scalarField& derivs,
            const scalar r,
            const scalarField& x
        );
        tmp<scalarField> q(const scalarField& derivs, const scalar r);

        //- g of all constraint functions
        //  Computed using the current set of design variables, updated during
        //  the solution of the subproblem
        tmp<scalarField> gConstr(const scalarField& vars);

        //- The rhs of the contraint approximations.
        //  Always computed using the current (frozen) set of design variables
        tmp<scalarField> b();

        //- p and q functions, using the dual lamda
        tmp<scalarField> pLamda();
        tmp<scalarField> qLamda();

        //- Zero the updates computed in the previous optimisation cycle
        void zeroUpdates();

        //- Computation of the approximate function
        scalar fTilda
        (
            const scalarField& xInit,
            const scalarField& x,
            const scalar f,
            const scalarField& dfdx,
            const scalar rho
        );

        //- Read in scalarField from self (potentially in binary), or allocate
        //- field with the size of the active design variables and given value
        tmp<scalarField> getOrDefaultScalarField
        (
            const word& name,
            const label size,
            const scalar value = Zero
        );

        //- Read in scalarField from self (potentially in binary), or allocate
        //- field with the size of the active design variables and given value
        void setOrDefaultScalarField
        (
            scalarField& field,
            const word& name,
            const label size,
            const scalarField& defaultField
        );


        // Functions related to the solution of the primal-dual subproblem

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

            //- Compute and return residuals based on the current solution
            tmp<scalarField> computeResiduals();


        //- Normalise the objective and constraints if necessary
        void normalise();


private:

    // Private Member Functions

        //- No copy construct
        MMA(const MMA&);

        //- No copy assignment
        void operator=(const MMA&);


public:

    //- Runtime type information
    TypeName("MMA");


    // Constructors

        //- Construct from components
        MMA
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~MMA() = default;


    // Member Functions

        //- Compute design variables correction
        void computeCorrection();

        //- Solve subproblem using a Newton optimiser
        void solveSubproblem();

        // Functions needed by GCMMA

            //- Compute objective/constraint values and their approximations
            //  Needed to check convergence and update rho, if necessary
            void updateValuesAndApproximations();

            //- Get objective/constraint values and their approximations
            const PtrList<scalarField>& getValuesAndApproximations() const;

            //- Update rho value if needed
            void updateRho();

            //- Get rho value if needed
            const scalarField& getRho() const;

            //- Set variable rho
            void setVariableRho(bool varRho = true);

            //- Checks convergence of GCMMA
            bool converged();


        //- Write useful quantities to files
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
