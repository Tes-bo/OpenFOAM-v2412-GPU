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

\*---------------------------------------------------------------------------*/

#include "MMA.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MMA, 2);
    addToRunTimeSelectionTable
    (
        updateMethod,
        MMA,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        MMA,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::MMA::updateSizes()
{
    const label n = activeDesignVars_.size();
    if (n != xNew_.size())
    {
        x0_.setSize(n, Zero);
        x00_.setSize(n, Zero);
        xNew_.setSize(n, Zero);
        lower_.setSize(n, Zero);
        upper_.setSize(n, Zero);
        a_.setSize(n, Zero);
        b_.setSize(n, Zero);
        ksi_.setSize(n, Zero);
        Eta_.setSize(n, Zero);
        deltaX_.setSize(n, Zero);
        deltaEta_.setSize(n, Zero);
        deltaKsi_.setSize(n, Zero);
    }
}


void Foam::MMA::initializeRho()
{
    rho_.setSize(cValues_.size() + 1, raa0_);

    // Store old objective and constraint values for GCMMA
    oldObjectiveValue_ = objectiveValue_;
    oldCValues_ = cValues_;

    if (variableRho_)
    {
        const scalarField x(designVars_().getVars(), activeDesignVars_);
        const scalarField xMin
            (designVars_().lowerBounds()(), activeDesignVars_);
        const scalarField xMax
            (designVars_().upperBounds()(), activeDesignVars_);
        const scalarField span(xMax - xMin);
        const scalarField activeObjDerivs
            (objectiveDerivatives_, activeDesignVars_);

        // Objective r
        rho_[0] =
            maxInitRhoMult_
           /globalSum(x.size())*globalSum(mag(activeObjDerivs)*span);

        // Constraints r
        forAll(constraintDerivatives_, i)
        {
            const scalarField activeDerivs
            (
                constraintDerivatives_[i], activeDesignVars_
            );
            rho_[i + 1] =
                 maxInitRhoMult_
                /globalSum(x.size())*globalSum(mag(activeDerivs)*span);

        }

        // Find a rho value that produces an approximate objective/constraint
        // combination, evaluated at the new point, that is lower than the
        // current one
        if (dynamicRhoInitialisation_)
        {
            Info<< "-----------------------------------------" << endl;
            Info<< "Solving sub problems for initializing rho" << endl;
            Info<< "-----------------------------------------" << endl;
            // Backup bounds and initialisation, to be restored after the rho
            // initialisation phase
            scalarField lowerBck = lower_;
            scalarField upperBck = upper_;
            scalarField aBck = a_;
            scalarField bBck = b_;
            scalarField x0Bck = x0_;
            scalarField x00Bck = x00_;

            // First, do a peudo-update of the the design variables
            designVars_().storeDesignVariables();
            x0_.map(designVars_().getVars(), activeDesignVars_);

            Info<< "Initial pseudo update of the design variables" << endl;
            updateBounds();
            initialize();
            solveSubproblem();

            designVars_().update(correction_);

            // Check if the approximate function is higher than the CURRENT
            // function and, if not, update the rho value until this happens.
            // Used as an inexpensive way to obtain a decent rho initialisation
            label iter(0);
            while (!converged() && iter++ < 10)
            {
                Info<< nl << "Dynamic rho initialisation iteration " << iter
                    << nl << endl;
                designVars_().resetDesignVariables();
                updateRho();
                solveSubproblem();
                designVars_().update(correction_);
            }

            Info<< "-----------------------------------------" << endl;
            Info<< "Dynamic rho initialisation converged in " << iter
                << " iterations " << endl;
            Info<< "-----------------------------------------" << endl;

            // Restore bounds and design variables
            lower_ = lowerBck;
            upper_ = upperBck;
            a_ = aBck;
            b_ = bBck;
            designVars_().resetDesignVariables();
            x0_ = x0Bck;
            x00_ = x00Bck;
        }

        rho_ *= dynamicRhoMult_;

        // Bound too small values
        if (boundRho_)
        {
            rho_ = max(rho_, scalar(1.e-6));
        }
        DebugInfo
            << "Computed r values " << rho_ << endl;
    }
}


void Foam::MMA::updateBounds()
{
    const scalarField x(designVars_().getVars(), activeDesignVars_);
    const scalarField xMin(designVars_().lowerBoundsRef(), activeDesignVars_);
    const scalarField xMax(designVars_().upperBoundsRef(), activeDesignVars_);
    const scalarField span(xMax - xMin);

    if (counter_ > 1)
    {
        forAll(x, i)
        {
            scalar gamma
            (
                (x[i] - x0_[i])*(x0_[i] - x00_[i]) < scalar(0)
              ? asymDecr_ : asymIncr_
            );
            lower_[i] = x[i] - gamma*(x0_[i] - lower_[i]);
            upper_[i] = x[i] + gamma*(upper_[i] - x0_[i]);
        }

        lower_ = min(lower_, x - 0.01*(xMax - xMin));
        lower_ = max(lower_, x - 10.0*(xMax - xMin));

        upper_ = max(upper_, x + 0.01*(xMax - xMin));
        upper_ = min(upper_, x + 10.0*(xMax - xMin));
    }
    else
    {
        lower_ =  x - sInit_*span;
        upper_ =  x + sInit_*span;
    }

    a_ = max(xMin, lower_ + 0.1*(x - lower_));
    a_ = max(a_, x - move_*span);
    b_ = min(xMax, upper_ - 0.1*(upper_ - x));
    b_ = min(b_, x + move_*span);
}


void Foam::MMA::initialize()
{
    const label m(cValues_.size());
    // Allocate correct sizes for array depending on the number of constraints
    if (!c_.good())
    {
        alpha_.setSize(m, Zero);
        if (coeffsDict(typeName).found("c"))
        {
            Info<< "Reading constraint penalty function type from dict" << endl;
            c_.reset(Function1<scalar>::New("c", coeffsDict(typeName)));
        }
        else
        {
            Info<< "Setting constant penalty factor" << endl;
            c_.reset(new Function1Types::Constant<scalar>("c", 100));
        }
        d_.setSize(m, coeffsDict(typeName).getOrDefault<scalar>("d", 1));
        deltaLamda_.setSize(m, Zero);
        deltaY_.setSize(m, Zero);
        deltaS_.setSize(m, Zero);
        deltaMu_.setSize(m, Zero);
    }
    const scalar t = mesh_.time().timeOutputValue();
    Info<< "Penalty multiplier (c):: " << c_->value(t) << endl;

    // Scalar quantities
    eps_ = 1.;
    z_ = 1.;
    zeta_ = 1.;

    // Quantities related to design variables
    xNew_ = 0.5*(a_ + b_);
    ksi_ = max(scalar(1), 1./(xNew_ - a_));
    Eta_ = max(scalar(1), 1./(b_ - xNew_));

    // Quantities related to constraints
    y_.setSize(m, scalar(1));
    lamda_.setSize(m, scalar(1));
    s_.setSize(m, scalar(1));
    mu_.setSize(m, Zero);
    mu_ = max(scalar(1), 0.5*c_->value(t));
}


void Foam::MMA::storeOldValues()
{
    // Store old solutions
    if (counter_ > 0)
    {
        x00_ = x0_;
    }
    x0_.map(designVars_().getVars(), activeDesignVars_);
}


Foam::tmp<Foam::scalarField> Foam::MMA::p
(
    const scalarField& derivs,
    const scalar r,
    const scalarField& x
)
{
    const scalarField lowerBounds
    (
        designVars_().lowerBoundsRef(), activeDesignVars_
    );
    const scalarField upperBounds
    (
        designVars_().upperBoundsRef(), activeDesignVars_
    );
    tmp<scalarField> tres(tmp<scalarField>::New(x.size(), Zero));
    scalarField& res = tres.ref();

    res =
        sqr(upper_ - x)
       *(
            1.001*max(derivs, scalar(0))
          + 0.001*max(-derivs, scalar(0))
          + r/(upperBounds - lowerBounds)
        );

    return tres;
}

Foam::tmp<Foam::scalarField> Foam::MMA::p
(
    const scalarField& derivs,
    const scalar r
)
{
    return
        p(derivs, r, scalarField(designVars_().getVars(), activeDesignVars_));
}


Foam::tmp<Foam::scalarField> Foam::MMA::q
(
    const scalarField& derivs,
    const scalar r,
    const scalarField& x
)
{
    const scalarField lowerBounds
    (
        designVars_().lowerBoundsRef(), activeDesignVars_
    );
    const scalarField upperBounds
    (
        designVars_().upperBoundsRef(), activeDesignVars_
    );
    tmp<scalarField> tres(tmp<scalarField>::New(x.size(), Zero));
    scalarField& res = tres.ref();

    res =
        sqr(x - lower_)
       *(
            0.001*max(derivs, scalar(0))
          + 1.001*max(-derivs, scalar(0))
          + r/(upperBounds - lowerBounds)
        );

    return tres;
}

Foam::tmp<Foam::scalarField> Foam::MMA::q
(
    const scalarField& derivs,
    const scalar r
)
{
    return
        q(derivs, r, scalarField(designVars_().getVars(), activeDesignVars_));
}


Foam::tmp<Foam::scalarField> Foam::MMA::gConstr
(
    const scalarField& vars
)
{
    tmp<scalarField> tres(tmp<scalarField>::New(cValues_.size(), Zero));
    scalarField& res = tres.ref();
    forAll(res, i)
    {
        const scalarField activeDerivs
        (
            constraintDerivatives_[i], activeDesignVars_
        );
        tmp<scalarField> pI(p(activeDerivs, rho_[i + 1]));
        tmp<scalarField> qI(q(activeDerivs, rho_[i + 1]));
        res[i] = globalSum(pI/(upper_ - vars) + qI/(vars - lower_));
    }

    return tres;
}


Foam::tmp<Foam::scalarField> Foam::MMA::b()
{
    const scalarField vars(designVars_().getVars(), activeDesignVars_);
    tmp<scalarField> tres( - cValues_ + gConstr(vars));

    return tres;
}


Foam::tmp<Foam::scalarField> Foam::MMA::pLamda()
{
    const scalarField activeObjDerivs(objectiveDerivatives_, activeDesignVars_);
    tmp<scalarField> tres(p(activeObjDerivs, rho_[0]));
    scalarField& res = tres.ref();
    forAll(constraintDerivatives_, cI)
    {
        const scalarField activeDerivs
        (
            constraintDerivatives_[cI], activeDesignVars_
        );
        res += lamda_[cI]*p(activeDerivs, rho_[cI + 1]);
    }

    return tres;
}


Foam::tmp<Foam::scalarField> Foam::MMA::qLamda()
{
    const scalarField activeObjDerivs(objectiveDerivatives_, activeDesignVars_);
    tmp<scalarField> tres(q(activeObjDerivs, rho_[0]));
    scalarField& res = tres.ref();
    forAll(constraintDerivatives_, cI)
    {
        const scalarField activeDerivs
        (
            constraintDerivatives_[cI], activeDesignVars_
        );
        res += lamda_[cI]*q(activeDerivs, rho_[cI + 1]);
    }

    return tres;
}


void Foam::MMA::zeroUpdates()
{
    deltaLamda_ = Zero;
    deltaZ_ = Zero;
    deltaX_ = Zero;
    deltaY_ = Zero;
    deltaS_ = Zero;
    deltaZeta_ = Zero;
    deltaMu_ = Zero;
    deltaEta_ = Zero;
    deltaKsi_ = Zero;
}


Foam::scalar Foam::MMA::fTilda
(
    const scalarField& xInit,
    const scalarField& x,
    const scalar f,
    const scalarField& dfdx,
    const scalar rho
)
{
    return scalar
        (
            f
          + globalSum
            (
                p(dfdx, rho, xInit)*(1./(upper_ - x) - 1./(upper_ - xInit))
              + q(dfdx, rho, xInit)*(1./(x - lower_) - 1./(xInit - lower_))
            )
        );
}


Foam::tmp<Foam::scalarField> Foam::MMA::getOrDefaultScalarField
(
    const word& name,
    const label size,
    const scalar value
)
{
    tmp<scalarField> tfield(tmp<scalarField>::New(size, value));
    scalarField& field = tfield.ref();
    if (found(name))
    {
        field = scalarField(name, *this, field.size());
    }
    return tfield;
}


void Foam::MMA::setOrDefaultScalarField
(
    scalarField& field,
    const word& name,
    const label size,
    const scalarField& defaultField
)
{
    if (found(name))
    {
        field = scalarField(name, *this, size);
    }
    else
    {
        field = defaultField;
    }
}


void Foam::MMA::computeNewtonDirection()
{
    // Zero the updates computed in the previous optimisation cycle
    zeroUpdates();

    // Easy access to number of design variables and constraints
    const label n(xNew_.size());
    const label m(cValues_.size());
    const scalarField x(designVars_().getVars(), activeDesignVars_);

    // Fields needed to compute Psi and its derivatives
    tmp<scalarField> tpL(pLamda());
    const scalarField& pL = tpL();
    tmp<scalarField> tqL(qLamda());
    const scalarField& qL = tqL();

    // Build parts of the LHS matrices
    scalarField DYInv(scalar(1.)/(d_ + mu_/y_));
    scalarField DLamdaY((s_/lamda_) + DYInv);
    scalarField DxInv
    (
        scalar(1.)/
        (
            scalar(2)*pL/pow3(upper_ - xNew_)
          + scalar(2)*qL/pow3(xNew_ - lower_)
          + ksi_/(xNew_ - a_)
          + Eta_/(b_ - xNew_)
        )
    );

    // Compute G
    scalarRectangularMatrix G(m, n, Zero);
    for (label i = 0; i < m; ++i)
    {
        const scalarField activeDerivs
        (
            constraintDerivatives_[i], activeDesignVars_
        );
        tmp<scalarField> tp = p(activeDerivs, rho_[i + 1]);
        tmp<scalarField> tq = q(activeDerivs, rho_[i + 1]);
        scalarField row(tp/sqr(upper_ - xNew_) - tq/sqr(xNew_ - lower_));
        forAll(row, j)
        {
            G[i][j] = row[j];
        }
    }

    // Compute parts of the RHS
    scalarField deltaXTilda
    (
        pL/sqr(upper_ - xNew_) - qL/sqr(xNew_ - lower_)
      - eps_/(xNew_ - a_)
      + eps_/(b_ - xNew_)
    );
    const scalar t = mesh_.time().timeOutputValue();
    scalarField deltaYTilda
    (
        c_->value(t) + d_*y_ - lamda_ - eps_/y_
    );
    scalar deltaZTilda = alpha0_ - sum(lamda_*alpha_) - eps_/z_;
    scalarField deltaLamdaYTilda
    (
        // part of deltaLamdaTilda
        gConstr(xNew_)
      - alpha_*z_
      - y_
      - b()
      + eps_/lamda_
      // part of deltaYTilda
      + DYInv*deltaYTilda
    );


    // Solve system for deltaLamda and deltaZ
    // The deltaZ equation is substituted into
    // the one of deltaLamda
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarSquareMatrix lhs(m, m, Zero);
    scalar zRatio(z_/zeta_);
    scalarField rhs(m, Zero);

    // Deal with parts of lhs and rhs that might need global reduction first
    for (label i = 0; i < m; ++i)
    {
        // Prepare LHS
        for (label j = 0; j < m; ++j)
        {
            for (label l = 0; l < n; ++l)
            {
                // Part from GDxInvGT
                lhs[i][j] += DxInv[l]*G[i][l]*G[j][l];
            }
        }
        // Prepare RHS
        for (label l = 0; l < n; ++l)
        {
            // Part from GDxInvDeltaXTilda
            rhs[i] -= DxInv[l]*G[i][l]*deltaXTilda[l];
        }
    }
    // If design variables are distributed to many processors, reduce
    if (globalSum_)
    {
        reduce(lhs, sumOp<scalarSquareMatrix>());
        Pstream::listCombineAllGather(rhs, plusEqOp<scalar>());
    }

    // Add remaining parts from deltaLamdaYTilda and the deltaZ eqn
    rhs += deltaLamdaYTilda + alpha_*deltaZTilda*zRatio;
    for (label i = 0; i < m; ++i)
    {
        // Prepare LHS
        for (label j = 0; j < m; ++j)
        {
            // Part from delta z
            lhs[i][j] += alpha_[i]*alpha_[j]*zRatio;
        }
        // Part from DLamdaY
        lhs[i][i] += DLamdaY[i];
    }
    solve(deltaLamda_, lhs, rhs);
    deltaZ_ = (sum(alpha_*deltaLamda_) - deltaZTilda)*zRatio;

    // Compute the rest of the corrections using backwards substitution

    // deltaX
    // No need for reduction since DxInv is a diagonal matrix
    deltaX_ = - DxInv*(G.Tmul(deltaLamda_) + deltaXTilda);

    // deltaY
    deltaY_ = DYInv*(deltaLamda_ - deltaYTilda);

    // deltaS
    deltaS_ = (-s_*deltaLamda_ + eps_)/lamda_ - s_;

    // deltaZeta
    deltaZeta_ = -zeta_/z_*deltaZ_ - zeta_ + eps_/z_;

    // deltaMu
    deltaMu_ = (-mu_*deltaY_ + eps_)/y_ - mu_;

    // deltaEta
    deltaEta_ = (Eta_*deltaX_ + eps_)/(b_ - xNew_) - Eta_;

    // deltaKsi
    deltaKsi_ = (-ksi_*deltaX_ + eps_)/(xNew_ - a_) - ksi_;
}


Foam::scalar Foam::MMA::lineSearch()
{
    const label n(xNew_.size());
    const label m(cValues_.size());
    scalar step(1.);

    // Perform bound checks and adjust step accordingly
    for (label i = 0; i < n; ++i)
    {
        if
        (
            xNew_[i] + step*deltaX_[i] - a_[i] - 0.01*(xNew_[i] - a_[i])
          < scalar(0)
        )
        {
            step = -0.99*(xNew_[i] - a_[i])/deltaX_[i];
        }

        if
        (
            -xNew_[i] - step*deltaX_[i] + b_[i] - 0.01*(b_[i] - xNew_[i])
          < scalar(0)
        )
        {
            step = 0.99*(b_[i] - xNew_[i])/deltaX_[i];
        }

        adjustStep(step, ksi_[i], deltaKsi_[i]);
        adjustStep(step, Eta_[i], deltaEta_[i]);
    }

    for (label i = 0; i < m; ++i)
    {
        adjustStep(step, y_[i], deltaY_[i]);
        adjustStep(step, lamda_[i], deltaLamda_[i]);
        adjustStep(step, mu_[i], deltaMu_[i]);
        adjustStep(step, s_[i], deltaS_[i]);
    }

    adjustStep(step, z_, deltaZ_);
    adjustStep(step, zeta_, deltaZeta_);

    // Each processor might have computed a different step, if design variables
    // are distributed. Get the global minimum
    if (globalSum_)
    {
        reduce(step, minOp<scalar>());
    }

    if (debug > 1)
    {
        Info<< "Step before line search is " << step << endl;
    }

    // Old residual
    scalar normResOld = sqrt(globalSum(magSqr(computeResiduals())));
    scalar maxRes(GREAT);

    for (label i = 0; i < maxLineSearchIters_ ; ++i)
    {
        // Update the solution with given step
        updateSolution(step);

        // Compute new residuals and their max value
        scalarField resNew(computeResiduals());
        scalar normResNew  = sqrt(globalSum(magSqr(resNew)));
        maxRes = gMax(mag(resNew));

        if (normResNew < normResOld)
        {
            DebugInfo
                << "Initial residual = " << normResOld << ", "
                << "Final residual = " << normResNew << ", "
                << "No of LineSearch Iterations = " << i + 1
                << endl;
            break;
        }
        else
        {
            // Return solution to previous and reduce step
            updateSolution(-step);
            step *= 0.5;

            // If line search could find an appropriate step, increase eps
            if (i == maxLineSearchIters_ - 1)
            {
                eps_ *= 1.5;
                Info<< "Line search could not find a step that reduced"
                    << " residuals while satisfying the constraints" << nl
                    << "Increasing eps to " << eps_ << endl;
            }
        }
    }

    if (debug > 1)
    {
        Info<< "Step after line search is " << step << nl <<  endl;
    }

    return maxRes;
}


void Foam::MMA::adjustStep
(
    scalar& step,
    const scalar value,
    const scalar update
)
{
    if (0.99*value + step*update < scalar(0))
    {
        step = -0.99*value/update;
    }
}


void Foam::MMA::updateSolution(const scalar step)
{
    xNew_ += step*deltaX_;
    y_ += step*deltaY_;
    z_ += step*deltaZ_;
    lamda_ += step*deltaLamda_;
    ksi_ += step*deltaKsi_;
    Eta_ += step*deltaEta_;
    mu_ += step*deltaMu_;
    zeta_ += step*deltaZeta_;
    s_ += step*deltaS_;
}


Foam::tmp<Foam::scalarField> Foam::MMA::computeResiduals()
{
    const label n(xNew_.size());
    const label m(cValues_.size());
    tmp<scalarField> tres(tmp<scalarField>::New(3*n + 4*m + 2, Zero));
    scalarField& res = tres.ref();

    label iRes(0);

    // dLdx
    scalarField dPsidx
       (pLamda()/sqr(upper_ - xNew_) - qLamda()/sqr(xNew_ - lower_));
    for (label i = 0; i < n; ++i)
    {
        res[iRes++] = dPsidx[i] - ksi_[i] + Eta_[i];
    }

    // dLdy
    const scalar t = mesh_.time().timeOutputValue();
    for (label i = 0; i < m; ++i)
    {
        res[iRes++] = c_->value(t) + d_[i]*y_[i] - lamda_[i] - mu_[i];
    }

    // dLdz
    res[iRes++] = alpha0_ - zeta_ - sum(lamda_*alpha_);

    // Primal feasibility (constraints)
    scalarField gb(gConstr(xNew_) - b());
    for (label i = 0; i < m; ++i)
    {
        res[iRes++] = gb[i] - alpha_[i]*z_ - y_[i] + s_[i];
    }

    // ksi slackness
    for (label i = 0; i < n; ++i)
    {
        res[iRes++] = ksi_[i]*(xNew_[i] - a_[i]) - eps_;
    }

    // Eta slackness
    for (label i = 0; i < n; ++i)
    {
        res[iRes++] = Eta_[i]*(b_[i] - xNew_[i]) - eps_;
    }

    // y slackness
    for (label i = 0; i < m; ++i)
    {
        res[iRes++] = mu_[i]*y_[i] - eps_;
    }

    // z slackness
    res[iRes++] = z_*zeta_ - eps_;

    // Constraint slackness
    for (label i = 0; i < m; ++i)
    {
        res[iRes++] = lamda_[i]*s_[i] - eps_;
    }

    return tres;
}


void Foam::MMA::normalise()
{
    if (normalise_)
    {
        // Compute normalisation factors
        if
        (
            !Jnorm_
         || (continuousNormalisation_ && counter_ < lastNormalisationStep_)
        )
        {
            scalarField activeSens(objectiveDerivatives_, activeDesignVars_);
            Jnorm_.reset(new scalar(Foam::sqrt(globalSum(sqr(activeSens)))));
            Cnorm_.reset(new scalarField(cValues_.size(), Zero));
            scalarField& Cnorm = Cnorm_.ref();
            forAll(constraintDerivatives_, cI)
            {
                scalarField activeConstrSens
                    (constraintDerivatives_[cI], activeDesignVars_);
                Cnorm[cI] = Foam::sqrt(globalSum(sqr(activeConstrSens)));
            }
            Info<< "Computed normalisation factors " << nl
                << "J: " << Jnorm_() << nl
                << "C: " << Cnorm_() << endl;
        }

        // Normalise objective values and gradients
        objectiveValue_ /= Jnorm_() + SMALL;
        objectiveDerivatives_ /= Jnorm_() + SMALL;

        // Normalise constraint values and gradients
        cValues_ *= cw_/(Cnorm_() + SMALL);
        forAll(constraintDerivatives_, cI)
        {
            constraintDerivatives_[cI] *= cw_/(Cnorm_()[cI] + SMALL);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MMA::MMA
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    constrainedOptimisationMethod(mesh, dict, designVars, nConstraints, type),
    updateMethod(mesh, dict, designVars, nConstraints, type),
    x0_(getOrDefaultScalarField("x0", activeDesignVars_.size())),
    x00_(getOrDefaultScalarField("x00", activeDesignVars_.size())),
    xNew_(activeDesignVars_.size()),
    oldObjectiveValue_(getOrDefault<scalar>("J0", Zero)),
    oldCValues_(getOrDefaultScalarField("C0", nConstraints_)),
    valsAndApproxs_(2),
    z_(coeffsDict(type).getOrDefault<scalar>("z", 1.)),
    alpha0_(coeffsDict(type).getOrDefault<scalar>("alpha0", 1.)),
    alpha_(0),
    y_(0),
    c_(nullptr),
    d_(0),
    lower_(xNew_.size(), -GREAT),
    upper_(xNew_.size(),  GREAT),
    a_(xNew_.size()),
    b_(xNew_.size()),
    rho_(0),
    boundRho_(coeffsDict(type).getOrDefault<bool>("boundRho", false)),
    correctDVs_(coeffsDict(type).getOrDefault<bool>("correct", true)),
    lamda_(0),
    ksi_(xNew_.size()),
    Eta_(xNew_.size()),
    mu_(0),
    zeta_(coeffsDict(type).getOrDefault<scalar>("zeta", 1.)),
    s_(0),
    deltaLamda_(0),
    deltaZ_(Zero),
    deltaX_(xNew_.size(), Zero),
    deltaY_(0),
    deltaS_(0),
    deltaZeta_(Zero),
    deltaMu_(0),
    deltaEta_(xNew_.size(), Zero),
    deltaKsi_(xNew_.size(), Zero),
    eps_(1),
    maxNewtonIters_(coeffsDict(type).getOrDefault<label>("maxIters", 6000)),
    maxLineSearchIters_
    (
        coeffsDict(type).getOrDefault<label>("maxLineSearchIters", 10)
    ),
    initializeEverySubproblem_
        (coeffsDict(type).getOrDefault<bool>("initializeEverySubproblem", false)),
    asymDecr_(coeffsDict(type).getOrDefault<scalar>("asymptoteDecrease", 0.7)),
    asymIncr_(coeffsDict(type).getOrDefault<scalar>("asymptoteIncrease", 1.2)),
    sInit_(coeffsDict(type).getOrDefault<scalar>("sInit", 0.5)),
    move_(coeffsDict(type).getOrDefault<scalar>("move", 0.5)),
    raa0_(coeffsDict(type).getOrDefault<scalar>("raa0", 1.e-05)),
    maxInitRhoMult_(coeffsDict(type).getOrDefault<scalar>("maxInitRhoMult", 0.1)),
    maxRhoMult_(coeffsDict(type).getOrDefault<scalar>("maxRhoMult", 10)),
    variableRho_(coeffsDict(type).getOrDefault<bool>("variableRho", false)),
    dynamicRhoInitialisation_
        (coeffsDict(type).getOrDefault<bool>("dynamicRhoInitialisation", false)),
    dynamicRhoMult_
        (coeffsDict(type).getOrDefault<scalar>("dynamicRhoMult", 0.1)),
    normalise_(coeffsDict(type).getOrDefault<bool>("normalise", false)),
    continuousNormalisation_
        (coeffsDict(type).getOrDefault<bool>("continuousNormalisation", false)),
    Jnorm_(nullptr),
    Cnorm_(nullptr),
    cw_(coeffsDict(type).getOrDefault<scalar>("constraintWeight", 1)),
    lastNormalisationStep_
    (
        coeffsDict(type).getOrDefault<label>("lastNormalisationStep", 20)
    )
{
    // Check that the design variables bounds have been set
    if (!designVars().lowerBounds() || !designVars().upperBounds())
    {
        FatalErrorInFunction
            << "Lower and upper design variable bounds must be set for MMA"
            << abort(FatalError);
    }

    // Set initial bounds
    setOrDefaultScalarField
    (
        lower_,
        "lower",
        activeDesignVars_.size(),
        scalarField(designVars().lowerBoundsRef(), activeDesignVars_)
    );
    setOrDefaultScalarField
    (
        upper_,
        "upper",
        activeDesignVars_.size(),
        scalarField(designVars().upperBoundsRef(), activeDesignVars_)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MMA::computeCorrection()
{
    if (correctDVs_)
    {
        // Update sizes of fields related to the active design variables
        updateSizes();

        // Perform normalisation of the objective and constraint values
        normalise();

        // Initialize rho values in the p,q functions
        // These determine how aggressive or conservative the method will be
        initializeRho();

        // Update the bounds associated with the design variables
        updateBounds();

        // Initiaze variables
        initialize();

        // Solve subproblem using a Newton optimiser
        solveSubproblem();

        // Store old solutions, bjective and constaint values
        storeOldValues();

        // Increase counter
        ++counter_;
    }
}


void Foam::MMA::solveSubproblem()
{
    // Initialize the solution
    if (initializeEverySubproblem_)
    {
        initialize();
    }

    scalar resMax(gMax(mag(computeResiduals())));
    label iter(0);
    do
    {
        DebugInfo
            << "Newton iter " << iter << nl << endl;

        // Decrease eps
        if (resMax < 0.9*eps_)
        {
            eps_ *= 0.1;
        }

        // Computes Newton direction for the subproblem
        computeNewtonDirection();

        // Perform line search and return max residual of the solution
        // satisfying the bound constraints and the residual reduction.
        // Upates solution.
        resMax = lineSearch();
        DebugInfo
            << "max residual = " << resMax << ", "
            << "eps = " << eps_ << nl << endl;

        mesh_.time().printExecutionTime(Info);
    } while (iter++ < maxNewtonIters_ && (eps_ > 1.e-07 || resMax > 0.9*eps_));

    Info<< "Solved the MMA Newton problem in " << iter << " iterations "
        << nl << endl;

    DebugInfo
        << "Constraint penalization variables (y) " << y_ << endl;

    // Pass update to correction field
    const scalarField& oldVars = designVars_().getVars();
    forAll(activeDesignVars_, avI)
    {
        const label vI(activeDesignVars_[avI]);
        correction_[vI] = xNew_[avI] - oldVars[vI];
    }
}


void Foam::MMA::updateValuesAndApproximations()
{
    // Return values, including the objective and constraint values and
    // their approximations
    label m(cValues_.size());
    valsAndApproxs_.set(0, new scalarField(m + 1, Zero));
    valsAndApproxs_.set(1, new scalarField(m + 1, Zero));
    scalarField& vals = valsAndApproxs_[0];
    scalarField& approx = valsAndApproxs_[1];

    // Objective value and approximation
    const scalarField activeObjDerivs(objectiveDerivatives_, activeDesignVars_);
    vals[0] = objectiveValue_;
    approx[0] = fTilda(x0_, xNew_, oldObjectiveValue_, activeObjDerivs, rho_[0]);

    // Constraint values and approximations
    forAll(constraintDerivatives_, i)
    {
        const scalarField activeDerivs
        (
            constraintDerivatives_[i], activeDesignVars_
        );
        vals[i + 1] = cValues_[i];
        approx[i + 1] =
            fTilda(x0_, xNew_, oldCValues_[i], activeDerivs, rho_[i + 1]);
    }
}


const Foam::PtrList<Foam::scalarField>&
Foam::MMA::getValuesAndApproximations() const
{
    return valsAndApproxs_;
}


void Foam::MMA::updateRho()
{
    // Objective/constraint values and approximations
    const scalarField& vals = valsAndApproxs_[0];
    const scalarField& approx = valsAndApproxs_[1];

    const scalarField xMin(designVars_().lowerBounds()(), activeDesignVars_);
    const scalarField xMax(designVars_().upperBounds()(), activeDesignVars_);

    scalar d = globalSum
        (
            (upper_ - lower_)*sqr(xNew_ - x0_)
           /(upper_ - xNew_)/(xNew_ - lower_)/(xMax - xMin)
        );

    // Update rho values
    forAll(approx, i)
    {
        const scalar delta = (vals[i] - approx[i])/d;
        if (delta > 0)
        {
            rho_[i] = min(1.1*(rho_[i] + delta), maxRhoMult_*rho_[i]);
        }
    }
    DebugInfo
        << "Updated rho values to " << rho_ << endl;
}


const Foam::scalarField& Foam::MMA::getRho() const
{
    return rho_;
}


void Foam::MMA::setVariableRho(bool varRho)
{
    variableRho_ = varRho;
}


bool Foam::MMA::converged()
{
    updateValuesAndApproximations();
    const scalarField& vals = valsAndApproxs_[0];
    const scalarField& approx = valsAndApproxs_[1];

    bool isConverged(true);
    forAll(vals, i)
    {
        DebugInfo
            << nl << "MMA, objective/constraint " << i << nl
            << "Approximation " << approx[i]
            <<  " | old value " << vals[i] << nl << endl;
        isConverged = isConverged && (approx[i] - vals[i] > 0.);
    }

    return isConverged;
}


bool Foam::MMA::writeData(Ostream& os) const
{
    x0_.writeEntry("x0", os);
    x00_.writeEntry("x00", os);
    lower_.writeEntry("lower", os);
    upper_.writeEntry("upper", os);
    os.writeEntry("J0", oldObjectiveValue_);
    oldCValues_.writeEntry("C0", os);

    return updateMethod::writeData(os);
}


// ************************************************************************* //
