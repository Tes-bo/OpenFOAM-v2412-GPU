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

#include "ISQP.H"
#include "IOmanip.H"
#include "Constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ISQP, 2);
    addToRunTimeSelectionTable
    (
        updateMethod,
        ISQP,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        ISQP,
        dictionary
    );
}


const Foam::Enum<Foam::ISQP::preconditioners>
Foam::ISQP::preconditionerNames
({
    { preconditioners::diag, "diag" },
    { preconditioners::invHessian, "invHessian" },
    { preconditioners::ShermanMorrison, "ShermanMorrison" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ISQP::updateSizes()
{
    const label n = activeDesignVars_.size();
    if (n != deltaP_.size())
    {
        // Correction fields
        p_.setSize(n, Zero);
        deltaP_.setSize(n, Zero);

        // Lagrange multipliers and slack variables for bound constraints
        if (includeBoundConstraints_)
        {
            lTilda_().setSize(n, Zero);
            ls_().setSize(n, Zero);
            uTilda_().setSize(n, Zero);
            us_().setSize(n, Zero);

            deltaLTilda_().setSize(n, Zero);
            deltaLs_().setSize(n, Zero);
            deltaUTilda_().setSize(n, Zero);
            deltaUs_().setSize(n, Zero);
        }

        // Fields used to compute the Hessian
        for (label i = 0; i < nPrevSteps_; ++i)
        {
            y_[i].setSize(n, Zero);
            s_[i].setSize(n, Zero);
        }
    }
}


void Foam::ISQP::allocateBoundMultipliers()
{
    if (includeBoundConstraints_)
    {
        // Number of constraints
        const label n(activeDesignVars_.size());

        if (!lTilda_)
        {
            lTilda_.reset(autoPtr<scalarField>::New(n, Zero));
        }
        ls_.reset(autoPtr<scalarField>::New(n, Zero));
        if (!uTilda_)
        {
            uTilda_.reset(autoPtr<scalarField>::New(n, Zero));
        }
        us_.reset(autoPtr<scalarField>::New(n, Zero));

        deltaLTilda_.reset(autoPtr<scalarField>::New(n, Zero));
        deltaLs_.reset(autoPtr<scalarField>::New(n, Zero));
        deltaUTilda_.reset(autoPtr<scalarField>::New(n, Zero));
        deltaUs_.reset(autoPtr<scalarField>::New(n, Zero));
    }
}


void Foam::ISQP::allocateLagrangeMultipliers()
{
    // Number of constraints
    const label m(nConstraints_);
    // Allocate the extra variables ensuring the feasibility
    if (includeExtraVars_)
    {
        extraVars_.reset(autoPtr<scalarField>::New(m, 1));
        const scalar t = mesh_.time().timeOutputValue();
        z_.reset(autoPtr<scalarField>::New(m, max(1, 0.5*c_->value(t))));

        deltaExtraVars_.reset(autoPtr<scalarField>::New(m, Zero));
        deltaZ_.reset(autoPtr<scalarField>::New(m, Zero));
    }

    doAllocateLagrangeMultipliers_ = false;
}


void Foam::ISQP::updateYS()
{
    // Compute the Lagranian and old Lagrangian derivatives
    scalarField LagrangianDerivativesOld(derivativesOld_);
    forAll(constraintDerivatives_, cI)
    {
        LagrangianDerivatives_ += lamdas_[cI]*constraintDerivatives_[cI];
        LagrangianDerivativesOld += lamdas_[cI]*constraintDerivativesOld_[cI];
    }

    if (includeBoundConstraints_)
    {
        forAll(activeDesignVars_, aI)
        {
            const label varI(activeDesignVars_[aI]);
            const scalar contr(uTilda_()[aI] - lTilda_()[aI]);
            LagrangianDerivatives_[varI] += contr;
            LagrangianDerivativesOld[varI] += contr;
        }
    }

    // Update vectors associated to the inverse Hessian matrix
    updateVectors(LagrangianDerivatives_, LagrangianDerivativesOld);
}


void Foam::ISQP::initialize()
{
    const scalarField x(designVars_().getVars(), activeDesignVars_);

    // Quantities related to design variables
    p_ = Zero;
    if (includeBoundConstraints_)
    {
        lTilda_() = scalar(1);
        uTilda_() = scalar(1);
        ls_() = scalar(1);
        us_() = scalar(1);
    }

    // Quantities related to constraints
    lamdas_ = scalar(1);
    gs_ = scalar(1);

    if (includeExtraVars_)
    {
        extraVars_() = scalar(1);
        const scalar c = c_->value(mesh_.time().timeOutputValue());
        z_() = max(1, 0.5*c);
        Info<< "Penalty multiplier (c):: " << c << endl;
    }

    // Reset eps
    eps_ = 1;
}


void Foam::ISQP::zeroUpdates()
{
    deltaP_ = Zero;
    deltaLamda_ = Zero;
    deltaGs_ = Zero;

    if (includeBoundConstraints_)
    {
        deltaLTilda_() = Zero;
        deltaLs_() = Zero;
        deltaUTilda_() = Zero;
        deltaUs_() = Zero;
    }

    if (includeExtraVars_)
    {
        deltaExtraVars_() = Zero;
        deltaZ_() = Zero;
    }
}


void Foam::ISQP::solveDeltaPEqn()
{
    addProfiling(ISQP, "ISQP::solveDeltaPEqn");
    // Explicit part of the right hand side of the deltaX equation
    scalarField FDx(-resFL());
    if (includeBoundConstraints_)
    {
        FDx +=
            (uTilda_()*resFus() + resFuTilda())/us_()
          - (lTilda_()*resFls() + resFlTilda())/ls_();
    }
    scalarField AMult(resFlamda()/lamdas_ - resFGs());
    scalarField mult(gs_/lamdas_);
    if (includeExtraVars_)
    {
        mult += extraVars_()/z_();
        AMult -= (extraVars_()*resFExtraVars() + resFz())/z_();
    }
    AMult /= mult;
    forAll(FDx, aI)
    {
        const label varI(activeDesignVars_[aI]);
        forAll(constraintDerivatives_, cI)
        {
            FDx[aI] += constraintDerivatives_[cI][varI]*AMult[cI];
        }
    }
    CGforDeltaP(FDx);
}


Foam::tmp<Foam::scalarField> Foam::ISQP::computeRHSForDeltaX
(
    const scalarField& FDx
)
{
    tmp<scalarField> trhs(tmp<scalarField>::New(-FDx));
    scalarField& rhs = trhs.ref();

    // Compute (Gs)^(-1)*Λ*A*Dp
    scalarField GsLADp(cValues_.size(), Zero);
    forAll(constraintDerivatives_, cI)
    {
        const scalarField& cDerivsI = constraintDerivatives_[cI];
        GsLADp[cI] +=
            globalSum(scalarField(cDerivsI, activeDesignVars_)*deltaP_);
    }
    GsLADp *= lamdas_/gs_;

    // Multiply with A^T
    forAll(rhs, aI)
    {
        const label varI(activeDesignVars_[aI]);
        forAll(constraintDerivatives_, cI)
        {
            rhs[aI] += constraintDerivatives_[cI][varI]*GsLADp[cI];
        }
    }

    // Contributions from bounds
    if (includeBoundConstraints_)
    {
        rhs += (lTilda_()/ls_() + uTilda_()/us_())*deltaP_;
    }

    rhs = -invHessianVectorProduct(rhs);

    rhs = 0.95*deltaP_ + 0.05*rhs;
    return trhs;
}


void Foam::ISQP::CGforDeltaP(const scalarField& FDx)
{
    addProfiling(ISQP, "ISQP::CGforDeltaP");
    autoPtr<scalarField> precon(nullptr);
    scalarField r(FDx - matrixVectorProduct(deltaP_));
    scalarField z(preconVectorProduct(r, precon));
    scalarField p(z);
    scalar res(sqrt(globalSum(r*r)));
    scalar resInit(res);
    scalar rz(globalSum(r*z));
    scalar rzOld(rz);
    label iter(0);
    do
    {
        scalarField Ap(matrixVectorProduct(p));
        scalar a = rz/globalSum(p*Ap);
        deltaP_ += a*p;
        r -= a*Ap;
        res = sqrt(globalSum(r*r));
        z = preconVectorProduct(r, precon);
        rz = globalSum(r*z);
        scalar beta = rz/rzOld;
        p = z + beta*p;
        rzOld = rz;
    } while (iter++ < maxDxIters_ && res > relTol_*resInit);
    DebugInfo
        << "CG, Solving for deltaP, Initial Residual " << resInit
        << ", Final Residual " << res
        << ", No Iterations " << iter << endl;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::matrixVectorProduct
(
    const scalarField& vector
)
{
    addProfiling(ISQP, "ISQP::MatrixVectorProduct");
    tmp<scalarField> tAp(HessianVectorProduct(vector));
    scalarField& Ap = tAp.ref();
    scalarField GsLAv(cValues_.size(), Zero);
    forAll(constraintDerivatives_, cI)
    {
        const scalarField& cDerivsI = constraintDerivatives_[cI];
        GsLAv[cI] =
            globalSum(scalarField(cDerivsI, activeDesignVars_)*vector);
    }
    scalarField mult(gs_/lamdas_);
    if (includeExtraVars_)
    {
        mult += extraVars_()/z_();
    }
    GsLAv /= mult;

    // Multiply with A^T
    forAll(Ap, aI)
    {
        const label varI(activeDesignVars_[aI]);
        forAll(constraintDerivatives_, cI)
        {
            Ap[aI] += constraintDerivatives_[cI][varI]*GsLAv[cI];
        }
    }

    // Contributions from bounds
    if (includeBoundConstraints_)
    {
        Ap += (lTilda_()/ls_() + uTilda_()/us_())*vector;
    }


    return tAp;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::preconVectorProduct
(
    const scalarField& vector,
    autoPtr<scalarField>& precon
)
{
    addProfiling(ISQP, "ISQP::preconVectorProduct");
    if (preconType_ == preconditioners::diag)
    {
        if (!precon)
        {
            precon.reset(diagPreconditioner().ptr());
        }
        return precon()*vector;
    }
    else if (preconType_ == preconditioners::invHessian)
    {
        return invHessianVectorProduct(vector);
    }
    else if (preconType_ == preconditioners::ShermanMorrison)
    {
        return ShermanMorrisonPrecon(vector);
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::diagPreconditioner()
{
    addProfiling(ISQP, "ISQP::deltaPDiagPreconditioner");
    tmp<scalarField> tpreconditioner(HessianDiag());
    scalarField& preconditioner = tpreconditioner.ref();

    // Part related to the constraints
    forAll(constraintDerivatives_, cI)
    {
        scalarField cDerivs(constraintDerivatives_[cI], activeDesignVars_);
        scalar mult(gs_[cI]/lamdas_[cI]);
        if (includeExtraVars_)
        {
            mult += extraVars_()[cI]/z_()[cI];
        }
        preconditioner += sqr(cDerivs)/mult;
    }

    if (includeBoundConstraints_)
    {
        preconditioner += lTilda_()/ls_() + uTilda_()/us_();
    }

    preconditioner = 1./preconditioner;

    return tpreconditioner;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::ShermanMorrisonPrecon
(
    const scalarField& vector
)
{
    // Recursively update the inv(LHS)*vector since the LHS consists of the
    // L-BFGS-based Hessian, computed with rank-2 updates, and the part related
    // to flow constraints, computed as rank-1 updates. In the inversion
    // process, the diagonal matrix related to bound constraints is treated as
    // the initial matrix of the L-BFGS update.

    // Contribution from bound constraints, treated as the seed of the
    // L-BFGS inverse
    refPtr<scalarField> tdiag(nullptr);
    if (includeBoundConstraints_)
    {
        tdiag.reset((lTilda_()/ls_() + uTilda_()/us_()).ptr());
    }

    // List of vectors to be used in the rank-1 updates related to the flow
    // constraitns
    PtrList<scalarField> r1Updates(cValues_.size());

    forAll(constraintDerivatives_, cI)
    {
        const scalarField& cDerivsI = constraintDerivatives_[cI];
        r1Updates.set(cI, new scalarField(cDerivsI, activeDesignVars_));
    }
    // Multipliers of the rank-1 updates
    scalarField mult(gs_/lamdas_);
    if (includeExtraVars_)
    {
        mult += extraVars_()/z_();
    }

    return
        ShermanMorrisonRank1Update(r1Updates, vector, tdiag, mult, mult.size());
}


Foam::tmp<Foam::scalarField> Foam::ISQP::ShermanMorrisonRank1Update
(
    const PtrList<scalarField>& r1Updates,
    const scalarField& p,
    const refPtr<scalarField>& diag,
    const scalarField& mult,
    label n
)
{
    auto tAp(tmp<scalarField>::New(activeDesignVars_.size(), Zero));
    scalarField& Ap = tAp.ref();

    if (n == 0)
    {
        Ap = invHessianVectorProduct(p, counter_, diag.shallowClone());
        return tAp;
    }

    do
    {
        --n;
        Ap = ShermanMorrisonRank1Update(r1Updates, p, diag, mult, n);
        tmp<scalarField> tAv =
            ShermanMorrisonRank1Update(r1Updates, r1Updates[n], diag, mult, n);
        scalarField& Av = tAv.ref();
        scalar yHs = globalSum(r1Updates[n]*Av)/mult[n];
        Ap -= Av*globalSum(r1Updates[n]*Ap)/(1 + yHs)/mult[n];
    } while (n > 0);
    return tAp;
}


void Foam::ISQP::computeNewtonDirection()
{
    // Zero the updates computed in the previous optimisation cycle
    //zeroUpdates();

    // Solve equation for deltaP_. The expensive part of the step. Everything
    // else can be computed based on this
    addProfiling(ISQP, "ISQP::computeNewtonDirection");
    solveDeltaPEqn();

    // deltaLamda
    forAll(constraintDerivatives_, cI)
    {
        const scalarField& cDerivsI = constraintDerivatives_[cI];
        deltaLamda_[cI]  =
            globalSum(scalarField(cDerivsI, activeDesignVars_)*deltaP_);
    }
    scalarField mult(gs_/lamdas_);
    if (includeExtraVars_)
    {
        mult += extraVars_()/z_();
        deltaLamda_ += (resFz() + extraVars_()*resFExtraVars())/z_();
    }
    deltaLamda_ += resFGs() - resFlamda()/lamdas_;
    deltaLamda_ /= mult;

    // deltaGs
    deltaGs_ = -(gs_*deltaLamda_ + resFlamda())/lamdas_;

    if (includeBoundConstraints_)
    {
        // deltaLs
        deltaLs_() = deltaP_ + resFls();

        // deltaUs
        deltaUs_() = -deltaP_ + resFus();

        // deltaLTilda
        deltaLTilda_() = -(lTilda_()*deltaLs_() + resFlTilda())/ls_();

        // deltaUTilda
        deltaUTilda_() = -(uTilda_()*deltaUs_() + resFuTilda())/us_();
    }

    if (includeExtraVars_)
    {
        deltaZ_() = -deltaLamda_ + resFExtraVars();
        deltaExtraVars_() = - (extraVars_()*deltaZ_() + resFz())/z_();
    }
}


Foam::scalar Foam::ISQP::lineSearch()
{
    const label n(p_.size());
    const label m(cValues_.size());
    scalar step(1.);

    if (includeBoundConstraints_)
    {
        for (label i = 0; i < n; ++i)
        {
            adjustStep(step, ls_()[i], deltaLs_()[i]);
            adjustStep(step, us_()[i], deltaUs_()[i]);
            adjustStep(step, lTilda_()[i], deltaLTilda_()[i]);
            adjustStep(step, uTilda_()[i], deltaUTilda_()[i]);
        }
    }

    // Perform bound checks and adjust step accordingly
    for (label i = 0; i < m; ++i)
    {
        adjustStep(step, lamdas_[i], deltaLamda_[i]);
        adjustStep(step, gs_[i], deltaGs_[i]);
        if (includeExtraVars_)
        {
            adjustStep(step, extraVars_()[i], deltaExtraVars_()[i]);
            adjustStep(step, z_()[i], deltaZ_()[i]);
        }
    }

    // Each processor might have computed a different step, if design variables
    // are distributed. Get the global minimum
    if (globalSum_)
    {
        reduce(step, minOp<scalar>());
    }

    step = min(1, step);

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
            if (i != maxLineSearchIters_ - 1)
            {
                updateSolution(-step);
                step *= 0.5;
            }
            else
            {
                Info<< tab << "Line search did not converge. "
                    << "Procceding with the last compute step" << endl;
            }
        }
    }

    if (debug > 1)
    {
        Info<< "Step after line search is " << step << nl <<  endl;
    }

    return maxRes;
}


void Foam::ISQP::adjustStep
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


void Foam::ISQP::updateSolution(const scalar step)
{
    p_ += step*deltaP_;
    lamdas_ += step*deltaLamda_;
    gs_ += step*deltaGs_;
    if (includeBoundConstraints_)
    {
        lTilda_() += step*deltaLTilda_();
        ls_() += step*deltaLs_();
        uTilda_() += step*deltaUTilda_();
        us_() += step*deltaUs_();
    }
    if (includeExtraVars_)
    {
        extraVars_() += step*deltaExtraVars_();
        z_() += step*deltaZ_();
    }
}


Foam::tmp<Foam::scalarField> Foam::ISQP::computeResiduals()
{
    const label n(activeDesignVars_.size());
    const label m(cValues_.size());
    label size(includeBoundConstraints_ ? 5*n + 2*m : n + 2*m);
    if (includeExtraVars_)
    {
        size += 2*m;
    }
    tmp<scalarField> tres(tmp<scalarField>::New(size, Zero));
    scalarField& res = tres.ref();

    label iRes(0);

    // Gradient of the Lagrangian
    res.rmap(resFL()(), identity(n));
    iRes = n;

    // Inequality constraints slacks
    res.rmap(resFGs()(), identity(m, iRes));
    iRes += m;

    // Inequality constraints complementarity slackness
    res.rmap(resFlamda()(), identity(m, iRes));
    iRes += m;

    if (includeBoundConstraints_)
    {
        // Lower bounds slacks
        res.rmap(resFls()(), identity(n, iRes));
        iRes += n;

        // Upper bounds slacks
        res.rmap(resFus()(), identity(n, iRes));
        iRes += n;

        // Lower bounds complementarity slackness
        res.rmap(resFlTilda()(), identity(n, iRes));
        iRes += n;

        // Upper bounds complementarity slackness
        res.rmap(resFuTilda()(), identity(n, iRes));
        iRes += n;
    }

    if (includeExtraVars_)
    {
        // Lagragian derivative wrt the extra variables
        res.rmap(resFExtraVars()(), identity(m, iRes));
        iRes += m;

        // Lagrange multipliers for the extra variables positiveness
        res.rmap(resFz(), identity(m, iRes));
        iRes += m;
    }

    return tres;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFL()
{
    tmp<scalarField> tgradL
        (tmp<scalarField>::New(objectiveDerivatives_, activeDesignVars_));
    scalarField& gradL = tgradL.ref();

    scalarField Hp(HessianVectorProduct(p_));
  //scalarField Hp = SR1HessianVectorProduct(p_);
    gradL += Hp;

    if (debug > 2)
    {
        scalarField H1Hp(invHessianVectorProduct(Hp));
        Info << "Diff H1Hp - p " << gSum(mag(H1Hp - p_)) << endl;
    }

    forAll(constraintDerivatives_, cI)
    {
        gradL +=
            lamdas_[cI]
           *scalarField(constraintDerivatives_[cI], activeDesignVars_);
    }

    if (includeBoundConstraints_)
    {
        gradL += uTilda_() - lTilda_();
    }

    return tgradL;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::invHFL()
{
    tmp<scalarField> tinvHFL
        (tmp<scalarField>::New(objectiveDerivatives_, activeDesignVars_));
    scalarField& invHFL = tinvHFL.ref();

    forAll(constraintDerivatives_, cI)
    {
        invHFL +=
            lamdas_[cI]
           *scalarField(constraintDerivatives_[cI], activeDesignVars_);
    }

    if (includeBoundConstraints_)
    {
        invHFL += uTilda_() - lTilda_();
    }

    invHFL = invHessianVectorProduct(invHFL);
    invHFL += p_;

    return tinvHFL;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFGs()
{
    tmp<scalarField> tFGs
    (
        tmp<scalarField>::New(gs_ + cValues_ - max((1 - cRed_)*cValues_, Zero))
    );
    scalarField& FGs = tFGs.ref();

    forAll(constraintDerivatives_, cI)
    {
        FGs[cI] +=
            globalSum
            (
                scalarField(constraintDerivatives_[cI], activeDesignVars_)*p_
            );
    }

    if (includeExtraVars_)
    {
        FGs -= extraVars_();
    }

    return tFGs;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFlamda()
{
    return (lamdas_*gs_ - eps_);
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFls()
{
    if (includeBoundConstraints_)
    {
        const scalarField x(designVars_().getVars(), activeDesignVars_);
        const scalarField xMin
            (designVars_().lowerBounds()(), activeDesignVars_);

        return (x + p_ - xMin - ls_());
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFus()
{
    if (includeBoundConstraints_)
    {
        const scalarField x(designVars_().getVars(), activeDesignVars_);
        const scalarField xMax
            (designVars_().upperBounds()(), activeDesignVars_);

        return (xMax - x - p_ - us_());
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFlTilda()
{
    if (includeBoundConstraints_)
    {
        return (lTilda_()*ls_() - eps_);
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFuTilda()
{
    if (includeBoundConstraints_)
    {
        return (uTilda_()*us_() - eps_);
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFExtraVars()
{
    if (includeExtraVars_)
    {
        return (c_->value(mesh_.time().timeOutputValue()) - lamdas_ - z_());
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::ISQP::resFz()
{
    if (includeExtraVars_)
    {
        return (z_()*extraVars_() - eps_);
    }
    return nullptr;
}


void Foam::ISQP::solveSubproblem()
{
    addProfiling(ISQP, "ISQP::solveSubproblem");
    zeroUpdates();
    if (includeBoundConstraints_ || !cValues_.empty())
    {
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
        } while
        (
            iter++ < maxNewtonIters_ && (eps_ > epsMin_ || resMax > 0.9*eps_)
        );
        Info<< "Finished solving the QP subproblem" << nl << endl;
        if (iter == maxNewtonIters_)
        {
            WarningInFunction
                << "Iterative solution of the QP problem did not converge"
                << endl;
        }
        if (debug > 2)
        {
            scalarField vars(designVars_().getVars(), activeDesignVars_);
            scalarField newVars(vars + p_);
            Info<< "Min of updated vars " << gMin(newVars) << endl;
            Info<< "Max of updated vars " << gMax(newVars) << endl;

            Info<< "Min of lamda " << gMin(lamdas_) << endl;
            Info<< "Max of gs " << gMax(gs_) << endl;
            Info<< "Max of lamda*gs " << gMax(lamdas_*gs_) << endl;

            if (includeBoundConstraints_)
            {
                Info<< "Min of lTilda " << gMin(lTilda_()) << endl;
                Info<< "Min of uTilda " << gMin(uTilda_()) << endl;
                Info<< "Min of ls " << gMin(ls_()) << endl;
                Info<< "Min of us " << gMin(us_()) << endl;
                Info<< "Max of lTilda*ls " << gMax(lTilda_()*ls_()) << endl;
                Info<< "Max of uTilda*us " << gMax(uTilda_()*us_()) << endl;
            }
            if (includeExtraVars_)
            {
                Info<< "Min of extraVars " << gMin(extraVars_()) << endl;
                Info<< "Max of extraVars*z " << gMax(extraVars_()*z_()) << endl;
            }
        }
        if (includeExtraVars_)
        {
            DebugInfo
                << "Constraint penalization variables (y) " << extraVars_()
                << endl;
        }
    }
    else
    {
        computeNewtonDirection();
        lineSearch();
    }

    // Pass update to correction field
    correction_.rmap(p_, activeDesignVars_);
    if (!counter_)
    {
        correction_ *= eta_;
    }
    else
    {
        correction_ *= etaHessian_;
    }
}


void Foam::ISQP::storeOldFields()
{
    derivativesOld_ = objectiveDerivatives_;
    if (constraintDerivativesOld_.empty())
    {
        constraintDerivativesOld_.setSize(constraintDerivatives_.size());
    }
    forAll(constraintDerivativesOld_, cI)
    {
        constraintDerivativesOld_[cI] = constraintDerivatives_[cI];
    }
    correctionOld_ = correction_;
}


Foam::scalar Foam::ISQP::meritFunctionConstraintPart() const
{
    // Assumes that all constraints are known by all processors
    // What about constraints directly imposed on distributed design variables?
    // These should be met in each iteration of the algorithm, so,
    // most probably, there is no problem
    return sum(pos(cValues_)*cValues_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ISQP::ISQP
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    LBFGS(mesh, dict, designVars, nConstraints, type),
    SQPBase(mesh, dict, designVars, *this, type),

    doAllocateLagrangeMultipliers_(true),
    includeBoundConstraints_
    (
        designVars->upperBounds() && designVars->lowerBounds()
    ),
    includeExtraVars_
    (
        coeffsDict(type).getOrDefault<bool>("includeExtraVars", false)
    ),
    p_(activeDesignVars_.size(), Zero),
    gs_(nConstraints_, 1),
    lTilda_
    (
        includeBoundConstraints_ && found("lTilda") ?
        new scalarField("lTilda", *this, activeDesignVars_.size()) :
        nullptr
    ),
    ls_(nullptr),
    uTilda_
    (
        includeBoundConstraints_ && found("uTilda") ?
        new scalarField("uTilda", *this, activeDesignVars_.size()) :
        nullptr
    ),
    us_(nullptr),
    extraVars_(nullptr),
    z_(nullptr),
    c_(nullptr),
    deltaP_(activeDesignVars_.size(), Zero),
    deltaLamda_(nConstraints_, Zero),
    deltaGs_(nConstraints_, Zero),
    deltaLTilda_(nullptr),
    deltaLs_(nullptr),
    deltaUTilda_(nullptr),
    deltaUs_(nullptr),
    deltaExtraVars_(nullptr),
    deltaZ_(nullptr),
    eps_(1),
    epsMin_(coeffsDict(type).getOrDefault<scalar>("epsMin", 1.e-07)),
    maxNewtonIters_(coeffsDict(type).getOrDefault<label>("maxIters", 1000)),
    maxLineSearchIters_
    (
        coeffsDict(type).getOrDefault<label>("maxLineSearchIters", 10)
    ),
    maxDxIters_(coeffsDict(type).getOrDefault<label>("maxDpIters", 1000)),
    relTol_(coeffsDict(type).getOrDefault<scalar>("relTol", 0.01)),
    preconType_
    (
        preconditionerNames.getOrDefault
        (
            "preconditioner", coeffsDict(type), preconditioners::diag
        )
    ),
    cRed_
        (coeffsDict(type).getOrDefault<scalar>("targetConstraintReduction", 1)),
    disableDamping_
        (coeffsDict(type).getOrDefault<bool>("disableDamping", false)),
    meritFunctionFile_(nullptr)
{
    Info<< "Preconditioner type of the SQP subproblem is ::"
        << preconditionerNames.names()[preconType_]
        << endl;
    if (!disableDamping_)
    {
        // Always apply damping of y in ISQP
        useYDamping_ = true;
        useSDamping_ = false;
    }

    // Determine c if necessary
    if (includeExtraVars_)
    {
        if (coeffsDict(type).found("c"))
        {
            Info<< "Reading constraint penalty function type from dict" << endl;
            c_.reset(Function1<scalar>::New("c", coeffsDict(type)));
        }
        else
        {
            Info<< "Setting constant penalty factor" << endl;
            c_.reset(new Function1Types::Constant<scalar>("c", 100));
        }
    }

    // Allocate multipliers and slack variables for the bound constraints
    allocateBoundMultipliers();
    allocateLagrangeMultipliers();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ISQP::computeCorrection()
{
    // Update sizes of fields related to the active design variables
    updateSizes();

    // The first iteration uses a unitary Hessian. No need to update
    LagrangianDerivatives_ = objectiveDerivatives_;
    if (counter_)
    {
        updateYS();
    }

    // Initiaze variables
    initialize();

    // Solve subproblem using a Newton optimiser
    solveSubproblem();

    // Store fields for the next iteration and write them to file
    storeOldFields();

    // Increase counter
    ++counter_;
}


Foam::scalar Foam::ISQP::computeMeritFunction()
{
    mu_ = max(pos(cValues_)*lamdas_) + delta_;
    scalar L = objectiveValue_ + mu_*sum(pos(cValues_)*cValues_);

    return L;
}


Foam::scalar Foam::ISQP::meritFunctionDirectionalDerivative()
{
    scalar deriv =
        globalSum(objectiveDerivatives_*correction_)
      - mu_*sum(pos(cValues_)*cValues_);

    return deriv;
}


void Foam::ISQP::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


bool Foam::ISQP::writeData(Ostream& os) const
{
    if (includeBoundConstraints_)
    {
        uTilda_().writeEntry("uTilda", os);
        lTilda_().writeEntry("lTilda", os);
    }

    return LBFGS::writeData(os) && SQPBase::addToFile(os);
}


bool Foam::ISQP::writeAuxiliaryData()
{
    return SQPBase::writeMeritFunction(*this);
}


// ************************************************************************* //
