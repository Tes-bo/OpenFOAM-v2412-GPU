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

\*---------------------------------------------------------------------------*/

#include "nullSpace.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nullSpace, 1);
    addToRunTimeSelectionTable
    (
        updateMethod,
        nullSpace,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        nullSpace,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nullSpace::initialise()
{
    eps_ = 1;
    updateViolatedConstraintsSubsets();

    const label m = iTildaEps_[0].size();
    mu_ = scalarField(m, 1);
    dualMu_ = scalarField(m, 1);
    deltaMu_ = scalarField(m, Zero);
    deltaDualMu_ = scalarField(m, Zero);

    const label nl = iTildaEps_[1].size();
    l_ = scalarField(nl, 1);
    dualL_ = scalarField(nl, 1);
    deltaL_ = scalarField(nl, Zero);
    deltaDualL_ = scalarField(nl, Zero);

    const label nu = iTildaEps_[2].size();
    u_ = scalarField(nu, 1);
    dualU_ = scalarField(nu, 1);
    deltaU_ = scalarField(nu, Zero);
    deltaDualU_ = scalarField(nu, Zero);
}


void Foam::nullSpace::zeroUpdates()
{
    deltaMu_ = Zero;
    deltaDualMu_ = Zero;
    deltaL_ = Zero;
    deltaDualL_ = Zero;
    deltaU_ = Zero;
    deltaDualU_ = Zero;
}


void Foam::nullSpace::solveDualProblem()
{
    label nConstraints =
        globalSum
        (
            iTildaEps_[0].size()
          + iTildaEps_[1].size()
          + iTildaEps_[2].size()
        );
    if (nConstraints && solveDualProblem_)
    {
        scalar resMax(gMax(mag(computeResiduals())));
        label iter(0);
        do
        {
            DebugInfo
                << "Dual problem Newton iter " << iter << nl << endl;

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
        } while
        (
            iter++ < maxNewtonIters_
         && (eps_ > dualTolerance_ || resMax > 0.9*eps_)
        );

        Info<< "Solved the dual Newton problem in " << iter << " iterations "
            << nl << endl;
        Info<< "fluid related Lagrange mults " << mu_ << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::nullSpace::computeResiduals()
{
    const labelList& iFlow = iTildaEps_[0];
    const labelList& iLower = iTildaEps_[1];
    const labelList& iUpper = iTildaEps_[2];
    const label m = iFlow.size();
    const label nl = iLower.size();
    const label nu = iUpper.size();
    tmp<scalarField> tres(tmp<scalarField>::New(2*(m + nl + nu), Zero));
    scalarField& res = tres.ref();

    label iRes = 0;

    // dLdMu
    scalarField dLdx(objectiveDerivatives_, activeDesignVars_);
    forAll(iFlow, i)
    {
        dLdx +=
            mu_[i]
           *scalarField(constraintDerivatives_[iFlow[i]], activeDesignVars_);
    }
    forAll(iLower, i)
    {
        dLdx[iLower[i]] -=  l_[i];
    }
    forAll(iUpper, i)
    {
        dLdx[iUpper[i]] += u_[i];
    }

    forAll(iFlow, i)
    {
        label ic = iFlow[i];
        scalar gradPart
        (
            2*globalSum
            (
                dLdx*scalarField(constraintDerivatives_[ic], activeDesignVars_)
            )
        );
        res[iRes++] = gradPart - dualMu_[i];
    }

    // mu complementarity slackness
    forAll(iFlow, i)
    {
        res[iRes++] = mu_[i]*dualMu_[i] - eps_;
    }

    // dKdl
    forAll(iLower, i)
    {
        res[iRes++] = - dualL_[i] - 2*dLdx[iLower[i]];
    }

    // dKdu
    forAll(iUpper, i)
    {
        res[iRes++] = - dualU_[i] + 2*dLdx[iUpper[i]];
    }

    // l complementarity slackness
    forAll(iLower, i)
    {
        res[iRes++] = l_[i]*dualL_[i] - eps_;
    }

    // u complementarity slackness
    forAll(iUpper, i)
    {
        res[iRes++] = u_[i]*dualU_[i] - eps_;
    }

    return tres;
}


void Foam::nullSpace::computeNewtonDirection()
{
    // Zero the updates computed in the previous optimisation cycle
    zeroUpdates();

    const scalarField res(computeResiduals());

    const labelList& iFlow = iTildaEps_[0];
    const labelList& iLower = iTildaEps_[1];
    const labelList& iUpper = iTildaEps_[2];
    const label m = iFlow.size();
    const label nl = iLower.size();
    const label nu = iUpper.size();

    // Update of the flow-related primal and dual Lagrange multipliers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarField diagKsi(nl, Zero);
    scalarField rhsKsi(nl, Zero);
    scalarField diagEta(nu, Zero);
    scalarField rhsEta(nu, Zero);
    if (nl)
    {
        diagKsi = scalar(1)/(2 + dualL_/l_);
        rhsKsi =
           -(
                SubField<scalar>(res, nl, 2*m)
              + SubField<scalar>(res, nl, 2*m + nl + nu)/l_
            )*diagKsi;
    }
    if (nu)
    {
        diagEta = scalar(1)/(2 + dualU_/u_);
        rhsEta =
           -(
                SubField<scalar>(res, nu, 2*m + nl)
              + SubField<scalar>(res, nu, 2*m + 2*nl + nu)/u_
            )*diagEta;
    }

    scalarSquareMatrix lhs(m, m, Zero);
    scalarField rhs(m, Zero);
    for (label i = 0; i < m; ++i)
    {
        const label ic = iFlow[i];
        scalarField dci(constraintDerivatives_[ic], activeDesignVars_);
        lhs[i][i] += 2*globalSum(dci*dci) + dualMu_[i]/mu_[i];
        rhs[i] -= res[i] + res[m + i]/mu_[i];

        scalarField dciKsi(scalarField(dci, iLower));
        scalarField dciEta(scalarField(dci, iUpper));
        lhs[i][i] -=
            4*globalSum(dciKsi*dciKsi*diagKsi)
         +  4*globalSum(dciEta*dciEta*diagEta);
        rhs[i] += 2*globalSum(dciKsi*rhsKsi) - 2*globalSum(dciEta*rhsEta);
        for (label j = i + 1; j < m; ++j)
        {
            const label jc = iFlow[j];
            scalarField dcj(constraintDerivatives_[jc], activeDesignVars_);
            scalar ij = 2*globalSum(dci*dcj);
            lhs[i][j] += ij;
            lhs[j][i] += ij;

            scalarField dciKsi(scalarField(dci, iLower));
            scalarField dcjKsi(scalarField(dcj, iLower));
            scalarField dciEta(scalarField(dci, iUpper));
            scalarField dcjEta(scalarField(dcj, iUpper));
            scalar ksiContr = 4*globalSum(dciKsi*dcjKsi*diagKsi);
            scalar etaContr = 4*globalSum(dciEta*dcjEta*diagEta);
            lhs[i][j] -= ksiContr + etaContr;
            lhs[j][i] -= ksiContr + etaContr;
        }
    }

    // Update flow related Lagrange multipliers
    if (m)
    {
        solve(deltaMu_, lhs, rhs);
        deltaDualMu_ = -(deltaMu_*dualMu_ + SubField<scalar>(res, m, m))/mu_;
    }

    // Compute the rest of the corrections using backwards substitution
    deltaL_ = Zero;
    deltaU_ = Zero;
    forAll(iFlow, i)
    {
        const label ic = iFlow[i];
        scalarField dci(constraintDerivatives_[ic], activeDesignVars_);
        deltaL_ += scalarField(dci, iLower)*deltaMu_[i];
        deltaU_ -= scalarField(dci, iUpper)*deltaMu_[i];
    }
    deltaL_ *= 2*diagKsi;
    deltaL_ += rhsKsi;

    deltaU_ *= 2*diagEta;
    deltaU_ += rhsEta;

    if (nl)
    {
        deltaDualL_ =
          - (dualL_*deltaL_ + SubField<scalar>(res, nl, 2*m + nl + nu))
           /l_;
    }

    if (nu)
    {
        deltaDualU_ =
          - (dualU_*deltaU_ + SubField<scalar>(res, nu, 2*m + 2*nl + nu))
           /u_;
    }
}


Foam::scalar Foam::nullSpace::lineSearch()
{
    scalar step(1.);

    // Perform bound checks and adjust step accordingly
    forAll(iTildaEps_[0], i)
    {
        adjustStep(step, mu_[i], deltaMu_[i]);
        adjustStep(step, dualMu_[i], deltaDualMu_[i]);
    }

    forAll(iTildaEps_[1], i)
    {
        adjustStep(step, l_[i], deltaL_[i]);
        adjustStep(step, dualL_[i], deltaDualL_[i]);
    }

    forAll(iTildaEps_[2], i)
    {
        adjustStep(step, u_[i], deltaU_[i]);
        adjustStep(step, dualU_[i], deltaDualU_[i]);
    }

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


void Foam::nullSpace::adjustStep
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


void Foam::nullSpace::updateSolution(const scalar step)
{
    mu_ += step*deltaMu_;
    dualMu_ += step*deltaDualMu_;

    l_ += step*deltaL_;
    dualL_ += step*deltaDualL_;

    u_ += step*deltaU_;
    dualU_ += step*deltaDualU_;
}


void Foam::nullSpace::updateViolatedConstraintsSubsets()
{
    updateViolatedIndices(0, cValues_);

    if (includeBoundConstraints_)
    {
        scalarField lowerBounds
            (designVars_->lowerBoundsRef() - designVars_(), activeDesignVars_);
        updateViolatedIndices(1, lowerBounds);

        scalarField upperBounds
            (designVars_() - designVars_->upperBoundsRef(), activeDesignVars_);
        updateViolatedIndices(2, upperBounds);
    }

    statistics(iTilda_, "violated");
    statistics(iTildaEps_, "violated-up-to-eps");
}


void Foam::nullSpace::updateNullAndRangeSpaceSubsets()
{
    if (solveDualProblem_)
    {
        updateCorrectionIndices(0, mu_, dualMu_);
        updateCorrectionIndices(1, l_, dualL_);
        updateCorrectionIndices(2, u_, dualU_);

        statistics(iHat_, "non-tangent,violated");
        statistics(iRangeSpace_, "to-be-reduced");
    }
    else
    {
        iHat_ = iTildaEps_;
        iRangeSpace_ = iTildaEps_;
    }
}


void Foam::nullSpace::updateViolatedIndices
(
    const label i,
    const scalarField& constraints
)
{
    // Set violated constraints
    labelList& subset = iTilda_[i];
    subset.setSize(constraints.size(), -1);
    label iViolated = Zero;
    forAll(constraints, i)
    {
        if (constraints[i] >= 0)
        {
            subset[iViolated++] = i;
        }
    }
    subset.setSize(iViolated);

    // Append violated constraints up to epsConstr_
    DynamicList<label> iTildaEps(subset);
    forAll(constraints, i)
    {
        if (constraints[i] >= -epsConstr_ && constraints[i] < 0)
        {
            iTildaEps.push_back(i);
        }
    }
    iTildaEps_[i].transfer(iTildaEps);
}


void Foam::nullSpace::updateCorrectionIndices
(
    const label i,
    const scalarField& LagrangeMults,
    const scalarField& dual
)
{
    // Subset with non-zero Lagrange multipliers
    const labelList& firstAddr = iTildaEps_[i];
    labelList& subset = iHat_[i];
    subset.setSize(LagrangeMults.size(), -1);
    label iViolated(Zero);

    // Loop over all indices included in iTilda
    forAll(iTilda_[i], j)
    {
        // The criterion of adding a contraint index to the subset
        // is based on whether the Lagrange multiplier is larger than
        // its dual (since their product should, theorerically be zero). This
        // avoids comparisons with arbitrary small numbers
        if (LagrangeMults[j] > dual[j])
        {
            subset[iViolated++] = firstAddr[j];
        }
    }
    // Loop over the remaining indices included in iTildaEps_ and append
    // elements with a positive Lagrange multiplier to iRangeSpace too
    DynamicList<label> iRangeSpace(iTilda_[i]);
    for (label j = iTilda_[i].size(); j < iTildaEps_[i].size(); ++j)
    {
        if (LagrangeMults[j] > dual[j])
        {
            subset[iViolated++] = firstAddr[j];
            iRangeSpace.push_back(firstAddr[j]);
        }
    }
    subset.setSize(iViolated);
    iRangeSpace_[i].transfer(iRangeSpace);
}


void Foam::nullSpace::correction()
{
    // Compute the null space part of the update
    const scalarField activeDerivs(objectiveDerivatives_, activeDesignVars_);
    scalarField rJ(Av(activeDerivs, iHat_));
    scalarField lamdaJ(constraintRelatedUpdate(rJ, iHat_));
    scalarField ksiJ(activeDerivs - ATv(lamdaJ, iHat_));

    // Compute the range space part of the update
    scalarField g(activeConstraints(iRangeSpace_));
    scalarField lamdaC(constraintRelatedUpdate(g, iRangeSpace_));
    scalarField ksiC(ATv(lamdaC, iRangeSpace_));

    if (debug)
    {
        scalar magKsiJ = sqrt(globalSum(sqr(ksiJ))) ;
        scalarField ksiJUnit(ksiJ/(magKsiJ + SMALL));
        for (const label i : iHat_[0])
        {
            scalarField cDerivs
                (constraintDerivatives_[i], activeDesignVars_);
            cDerivs /= sqrt(globalSum(sqr(cDerivs))) + SMALL;
            Info<< "\tNull space update projected to the " << i
                << "-th flow constraint gradient "
                << globalSum(ksiJUnit*cDerivs)
                << endl;
        }

        scalar sumLowConsts(Zero);
        scalar sumUppConsts(Zero);
        for (const label il : iHat_[1])
        {
            sumLowConsts += ksiJUnit[il];
        }
        for (const label iu : iHat_[2])
        {
            sumUppConsts += ksiJUnit[iu];
        }
        if (globalSum_)
        {
            reduce(sumLowConsts, sumOp<scalar>());
            reduce(sumUppConsts, sumOp<scalar>());
        }
        Info<< "\tSum of projections to the lower bound constraints "
            << sumLowConsts << endl;
        Info<< "\tSum of projections to the upper bound constraints "
            << sumUppConsts << endl;
    }

    scalar maxKsiJ = gMax(mag(ksiJ));
    if (adaptiveStep_ && maxKsiJ > VSMALL)
    {
        if (counter_ < lastAcceleratedCycle_)
        {
            aJ_ = maxDVChange_()/eta_/maxKsiJ;
        }
        else if (strictMaxDVChange_)
        {
            aJ_ = min(aJ_, maxDVChange_()/eta_/maxKsiJ);
        }

        /*
        // Can become unstable close to the optimum
        aJ_ =
            targetObjectiveReduction_*objectiveValue_
           /(eta_*mag(globalSum(activeDerivs*ksiJ)));
        */
        DebugInfo
            << "aJ set to " << aJ_ << endl;
    }

    scalar maxKsiC = gMax(mag(ksiC));
    if (adaptiveStep_ && maxKsiC > VSMALL)
    {
        aC_ =
            min
            (
                targetConstraintReduction_/eta_,
                maxDVChange_()/eta_/maxKsiC
            );
    }

    if (debug)
    {
        Info<< "Mag of ksiJ and ksiC "
            << Foam::sqrt(globalSum(ksiJ*ksiJ)) << " "
            << Foam::sqrt(globalSum(ksiC*ksiC))
            << endl;

        Info<< "Inner product of ksiJ and ksiC " << globalSum(ksiC*ksiJ)
            << endl;
        Info<< "max eta*aJ*ksiJ " << gMax(mag(eta_*aJ_*ksiJ)) << endl;
        Info<< "max eta*aC*ksiC " << gMax(mag(eta_*aC_*ksiC)) << endl;
    }

    correction_.rmap(-eta_*(aJ_*ksiJ + aC_*ksiC), activeDesignVars_);
}


Foam::tmp<Foam::scalarField> Foam::nullSpace::Av
(
    const scalarField& v,
    const labelListList& subsets
)
{
    const labelList& iFlow = subsets[0];
    const labelList& iLower = subsets[1];
    const labelList& iUpper = subsets[2];
    label m = iFlow.size();
    label nl = iLower.size();
    label nu = iUpper.size();
    if (v.size() != activeDesignVars_.size())
    {
        FatalErrorInFunction
            << "Input vector size not equal to the active design variables"
            << exit(FatalError);
    }
    auto tmvp = tmp<scalarField>::New(m + nl + nu, Zero);
    scalarField& mvp = tmvp.ref();

    // Flow-related constraints
    forAll(iFlow, ic)
    {
        scalarField di(constraintDerivatives_[iFlow[ic]], activeDesignVars_);
        mvp[ic] += globalSum(di*v);
    }

    // Lower bounds constraints
    forAll(iLower, il)
    {
        mvp[m + il] = -v[iLower[il]];
    }

    // Flow-upper bounds constraints interaction
    forAll(iUpper, iu)
    {
        mvp[m + nl + iu] = v[iUpper[iu]];
    }

    return tmvp;
}


Foam::tmp<Foam::scalarField> Foam::nullSpace::ATv
(
    const scalarField& v,
    const labelListList& subsets
)
{
    const labelList& iFlow = subsets[0];
    const labelList& iLower = subsets[1];
    const labelList& iUpper = subsets[2];
    label m = iFlow.size();
    label nl = iLower.size();
    label nu = iUpper.size();
    if (v.size() != m + nl + nu)
    {
        FatalErrorInFunction
            << "Input vector size not equal to the active constraints"
            << exit(FatalError);
    }
    auto tmvp = tmp<scalarField>::New(activeDesignVars_.size(), Zero);
    scalarField& mvp = tmvp.ref();

    // Flow-related constraints
    forAll(iFlow, ic)
    {
        scalarField di(constraintDerivatives_[iFlow[ic]], activeDesignVars_);
        mvp += di*v[ic];
    }

    // Lower bounds constraints
    forAll(iLower, il)
    {
        mvp[iLower[il]] -= v[m + il];
    }

    // Upper bounds constraints
    forAll(iUpper, iu)
    {
        mvp[iUpper[iu]] += v[m + nl + iu];
    }

    return tmvp;
}


Foam::tmp<Foam::scalarField> Foam::nullSpace::activeConstraints
(
    const labelListList& subsets
)
{
    const labelList& iFlow = subsets[0];
    const labelList& iLower = subsets[1];
    const labelList& iUpper = subsets[2];
    label m = iFlow.size();
    label nl = iLower.size();
    label nu = iUpper.size();
    auto tcons = tmp<scalarField>::New(m + nl + nu, Zero);
    scalarField& cons = tcons.ref();

    label ic = 0;
    for (const label i : iFlow)
    {
        cons[ic++] = cValues_[i];
    }

    const autoPtr<scalarField>& lowerBounds = designVars_->lowerBounds();
    for (const label i : iLower)
    {
        const label iActive = activeDesignVars_[i];
        cons[ic++] = bcMult_*(lowerBounds()[iActive] - designVars_()[iActive]);
    }

    const autoPtr<scalarField>& upperBounds = designVars_->upperBounds();
    for (const label i : iUpper)
    {
        const label iActive = activeDesignVars_[i];
        cons[ic++] = bcMult_*(designVars_()[iActive] - upperBounds()[iActive]);
    }

    return tcons;
}


Foam::tmp<Foam::scalarField> Foam::nullSpace::constraintRelatedUpdate
(
    const scalarField& rhs,
    const labelListList& subsets
)
{
    const labelList& iFlow = subsets[0];
    const labelList& iLower = subsets[1];
    const labelList& iUpper = subsets[2];
    label m = iFlow.size();
    label nl = iLower.size();
    label nu = iUpper.size();
    if (rhs.size() != m + nl + nu)
    {
        FatalErrorInFunction
            << "rhs should have the dimensions of the active constraints"
            << exit(FatalError);
    }
    auto tres = tmp<scalarField>::New(rhs.size(), Zero);
    scalarField& res = tres.ref();

    scalarSquareMatrix lhs(m, m, Zero);
    scalarField r(m, Zero);

    forAll(iFlow, i)
    {
        const label ic = iFlow[i];
        scalarField dci(constraintDerivatives_[ic], activeDesignVars_);

        // Diagonal part of lhs
        lhs[i][i] += globalSum(dci*dci) ;
        scalar lowerContr(Zero);
        scalar upperContr(Zero);
        for (const label il : iLower)
        {
            lowerContr -= dci[il]*dci[il];
        }
        for (const label iu : iUpper)
        {
            upperContr -= dci[iu]*dci[iu];
        }
        if (globalSum_)
        {
            reduce(lowerContr, sumOp<scalar>());
            reduce(upperContr, sumOp<scalar>());
        }
        lhs[i][i] += lowerContr + upperContr;

        // Off diagonal part of lhs
        for (label j = i + 1; j < m; ++j)
        {
            const label jc = iFlow[j];
            scalarField dcj(constraintDerivatives_[jc], activeDesignVars_);
            scalar ij = globalSum(dci*dcj);
            lhs[i][j] += ij;
            lhs[j][i] += ij;

            lowerContr = Zero;
            for (const label il : iLower)
            {
                lowerContr -= dci[il]*dcj[il];
            }
            upperContr = Zero;
            for (const label iu : iUpper)
            {
                upperContr -= dci[iu]*dcj[iu];
            }

            if (globalSum_)
            {
                reduce(lowerContr, sumOp<scalar>());
                reduce(upperContr, sumOp<scalar>());
            }
            lhs[i][j] += lowerContr + upperContr;
            lhs[j][i] += lowerContr + upperContr;
        }

        // rhs
        r[i] = rhs[i];
        lowerContr = Zero;
        forAll(iLower, il)
        {
            lowerContr += dci[iLower[il]]*rhs[m + il];
        }
        upperContr = Zero;
        forAll(iUpper, iu)
        {
            upperContr -= dci[iUpper[iu]]*rhs[m + nl + iu];
        }

        if (globalSum_)
        {
            reduce(lowerContr, sumOp<scalar>());
            reduce(upperContr, sumOp<scalar>());
        }
        r[i] += lowerContr + upperContr;
    }

    // Compute solution
    SubField<scalar>(res, nl, m) = SubField<scalar>(rhs, nl, m);
    SubField<scalar>(res, nu, nl + m) = SubField<scalar>(rhs, nu, nl + m);
    if (m)
    {
        scalarField solF(m, Zero);
        solve(solF, lhs, r);
        SubField<scalar>(res, m, 0) = solF;

        forAll(iFlow, i)
        {
            scalarField dci(constraintDerivatives_[iFlow[i]], activeDesignVars_);
            forAll(iLower, il)
            {
                res[m + il] += dci[iLower[il]]*solF[i];
            }
            forAll(iUpper, iu)
            {
                res[m + nl + iu] -= dci[iUpper[iu]]*solF[i];
            }
        }
    }
    return tres;
}


void Foam::nullSpace::statistics
(
    const labelListList& subset,
    const word& description
)
{
    DebugInfo
        << "Number of flow constraints (" << description << ") "
        << subset[0].size() << nl;
    if (includeBoundConstraints_)
    {
        DebugInfo
            << "Number of lower bounds constraints (" << description << ") "
            << globalSum(subset[1].size())
            << nl;
        DebugInfo
            << "Number of upper bounds constraints (" << description << ") "
            << globalSum(subset[2].size())
            << nl;
    }
    DebugInfo<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nullSpace::nullSpace
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
    includeBoundConstraints_
    (
        designVars->upperBounds() && designVars->lowerBounds()
    ),
    mu_(nConstraints_, Zero),
    l_(),
    u_(),
    solveDualProblem_
        (coeffsDict(type).getOrDefault<bool>("solveDualProblem", true)),
    dualMu_(nConstraints_, Zero),
    dualL_(),
    dualU_(),
    iTilda_(3),
    iTildaEps_(3),
    iHat_(3),
    iRangeSpace_(3),
    deltaMu_(nConstraints_, Zero),
    deltaL_(),
    deltaU_(),
    deltaDualMu_(nConstraints_, Zero),
    deltaDualL_(),
    deltaDualU_(),
    eps_(1),
    maxNewtonIters_(coeffsDict(type).getOrDefault<label>("maxIters", 6000)),
    maxLineSearchIters_
    (
        coeffsDict(type).getOrDefault<label>("maxLineSearchIters", 10)
    ),
    maxCGIters_
        (coeffsDict(type).getOrDefault<label>("maxCGIters", maxNewtonIters_)),
    dualTolerance_
        (coeffsDict(type).getOrDefault<scalar>("dualTolerance", 1.e-12)),
    epsConstr_
    (
        coeffsDict(type).
            getOrDefault<scalar>("violatedConstraintsThreshold", 1.e-3)
    ),
    epsPerturb_(coeffsDict(type).  getOrDefault<scalar>("perturbation", 1.e-2)),
    aJ_(coeffsDict(type).getOrDefault<scalar>("aJ", 1)),
    aC_(coeffsDict(type).getOrDefault<scalar>("aC", 1)),
    adaptiveStep_(coeffsDict(type).getOrDefault<bool>("adaptiveStep", true)),
    lastAcceleratedCycle_
        (coeffsDict(type).getOrDefault<label>("lastAcceleratedCycle", 20)),
    maxDVChange_(nullptr),
    strictMaxDVChange_
        (coeffsDict(type).getOrDefault<bool>("strictMaxDVChange", false)),
    targetConstraintReduction_
    (
        coeffsDict(type).
            getOrDefault<scalar>("targetConstraintReduction", 0.5)
    ),
    bcMult_(coeffsDict(type).getOrDefault<scalar>("boundConstraintsMult", 0.5))
{
    if (coeffsDict(type).found("maxDVChange"))
    {
        maxDVChange_.reset
            (new scalar(coeffsDict(type).get<scalar>("maxDVChange")));
    }
    else if (designVars_().isMaxInitChangeSet())
    {
        maxDVChange_.reset(new scalar(designVars_().getMaxInitChange()()));
    }
    else if (adaptiveStep_)
    {
        FatalErrorInFunction
            << "Either maxInitChange or maxDVChange has to be setup when using "
            << "an adaptive step"
            << exit(FatalError);
    }

    // Sanity checks for the bounds of the design variables.
    // If all design variables start from either their lower or upper bound
    // and there is at least one active flow-related constraint, the systems
    // providing the update of the desgin varaibles become singular.
    // Issue a warning and change the initial values of the design variables
    // by a small ammount, to kick start the optimisation
    if (designVars_->lowerBounds() && designVars_->upperBounds())
    {
        const scalarField& lowerBounds = designVars_->lowerBoundsRef();
        const scalarField& upperBounds = designVars_->upperBoundsRef();
        scalarField& designVars = designVars_();
        boolList isOnLowerBound(designVars.size(), false);
        boolList isOnUpperBound(designVars.size(), false);
        bool existsNonBoundVar = false;
        for (const label activeVarI : designVars_->activeDesignVariables())
        {
            isOnLowerBound[activeVarI] =
                mag(designVars[activeVarI] - lowerBounds[activeVarI]) < SMALL;
            isOnUpperBound[activeVarI] =
                mag(designVars[activeVarI] - upperBounds[activeVarI]) < SMALL;
            existsNonBoundVar =
                existsNonBoundVar
            || (!isOnLowerBound[activeVarI] && !isOnUpperBound[activeVarI]);
        }
        if (globalSum_)
        {
            reduce(existsNonBoundVar, orOp<bool>());
        }
        if (!existsNonBoundVar)
        {
            WarningInFunction
                << "All design variables lay on their bound values." << nl
                << tab << "This will lead to singular matrices in nullSpace"
                << nl
                << tab << "Perturbing the design variables by "
                << epsPerturb_ << "*(upperBound - lowerBounds)"
                << nl
                << endl;
            scalarField update(designVars.size(), Zero);
            for (const label activeVarI : designVars_->activeDesignVariables())
            {
                if (isOnLowerBound[activeVarI])
                {
                    update[activeVarI] =
                        epsPerturb_
                       *(upperBounds[activeVarI] - lowerBounds[activeVarI]);
                }
                if (isOnUpperBound[activeVarI])
                {
                    update[activeVarI] =
                      - epsPerturb_
                       *(upperBounds[activeVarI] - lowerBounds[activeVarI]);
                }
            }
            designVars_->update(update);
        }
    }

    // Read-in aJ in restarts
    dictionary::readIfPresent("aJ", aJ_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nullSpace::computeCorrection()
{
    // Update sizes of fields related to the active or saturated
    // constraints and initialise values
    initialise();

    // Solve the dual problem using a Newton optimiser
    solveDualProblem();

    // Update subsets related to the null space and range space updates
    updateNullAndRangeSpaceSubsets();

    // Compute the update of the design variables
    correction();

    // Increase counter
    ++counter_;
}


Foam::scalar Foam::nullSpace::computeMeritFunction()
{
    scalar LagrangePart = objectiveValue_;
    const autoPtr<scalarField>& lBounds = designVars_->lowerBounds();
    const autoPtr<scalarField>& uBounds = designVars_->upperBounds();
    forAll(iTildaEps_[0], i)
    {
        LagrangePart += mu_[i]*cValues_[iTildaEps_[0][i]];
    }
    scalar lowerPart = Zero;
    forAll(iTildaEps_[1], i)
    {
        const label iActive = activeDesignVars_[iTildaEps_[1][i]];
        lowerPart += l_[i]*(lBounds()[iActive] - designVars_()[iActive]);
    }
    scalar upperPart = Zero;
    forAll(iTildaEps_[2], i)
    {
        const label iActive = activeDesignVars_[iTildaEps_[2][i]];
        upperPart += u_[i]*(designVars_()[iActive] - uBounds()[iActive]);
    }
    if (globalSum_)
    {
        reduce(lowerPart, sumOp<scalar>());
        reduce(upperPart, sumOp<scalar>());
    }
    LagrangePart += lowerPart + upperPart;
    LagrangePart *= aJ_;

    scalarField g(activeConstraints(iTildaEps_));
    scalarField invAAT_g(constraintRelatedUpdate(g, iTildaEps_));

    scalar constraintPart = Zero;
    const label gSize = g.size();
    const label flowSize = iTildaEps_[0].size();
    forAll(iTildaEps_[0], i)
    {
        constraintPart += g[i]*invAAT_g[i];
    }

    scalarField a = SubField<scalar>(g, gSize - flowSize, flowSize);
    scalarField b = SubField<scalar>(invAAT_g, gSize - flowSize, flowSize);
    constraintPart += globalSum(a*b);
    constraintPart *= 0.5*aC_;

    return LagrangePart + constraintPart;
}


Foam::scalar Foam::nullSpace::meritFunctionDirectionalDerivative()
{
    scalarField LagrangeMults(mu_.size() + l_.size() + u_.size(), Zero);
    SubField<scalar>(LagrangeMults, mu_.size(), 0) = mu_;
    SubField<scalar>(LagrangeMults, l_.size(), mu_.size()) = l_;
    SubField<scalar>(LagrangeMults, u_.size(), mu_.size() + l_.size()) = u_;
    scalarField deriv(objectiveDerivatives_, activeDesignVars_);
    deriv += ATv(LagrangeMults, iTildaEps_);
    deriv *= aJ_;

    scalarField g(activeConstraints(iTildaEps_));
    scalarField lamdaC(constraintRelatedUpdate(g, iTildaEps_));
    scalarField ksiC(ATv(lamdaC, iTildaEps_));

    deriv += aC_*ksiC;

    scalarField derivAll(designVars_().size(), Zero);
    derivAll.rmap(deriv, activeDesignVars_);

    return globalSum(sqr(derivAll*correction_));
}


bool Foam::nullSpace::writeData(Ostream& os) const
{
    os.writeEntry("aJ", aJ_);
    return updateMethod::writeData(os);
}


// ************************************************************************* //
