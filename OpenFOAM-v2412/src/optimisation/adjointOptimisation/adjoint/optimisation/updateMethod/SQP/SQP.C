/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "SQP.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SQP, 1);
    addToRunTimeSelectionTable
    (
        updateMethod,
        SQP,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        SQP,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SQP::updateHessian()
{
    // Vectors needed to construct the (inverse) Hessian matrix
    scalarField y(activeDesignVars_.size(), Zero);
    scalarField s(activeDesignVars_.size(), Zero);
    scalarField LagrangianDerivativesOld = derivativesOld_;
    forAll(constraintDerivatives_, cI)
    {
        LagrangianDerivatives_ -= lamdas_[cI] * constraintDerivatives_[cI];
        LagrangianDerivativesOld -= lamdas_[cI] * constraintDerivativesOld_[cI];
    }
    y.map(LagrangianDerivatives_ - LagrangianDerivativesOld, activeDesignVars_);
    s.map(correctionOld_, activeDesignVars_);

    scalar ys = globalSum(s*y);
    if (counter_ == 1 && scaleFirstHessian_)
    {
        if (ys > scalar(0))
        {
            scalar scaleFactor = ys/globalSum(y*y);
            Info<< "Scaling Hessian with factor " << scaleFactor << endl;
            forAll(activeDesignVars_, varI)
            {
                Hessian_()[varI][varI] /= scaleFactor;
            }
        }
        else
        {
            WarningInFunction
                << " y*s is negative. Skipping the scaling of the first Hessian"
                << endl;
        }
    }
    scalar sBs = globalSum(leftMult(s, Hessian_())*s);

    // Check curvature condition
    scalar theta(1);
    if (ys < dumpingThreshold_*sBs)
    {
        WarningInFunction
            << " y*s is below threshold. Using damped form" << endl;
        theta = (1 - dumpingThreshold_)*sBs/(sBs - ys);
    }
    scalarField r(theta*y + (scalar(1) - theta)*rightMult(Hessian_(), s));
    DebugInfo
        << "Unmodified Hessian curvature index " << ys << endl;
    DebugInfo
        << "Modified Hessian curvature index " << globalSum(r*s) << endl;

    // Update the Hessian
    Hessian_() +=
      - outerProd(rightMult(Hessian_(), s), leftMult(s/sBs, Hessian_()))
      + outerProd(r, r/globalSum(s*r));
}


void Foam::SQP::update()
{
    // Also denoted below as W
    SquareMatrix<scalar> HessianInv = inv(Hessian_());
    if (debug > 1)
    {
        Info<< "Hessian " << Hessian_() << endl;
        Info<< "HessianInv " << HessianInv << endl;
        label n = Hessian_().n();
        SquareMatrix<scalar> test(n, Zero);
        for (label k = 0; k < n; k++)
        {
            for (label l = 0; l < n; l++)
            {
                scalar elem(Zero);
                for (label i = 0; i < n; i++)
                {
                    elem += Hessian_()[k][i] * HessianInv[i][l];
                }
                test[k][l]=elem;
            }
        }
        Info<< "Validation " << test << endl;
    }

    // Compute new Lagrange multipliers
    label nc = constraintDerivatives_.size();
    scalarField activeDerivs(activeDesignVars_.size(), Zero);

    // activeDerivs.map(objectiveDerivatives_, activeDesignVars_);
    activeDerivs.map(LagrangianDerivatives_, activeDesignVars_);
    scalarField WgradL = rightMult(HessianInv, activeDerivs);

    scalarField lamdaRHS(nc, Zero);
    forAll(lamdaRHS, cI)
    {
        scalarField activeConsDerivs(activeDesignVars_.size(), Zero);
        activeConsDerivs.map(constraintDerivatives_[cI], activeDesignVars_);
        lamdaRHS[cI] = globalSum(activeConsDerivs * WgradL) - cValues_[cI];
        if (debug > 1)
        {
            Info<< "lamdaRHS total|deriv part|constraint part "
                << lamdaRHS[cI] << " " << globalSum(activeConsDerivs * WgradL)
                << " " << cValues_[cI] << endl;
        }
    }

    // lhs for the lamda system
    SquareMatrix<scalar> AWA(nc, Zero);
    PtrList<scalarField> WA(nc);
    for (label j = 0; j < nc; j++)
    {
        scalarField gradcJ(activeDesignVars_.size(), Zero);
        gradcJ.map(constraintDerivatives_[j], activeDesignVars_);
        WA.set(j, new scalarField(rightMult(HessianInv, gradcJ)));
        for (label i = 0; i < nc; i++)
        {
            scalarField gradcI(activeDesignVars_.size(), Zero);
            gradcI.map(constraintDerivatives_[i], activeDesignVars_);
            AWA[i][j] = globalSum(gradcI * WA[j]);
        }
    }
    SquareMatrix<scalar> invAWA = inv(AWA);
    scalarField deltaLamda = rightMult(invAWA, lamdaRHS);
    if (debug > 1)
    {
        Info<< "AWA " << AWA << endl;
        Info<< "AWAInv " << invAWA << endl;
        Info<< "lamda update " << deltaLamda << endl;
    }
    lamdas_ += deltaLamda;

    // Compute design variables correction
    scalarField activeCorrection(-WgradL);
    forAll(WA, cI)
    {
        activeCorrection += WA[cI]*deltaLamda[cI];
    }
    activeCorrection *= etaHessian_;
    // Transfer correction to the global list
    correction_ = Zero;
    forAll(activeDesignVars_, varI)
    {
        correction_[activeDesignVars_[varI]] = activeCorrection[varI];
    }
    if (counter_ == 0)
    {
        correction_ *= eta_;
    }
}


void Foam::SQP::storeOldFields()
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


Foam::scalar Foam::SQP::meritFunctionConstraintPart() const
{
    // Assumes that all constraints are known by all processors
    // What about constraints directly imposed on distributed design variables?
    // These should be met in each iteration of the algorithm, so,
    // most probably, there is no problem
    return sum(mag(cValues_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SQP::SQP
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    quasiNewton(mesh, dict, designVars, nConstraints, type),
    SQPBase(mesh, dict, designVars, *this, type),
    dumpingThreshold_
    (
        coeffsDict(type).getOrDefault<scalar>("dumpingThreshold", 0.2)
    )
{
    allocateHessian();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SQP::computeCorrection()
{
    LagrangianDerivatives_ = objectiveDerivatives_;
    quasiNewton::computeCorrection();

    // Store fields for the next iteration and write them to file
    storeOldFields();
}


Foam::scalar Foam::SQP::computeMeritFunction()
{
    // If condition is not met, update mu value
    if (mu_ < max(mag(lamdas_)) + delta_)
    {
        mu_ = max(mag(lamdas_)) + 2*delta_;
        if (debug > 1)
        {
            Info<< "Updated mu value to " << mu_ << endl;
        }
    }
    scalar L = objectiveValue_ + mu_*sum(mag(cValues_));

    return L;
}


Foam::scalar Foam::SQP::meritFunctionDirectionalDerivative()
{
    scalar deriv =
        globalSum(objectiveDerivatives_*correction_)
      - mu_*sum(mag(cValues_));

    return deriv;
}


bool Foam::SQP::writeData(Ostream& os) const
{
    return quasiNewton::writeData(os) && SQPBase::addToFile(os);
}


bool Foam::SQP::writeAuxiliaryData()
{
    return SQPBase::writeMeritFunction(*this);
}


// ************************************************************************* //
