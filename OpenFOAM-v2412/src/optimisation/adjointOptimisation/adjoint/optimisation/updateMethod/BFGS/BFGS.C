/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "BFGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BFGS, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        BFGS,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::BFGS::updateHessian()
{
    // Vectors needed to construct the inverse HessianInv matrix
    scalarField y(activeDesignVars_.size(), Zero);
    scalarField s(activeDesignVars_.size(), Zero);
    y.map(objectiveDerivatives_ - derivativesOld_, activeDesignVars_);
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
                Hessian_()[varI][varI] *= scaleFactor;
            }
        }
        else
        {
            WarningInFunction
                << "y*s is negative. Skipping the scaling of the first Hessian"
                << endl;
        }
    }

    // Check curvature condition
    if (ys < curvatureThreshold_)
    {
        WarningInFunction
            << " y*s is below threshold! y*s=" << ys << endl;
    }

    DebugInfo
        << "Hessian curvature index " << ys << endl;

    // Construct the inverse HessianInv
    Hessian_() +=
        (ys + globalSum(leftMult(y, Hessian_())*y))/sqr(ys)*outerProd(s, s)
      - (scalar(1)/ys)*
        (
           outerProd(rightMult(Hessian_(), y), s)
         + outerProd(s, leftMult(y, Hessian_()))
        );
}


void Foam::BFGS::update()
{
    // In the first few iterations, use steepest descent but update the Hessian
    // matrix
    if (counter_ < nSteepestDescent_)
    {
        Info<< "Using steepest descent to update design variables" << endl;
        for (const label varI : activeDesignVars_)
        {
            correction_[varI] = -eta_*objectiveDerivatives_[varI];
        }
    }
    // Else use BFGS formula to update the design variables
    else
    {
        // Compute correction for active design variables
        scalarField activeDerivs(activeDesignVars_.size(), Zero);
        activeDerivs.map(objectiveDerivatives_, activeDesignVars_);
        scalarField activeCorrection
        (
            -etaHessian_*rightMult(Hessian_(), activeDerivs)
        );

        // Transfer correction to the global list
        correction_ = Zero;
        forAll(activeDesignVars_, varI)
        {
            correction_[activeDesignVars_[varI]] = activeCorrection[varI];
        }
    }

    // Store fields for the next iteration
    derivativesOld_ = objectiveDerivatives_;
    correctionOld_ = correction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BFGS::BFGS
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    quasiNewton(mesh, dict, designVars, nConstraints, type),
    curvatureThreshold_
    (
        coeffsDict(type).getOrDefault<scalar>("curvatureThreshold", 1e-10)
    )
{
    allocateHessian();
}


// ************************************************************************* //
