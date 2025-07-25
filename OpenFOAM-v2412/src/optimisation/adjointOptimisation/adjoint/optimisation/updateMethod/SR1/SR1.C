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

#include "SR1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SR1, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        SR1,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::SR1::updateHessian()
{
    // Vectors needed to construct the inverse HessianInv matrix
    scalarField y(activeDesignVars_.size(), Zero);
    scalarField s(activeDesignVars_.size(), Zero);
    y.map(objectiveDerivatives_ - derivativesOld_, activeDesignVars_);
    s.map(correctionOld_, activeDesignVars_);

    scalarField temp(s - rightMult(Hessian_(), y));

    // Construct the inverse HessianInv
    scalar tempMag = sqrt(globalSum(sqr(temp)));
    scalar yMag = sqrt(globalSum(sqr(y)));
    scalar HessYMag = sqrt(globalSum(sqr(rightMult(Hessian_(), y))));

    // Stability check
    if (tempMag > ratioThreshold_ * yMag * HessYMag)
    {
        Hessian_() += (scalar(1)/(globalSum(temp*y)))*outerProd(temp, temp);
    }
    else
    {
        WarningInFunction
            << "Denominator of update too small. Keeping old Hessian" << endl;
    }
}


void Foam::SR1::update()
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
    else
    {
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

Foam::SR1::SR1
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    quasiNewton(mesh, dict, designVars, nConstraints, type),
    ratioThreshold_
    (
        coeffsDict(type).getOrDefault<scalar>("ratioThreshold", 1e-08)
    )
{
    allocateHessian();
}


// ************************************************************************* //
