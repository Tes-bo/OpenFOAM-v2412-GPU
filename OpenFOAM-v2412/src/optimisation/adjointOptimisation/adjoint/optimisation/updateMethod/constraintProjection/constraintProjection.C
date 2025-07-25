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

#include "constraintProjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constraintProjection, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        constraintProjection,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        constraintProjection,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constraintProjection::constraintProjection
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
    useCorrection_(coeffsDict(type).getOrDefault<bool>("useCorrection", true)),
    delta_(coeffsDict(type).getOrDefault<scalar>("delta", 0.1))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constraintProjection::computeCorrection()
{
    // Reset to zero
    const label n = objectiveDerivatives_.size();
    const label m = constraintDerivatives_.size();

    // Matrix with constraint derivatives and its inverse
    scalarSquareMatrix MMT(m, Zero);
    forAll(constraintDerivatives_, cI)
    {
        forAll(constraintDerivatives_, cJ)
        {
            MMT[cI][cJ] =
                globalSum
                (
                    constraintDerivatives_[cI] * constraintDerivatives_[cJ]
                );
        }
    }
    scalarSquareMatrix invM(inv(MMT));

    // Contribution from constraints
    scalarField constraintContribution(n, Zero);
    scalarField nonLinearContribution(n, Zero);
    forAll(constraintDerivatives_, cI)
    {
        forAll(constraintDerivatives_, cJ)
        {
            constraintContribution +=
                constraintDerivatives_[cI]
               *invM[cI][cJ]
               *globalSum(constraintDerivatives_[cJ] * objectiveDerivatives_);

            // Correction to take non-linearities into consideration
            if (useCorrection_)
            {
                nonLinearContribution +=
                    constraintDerivatives_[cI]
                   *invM[cI][cJ]
                   *cValues_[cJ];
            }
        }
    }

    // Final correction
    scalarField correction(objectiveDerivatives_ - constraintContribution);
    correction *= -eta_;
    if (useCorrection_)
    {
        correction -= nonLinearContribution;
    }

    for (const label varI : activeDesignVars_)
    {
        correction_[varI] = correction[varI];
    }
}


Foam::scalar Foam::constraintProjection::computeMeritFunction()
{
    return objectiveValue_ + delta_*sum(mag(cValues_));
}


// ************************************************************************* //
