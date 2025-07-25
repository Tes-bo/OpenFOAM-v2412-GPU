/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "quasiNewton.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quasiNewton, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::quasiNewton::allocateHessian()
{
    Hessian_.reset(new SquareMatrix<scalar>(activeDesignVars_.size(), I));
    // Read in Hessian Matrix or initialiase
    const label nDVs(designVars_().activeDesignVariables().size());
    forAll(designVars_().activeDesignVariables(), iDV)
    {
        if (found("Hessian" + Foam::name(iDV)))
        {
            Hessian_().subColumn(iDV) =
                scalarField("Hessian" + Foam::name(iDV), *this, nDVs);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quasiNewton::quasiNewton
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
        const word& type
)
:
    updateMethod(mesh, dict, designVars, nConstraints, type),
    etaHessian_(coeffsDict(type).getOrDefault<scalar>("etaHessian", 1)),
    nSteepestDescent_
    (
        coeffsDict(type).getOrDefault<label>("nSteepestDescent", 1)
    ),
    scaleFirstHessian_
    (
        coeffsDict(type).getOrDefault<bool>("scaleFirstHessian", false)
    ),
    Hessian_(nullptr),
    derivativesOld_
    (
        readOrZeroField("derivativesOld", objectiveDerivatives_.size())
    ),
    correctionOld_(readOrZeroField("correctionOld", correction_.size()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::quasiNewton::computeCorrection()
{
    // The first iteration uses a unitary Hessian. No need to update
    if (counter_ != 0)
    {
        updateHessian();
    }

    update();
    ++counter_;
}


void Foam::quasiNewton::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


bool Foam::quasiNewton::writeData(Ostream& os) const
{
    if (Hessian_)
    {
        // Matrices cannot be written/read in binary.
        // Circumvent this by writing separate columns as scalarFields
        forAll(designVars_().activeDesignVariables(), iDV)
        {
            Hessian_().subColumn(iDV).operator Field<scalar>().
                writeEntry("Hessian" + Foam::name(iDV), os);
        }
    }
    derivativesOld_.writeEntry("derivativesOld", os);
    correctionOld_.writeEntry("correctionOld", os);

    return updateMethod::writeData(os);
}


// ************************************************************************* //
