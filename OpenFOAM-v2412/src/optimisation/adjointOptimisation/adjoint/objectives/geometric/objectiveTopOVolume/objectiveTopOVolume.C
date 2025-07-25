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

#include "objectiveTopOVolume.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveTopOVolume, 1);
addToRunTimeSelectionTable
(
    objectiveGeometric,
    objectiveTopOVolume,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveTopOVolume::objectiveTopOVolume
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveGeometric(mesh, dict, adjointSolverName, primalSolverName),
    targetPercentage_(Function1<scalar>::New("percentage", dict)),
    percentInDenom_(dict.getOrDefault<bool>("percentInDenom", true))
{
    // Allocate boundary field pointers
    dJdbPtr_.reset(createZeroFieldPtr<scalar>(mesh_, "dJdb", dimless));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar objectiveTopOVolume::J()
{
    J_ = Zero;
    if (mesh_.foundObject<volScalarField>("beta"))
    {
        const volScalarField& beta = mesh_.lookupObject<volScalarField>("beta");
        const DimensionedField<scalar, volMesh>& V = mesh_.V();
        const scalar time = mesh_.time().timeOutputValue();
        J_ =
            scalar(1) - gSum(beta.primitiveField()*V)/gSum(V)
          - targetPercentage_->value(time);
        if (percentInDenom_)
        {
            J_ /= targetPercentage_->value(time);
        }
    }
    else
    {
        WarningInFunction
            << "Beta field not yet registered in database. OK for start-up"
            << endl;
    }
    return J_;
}


void objectiveTopOVolume::update_dJdb()
{
    const scalar time = mesh_.time().timeOutputValue();
    dJdbPtr_().primitiveFieldRef() = -scalar(1)/gSum(mesh_.V());
    if (percentInDenom_)
    {
        dJdbPtr_().primitiveFieldRef() /= targetPercentage_->value(time);
    }
}


void objectiveTopOVolume::addHeaderColumns() const
{
    objFunctionFilePtr_()
        << setw(width_) << "TargetVolume" <<  " ";
}


void objectiveTopOVolume::addColumnValues() const
{
    const scalar time = mesh_.time().timeOutputValue();
    objFunctionFilePtr_()
        << setw(width_) << targetPercentage_->value(time) << " ";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
