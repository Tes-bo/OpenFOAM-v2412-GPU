/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "noPyrolysis.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noPyrolysis, 0);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, mesh);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void noPyrolysis::constructThermoChemistry()
{
    solidThermo_.reset
    (
        solidReactionThermo::New(regionMesh()).ptr()
    );

    solidChemistry_.reset
    (
        basicSolidChemistryModel::New(solidThermo_()).ptr()
    );

    radiation_.reset(radiation::radiationModel::New
    (
        solidChemistry_->solidThermo().T()
    ).ptr());
}


bool noPyrolysis::read()
{
    // No additional info to read
    return pyrolysisModel::read();
}


bool noPyrolysis::read(const dictionary& dict)
{
    // No additional info to read
    return pyrolysisModel::read(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPyrolysis::noPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(mesh, regionType),
    solidThermo_(nullptr),
    solidChemistry_(nullptr),
    radiation_(nullptr)
{
    if (active())
    {
        constructThermoChemistry();
    }
}


noPyrolysis::noPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(mesh, regionType),
    solidThermo_(nullptr),
    solidChemistry_(nullptr),
    radiation_(nullptr)
{
    if (active())
    {
        constructThermoChemistry();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noPyrolysis::~noPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noPyrolysis::preEvolveRegion()
{}


void noPyrolysis::evolveRegion()
{}


const volScalarField& noPyrolysis::rho() const
{
    return solidChemistry_->solidThermo().rho();
}


const volScalarField& noPyrolysis::T() const
{
    return solidChemistry_->solidThermo().T();
}


const tmp<volScalarField> noPyrolysis::Cp() const
{
    return solidChemistry_->solidThermo().Cp();
}


tmp<volScalarField> noPyrolysis::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> noPyrolysis::kappa() const
{
    return solidChemistry_->solidThermo().kappa();
}


const surfaceScalarField& noPyrolysis::phiGas() const
{
    FatalErrorInFunction
        << "phiGas field not available for " << type() << abort(FatalError);
    return surfaceScalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
