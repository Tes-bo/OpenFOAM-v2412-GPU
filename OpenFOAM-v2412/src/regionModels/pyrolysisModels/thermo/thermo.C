/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd
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

#include "thermo.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "absorptionEmissionModel.H"
#include "fvm.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermo, 0);
addToRunTimeSelectionTable(pyrolysisModel, thermo, mesh);
addToRunTimeSelectionTable(pyrolysisModel, thermo, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void thermo::readControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.readEntry("nNonOrthCorr", nNonOrthCorr_);
    time().controlDict().readEntry("maxDi", maxDiff_);
}


bool thermo::read()
{
    if (pyrolysisModel::read())
    {
        readControls();
        return true;
    }

    return false;
}


bool thermo::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readControls();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermo::thermo
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    solidThermo_(solidThermo::New(regionMesh())),
    radiation_(radiation::radiationModel::New(solidThermo_->T())),
    nNonOrthCorr_(-1),
    maxDiff_(10)
{
    if (active())
    {
         readControls();
    }
}


thermo::thermo
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    solidThermo_(solidThermo::New(regionMesh())),
    radiation_(radiation::radiationModel::New(solidThermo_->T())),
    nNonOrthCorr_(-1),
    maxDiff_(10)
{
    if (active())
    {
        readControls();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermo::~thermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermo::preEvolveRegion()
{
     if (active())
     {
        pyrolysisModel::preEvolveRegion();
     }
}


void thermo::evolveRegion()
{
    if (active())
    {
        Info<< "\nEvolving energy in region: " << regionMesh().name() << endl;

        volScalarField& h = solidThermo_->he();

        for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
        {
            tmp<volScalarField> alpha(solidThermo_->alpha());
            fvScalarMatrix hEqn
            (
                fvm::ddt(rho(), h)
              - fvm::laplacian(alpha, h)
              + fvc::laplacian(alpha, h)
              - fvc::laplacian(solidThermo_->kappa(), T())
            );

            hEqn.relax();
            hEqn.solve();
        }

        solidThermo_->correct();
    }
}


const volScalarField& thermo::rho() const
{
    return solidThermo_->rho();
}


const volScalarField& thermo::T() const
{
    return solidThermo_->T();
}


const tmp<volScalarField> thermo::Cp() const
{
    return solidThermo_->Cp();
}


tmp<volScalarField> thermo::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> thermo::kappa() const
{
    return solidThermo_->kappa();
}


const surfaceScalarField& thermo::phiGas() const
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
