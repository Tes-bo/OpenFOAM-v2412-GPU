/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "ReynoldsAnalogy.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{
    defineTypeNameAndDebug(ReynoldsAnalogy, 0);
    addToRunTimeSelectionTable
    (
        heatTransferCoeffModel,
        ReynoldsAnalogy,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::rho(const label patchi) const
{
    if (rhoName_ == "rhoInf")
    {
        const label n = mesh_.boundary()[patchi].size();
        return tmp<scalarField>::New(n, rhoRef_);
    }
    else if (mesh_.foundObject<volScalarField>(rhoName_, false))
    {
        const auto& rho = mesh_.lookupObject<volScalarField>(rhoName_);
        return rho.boundaryField()[patchi];
    }

    FatalErrorInFunction
        << "Unable to set rho for patch " << patchi
        << exit(FatalError);

    return nullptr;
}


Foam::tmp<Foam::scalarField>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::Cp(const label patchi) const
{
    if (CpName_ == "CpInf")
    {
        const label n = mesh_.boundary()[patchi].size();
        return tmp<scalarField>::New(n, CpRef_);
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const scalarField& pp = thermo.p().boundaryField()[patchi];
        const scalarField& Tp = thermo.T().boundaryField()[patchi];

        return thermo.Cp(pp, Tp, patchi);
    }

    FatalErrorInFunction
        << "Unable to set Cp for patch " << patchi
        << exit(FatalError);

    return nullptr;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::devReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff()/turb.rho();
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return turb.devReff();
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const auto& U = mesh_.lookupObject<volVectorField>(UName_);

        return -thermo.nu()*devTwoSymm(fvc::grad(U));
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            mesh_.lookupObject<transportModel>("transportProperties");

        const auto& U = mesh_.lookupObject<volVectorField>(UName_);

        return -laminarT.nu()*devTwoSymm(fvc::grad(U));
    }
    else if (mesh_.foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            mesh_.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        const auto& U = mesh_.lookupObject<volVectorField>(UName_);

        return -nu*devTwoSymm(fvc::grad(U));
    }

    FatalErrorInFunction
        << "No valid model for viscous stress calculation"
        << exit(FatalError);

    return nullptr;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::Cf() const
{
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const volVectorField::Boundary& Ubf = U.boundaryField();

    auto tCf = tmp<FieldField<Field, scalar>>::New(Ubf.size());
    auto& Cf = tCf.ref();

    forAll(Cf, patchi)
    {
        Cf.set(patchi, new Field<scalar>(Ubf[patchi].size(), Zero));
    }

    const volSymmTensorField R(devReff());
    const volSymmTensorField::Boundary& Rbf = R.boundaryField();

    for (const label patchi : patchIDs_)
    {
        const fvPatchVectorField& Up = Ubf[patchi];

        const symmTensorField& Rp = Rbf[patchi];

        tmp<vectorField> tnHat = Up.patch().nf();

        tmp<scalarField> ttauByRhop = mag(tnHat & Rp);

        Cf[patchi] = 2*ttauByRhop/magSqr(URef_);
    }

    return tCf;
}


void Foam::heatTransferCoeffModels::ReynoldsAnalogy::htc
(
    volScalarField& htc,
    const FieldField<Field, scalar>& q
)
{
    const FieldField<Field, scalar> CfBf(Cf());
    const scalar magU = mag(URef_);

    volScalarField::Boundary& htcBf = htc.boundaryFieldRef();

    for (const label patchi : patchIDs_)
    {
        tmp<scalarField> trhop = rho(patchi);
        tmp<scalarField> tCpp = Cp(patchi);

        htcBf[patchi] = 0.5*trhop*tCpp*magU*CfBf[patchi];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModels::ReynoldsAnalogy::ReynoldsAnalogy
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    heatTransferCoeffModel(dict, mesh, TName),
    UName_("U"),
    URef_(Zero),
    rhoName_("rho"),
    rhoRef_(0),
    CpName_("Cp"),
    CpRef_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModels::ReynoldsAnalogy::read
(
    const dictionary& dict
)
{
    if (!heatTransferCoeffModel::read(dict))
    {
        return false;
    }

    dict.readIfPresent("U", UName_);
    dict.readEntry("UInf", URef_);

    dict.readIfPresent("Cp", CpName_);
    if (CpName_ == "CpInf")
    {
        dict.readEntry("CpInf", CpRef_);
    }

    dict.readIfPresent("rho", rhoName_);
    if (rhoName_ == "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }

    return true;
}


// ************************************************************************* //
