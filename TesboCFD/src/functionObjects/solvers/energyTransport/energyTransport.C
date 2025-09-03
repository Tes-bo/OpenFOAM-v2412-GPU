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

#include "energyTransport.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(energyTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        energyTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::energyTransport::transportedField()
{
    if (!foundObject<volScalarField>(fieldName_))
    {
        auto tfldPtr = tmp<volScalarField>::New
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_
        );
        store(fieldName_, tfldPtr);
    }

    return lookupObjectRef<volScalarField>(fieldName_);
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::energyTransport::kappaEff() const
{
    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* turbPtr = findObject<turbType>
        (
            turbulenceModel::propertiesName
        );

        if (turbPtr)
        {
            return volScalarField::New
            (
                "kappaEff",
                IOobject::NO_REGISTER,
                (
                    kappa() + Cp()*turbPtr->nut()*rho()/Prt_
                )
            );
        }
    }

    FatalErrorInFunction
        << "Turbulence model not found" << exit(FatalError);

    return nullptr;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::energyTransport::rho() const
{
    auto trho = volScalarField::New
    (
        "trho",
        IOobject::NO_REGISTER,
        mesh_,
        rho_
    );

    if (phases_.size())
    {
        trho.ref() = lookupObject<volScalarField>(rhoName_);
    }
    return trho;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::energyTransport::Cp() const
{
    if (phases_.size())
    {
        tmp<volScalarField> tCp(phases_[0]*Cps_[0]);

        for (label i = 1; i < phases_.size(); i++)
        {
            tCp.ref() += phases_[i]*Cps_[i];
        }
        return tCp;
    }

    return volScalarField::New
    (
        "tCp",
        IOobject::NO_REGISTER,
        mesh_,
        Cp_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::energyTransport::kappa() const
{
    if (phases_.size())
    {
        tmp<volScalarField> tkappa(phases_[0]*kappas_[0]);

        for (label i = 1; i < phases_.size(); i++)
        {
            tkappa.ref() += phases_[i]*kappas_[i];
        }
        return tkappa;
    }

    return volScalarField::New
    (
        "tkappa",
        IOobject::NO_REGISTER,
        mesh_,
        kappa_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::energyTransport::energyTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimTemperature/dimVolume, Zero)
    ),
    fvOptions_(mesh_),
    multiphaseThermo_(dict.subOrEmptyDict("phaseThermos")),
    Cp_("Cp", dimEnergy/dimMass/dimTemperature, 0, dict),
    kappa_
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        0,
        dict
    ),
    rho_("rhoInf", dimDensity, 0, dict),
    Prt_("Prt", dimless, 1, dict),
    fieldName_(dict.getOrDefault<word>("field", "T")),
    schemesField_("unknown-schemesField"),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    tol_(1),
    nCorr_(0)
{
    read(dict);

    // If the flow is multiphase
    if (!multiphaseThermo_.empty())
    {
        Cps_.setSize(multiphaseThermo_.size());
        kappas_.setSize(Cps_.size());
        phaseNames_.setSize(Cps_.size());

        label phasei = 0;
        forAllConstIters(multiphaseThermo_, iter)
        {
            const word& key = iter().keyword();

            const dictionary* subDictPtr = multiphaseThermo_.findDict(key);

            if (!subDictPtr)
            {
                FatalErrorInFunction
                    << "Found non-dictionary entry " << iter()
                    << " in top-level dictionary " << multiphaseThermo_
                    << exit(FatalError);
            }

            const dictionary& dict = *subDictPtr;

            phaseNames_[phasei] = key;

            Cps_.set
            (
                phasei,
                new dimensionedScalar
                (
                    "Cp",
                    dimEnergy/dimMass/dimTemperature,
                    dict
                )
            );

            kappas_.set
            (
                phasei,
                new dimensionedScalar //[J/m/s/K]
                (
                    "kappa",
                    dimEnergy/dimTime/dimLength/dimTemperature,
                    dict
                )
            );

            ++phasei;
        }

        phases_.setSize(phaseNames_.size());
        forAll(phaseNames_, i)
        {
            phases_.set
            (
                i,
                mesh_.getObjectPtr<volScalarField>(phaseNames_[i])
            );
        }

        rhoCp_ = rho()*Cp();
        rhoCp_.oldTime();
    }
    else
    {
        if (Cp_.value() == 0.0 || kappa_.value() == 0.0)
        {
            FatalErrorInFunction
                << " Multiphase thermo dictionary not found and Cp/kappa "
                << " for single  phase are zero. Please entry either"
                << exit(FatalError);
        }

    }

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& s = transportedField();
    s.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::energyTransport::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    dict.readIfPresent("phi", phiName_);
    dict.readIfPresent("rho", rhoName_);

    schemesField_ = dict.getOrDefault("schemesField", fieldName_);

    dict.readIfPresent("nCorr", nCorr_);
    dict.readIfPresent("tolerance", tol_);

    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    return true;
}


bool Foam::functionObjects::energyTransport::execute()
{
    volScalarField& s = transportedField();

    Log << type() << " execute: " << s.name() << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Calculate the diffusivity
    const volScalarField kappaEff("kappaEff", this->kappaEff());

    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme
    (
        "laplacian(kappaEff," + schemesField_ + ")"
    );

    // Set under-relaxation coeff
    scalar relaxCoeff = 0;
    mesh_.relaxEquation(schemesField_, relaxCoeff);

    // Convergence monitor parameters
    bool converged = false;
    int iter = 0;

    if (phi.dimensions() == dimMass/dimTime)
    {
        rhoCp_ = rho()*Cp();
        const surfaceScalarField rhoCpPhi(fvc::interpolate(Cp())*phi);

        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rhoCp_, s)
              + fvm::div(rhoCpPhi, s, divScheme)
              - fvm::Sp(fvc::ddt(rhoCp_) + fvc::div(rhoCpPhi), s)
              - fvm::laplacian(kappaEff, s, laplacianScheme)
             ==
                fvOptions_(rhoCp_, s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            ++iter;
            converged = (sEqn.solve(schemesField_).initialResidual() < tol_);
            if (converged) break;
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        dimensionedScalar rhoCp(rho_*Cp_);

        const surfaceScalarField CpPhi(rhoCp*phi);

        auto trhoCp = volScalarField::New
        (
            "trhoCp",
            IOobject::NO_REGISTER,
            mesh_,
            rhoCp
        );

        for (int i = 0; i <= nCorr_; ++i)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rhoCp, s)
              + fvm::div(CpPhi, s, divScheme)
              - fvm::laplacian(kappaEff, s, laplacianScheme)
             ==
                fvOptions_(trhoCp.ref(), s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            ++iter;
            converged = (sEqn.solve(schemesField_).initialResidual() < tol_);
            if (converged) break;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    if (converged)
    {
        Log << type() << ": " << name() << ": "
            << s.name() << " is converged." << nl
            << tab << "initial-residual tolerance: " << tol_ << nl
            << tab << "outer iteration: " << iter << nl;
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::energyTransport::write()
{
    return true;
}


// ************************************************************************* //
