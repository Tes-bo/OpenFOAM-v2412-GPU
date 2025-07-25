/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "viscousDissipation.H"
#include "fvMatrices.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(viscousDissipation, 0);
    addToRunTimeSelectionTable(option, viscousDissipation, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::viscousDissipation::rho() const
{
    auto trho = volScalarField::New
    (
        "trho",
        IOobject::NO_REGISTER,
        mesh_,
        rho_
    );

    if (rho_.value() > 0)
    {
        return trho;
    }
    else if (rhoName_ != "none")
    {
        trho.ref() = mesh_.lookupObject<volScalarField>(rhoName_);
        return trho;
    }

    FatalErrorInFunction
        << "Neither rhoName nor rho are specified."
        << exit(FatalError);

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::viscousDissipation::viscousDissipation
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    rhoName_(coeffs_.getOrDefault<word>("rho", "none")),
    rho_
    (
        coeffs_.getOrDefault
        (
            "rhoInf",
            dimensionedScalar(dimDensity, 0)
        )
    )
{
    const auto* thermoPtr =
        mesh_.findObject<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        fieldNames_.resize(1, thermoPtr->he().name());
    }

    if (fieldNames_.empty())
    {
        coeffs_.readEntry("fields", fieldNames_);
    }

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::fv::viscousDissipation::devRhoReff() const
{
    // Incompressible
    {
        const auto* turbPtr =
            mesh_.findObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(rho()*turbPtr->devRhoReff());
        }
    }

    // Compressible
    {
        const auto* turbPtr =
            mesh_.findObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(turbPtr->devRhoReff());
        }
    }

    FatalErrorInFunction
        << " The turbulence model is not found in the database."
        << exit(FatalError);

    return nullptr;
}


void Foam::fv::viscousDissipation::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    typedef typename outerProduct<vector, vector>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const word gradUName("grad(" + UName_ + ')');

    auto tgradU = GradFieldType::New
    (
        "gradU",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedTensor(inv(dimTime), Zero)
    );

    // Cached?
    const auto* gradUPtr = mesh_.findObject<GradFieldType>(gradUName);

    if (gradUPtr)
    {
        tgradU.ref() = *gradUPtr;
    }
    else
    {
        const auto& U = mesh_.lookupObject<volVectorField>(UName_);
        tgradU.ref() = fvc::grad(U);
    }

    const volScalarField D("D", devRhoReff() && tgradU.ref());

    eqn -= D;
}


// ************************************************************************* //
