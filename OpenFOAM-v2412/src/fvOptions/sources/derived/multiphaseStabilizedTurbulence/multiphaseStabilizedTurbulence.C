/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "multiphaseStabilizedTurbulence.H"
#include "fvMatrices.H"
#include "turbulentTransportModel.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphaseStabilizedTurbulence, 0);
    addToRunTimeSelectionTable
    (
        option,
        multiphaseStabilizedTurbulence,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphaseStabilizedTurbulence::multiphaseStabilizedTurbulence
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            coeffs_,
            0.09
        )
    ),
    C_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C",
            coeffs_,
            1.51
        )
    ),
    lambda2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "lambda2",
            coeffs_,
            0.1
        )
    ),
    alpha_
    (
        dimensionedScalar::getOrAddToDict
        (
            "alpha",
            coeffs_,
            1.36
        )
    )
{
    fieldNames_.resize(2);

    // Note: incompressible only
    const auto* turbPtr =
        mesh_.findObject<incompressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (turbPtr)
    {
        const tmp<volScalarField>& tk = turbPtr->k();
        fieldNames_[0] = tk().name();

        const tmp<volScalarField>& tnut = turbPtr->nut();
        fieldNames_[1] = tnut().name();

        Log << "    Applying model to " << fieldNames_[0]
            << " and " << fieldNames_[1] << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find incompressible turbulence model"
            << exit(FatalError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphaseStabilizedTurbulence::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // Not applicable to compressible cases
    NotImplemented;
}


void Foam::fv::multiphaseStabilizedTurbulence::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi != 0)
    {
        return;
    }

    Log << this->name() << ": applying buoyancy production term to "
        << eqn.psi().name() << endl;

    // Buoyancy production in k eqn

    const auto* turbPtr =
        mesh_.findObject<incompressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find incompressible turbulence model"
            << exit(FatalError);
    }


    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    const volScalarField& epsilon = tepsilon();
    const volScalarField& k = eqn.psi();

    // Note: using solver density field for incompressible multiphase cases
    const auto& rho = mesh_.lookupObject<volScalarField>(rhoName_);

    const auto& g = meshObjects::gravity::New(mesh_.time());

    const dimensionedScalar eps0("eps0", epsilon.dimensions(), SMALL);

    // Note: differing from reference by replacing nut/k by Cmu*k/epsilon
    const volScalarField GbyK
    (
        "GbyK",
        Cmu_*k/(epsilon + eps0)*alpha_*(g & fvc::grad(rho))/rho
    );

    return eqn -= fvm::SuSp(GbyK, k);
}


void Foam::fv::multiphaseStabilizedTurbulence::correct(volScalarField& field)
{
    if (field.name() != fieldNames_[1])
    {
        return;
    }

    Log << this->name() << ": correcting " << field.name() << endl;

    const auto* turbPtr =
        mesh_.findObject<incompressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    // nut correction
    const auto& U = turbPtr->U();
    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    const auto& epsilon = tepsilon();
    tmp<volScalarField> tk = turbPtr->k();
    const auto& k = tk();

    tmp<volTensorField> tgradU = fvc::grad(U);
    const auto& gradU = tgradU();
    const dimensionedScalar pSmall("pSmall", dimless/sqr(dimTime), SMALL);
    const volScalarField pRat
    (
        magSqr(symm(gradU))/(magSqr(skew(gradU)) + pSmall)
    );

    const volScalarField epsilonTilde
    (
        max
        (
            epsilon,
            lambda2_*C_*pRat*epsilon
        )
    );

    field = Cmu_*sqr(k)/epsilonTilde;
    field.correctBoundaryConditions();
}


// ************************************************************************* //
