/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "P1.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(P1, 0);
        addToRadiationRunTimeSelectionTables(P1);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::P1::P1(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    )
{}


Foam::radiation::P1::P1(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::P1::read()
{
    if (radiationModel::read())
    {
        // nothing to read

        return true;
    }

    return false;
}


void Foam::radiation::P1::calculate()
{
    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0("a0", a_.dimensions(), ROOTVSMALL);

    // Construct diffusion
    const auto tgamma = volScalarField::New
    (
        "gammaRad",
        IOobject::REGISTER,  // used by boundary conditions
        1.0/(3.0*a_ + sigmaEff + a0)
    );
    const auto& gamma = tgamma();

    // Solve G transport equation
    solve
    (
        fvm::laplacian(gamma, G_)
      - fvm::Sp(a_, G_)
     ==
      - 4.0*(e_*physicoChemical::sigma*pow4(T_)) - E_
    );

    // Calculate radiative heat flux on boundaries.
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();
    const volScalarField::Boundary& GBf = G_.boundaryField();
    const volScalarField::Boundary& gammaBf = gamma.boundaryField();

    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (!GBf[patchi].coupled())
        {
            qrBf[patchi] = -gammaBf[patchi]*GBf[patchi].snGrad();
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::P1::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        IOobject::NO_REGISTER,
        4.0*absorptionEmission_->eCont()*physicoChemical::sigma
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::P1::Ru() const
{
    const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    const volScalarField::Internal a = absorptionEmission_->aCont()()();

    return a*G - E;
}


Foam::label Foam::radiation::P1::nBands() const
{
    return absorptionEmission_->nBands();
}


// ************************************************************************* //
