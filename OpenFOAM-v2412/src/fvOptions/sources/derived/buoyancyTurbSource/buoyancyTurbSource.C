/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "buoyancyTurbSource.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(buoyancyTurbSource, 0);
    addToRunTimeSelectionTable(option, buoyancyTurbSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal> Foam::fv::buoyancyTurbSource::
B() const
{
    const auto& alphat = mesh_.lookupObjectRef<volScalarField>(alphatName_);
    const auto& T = mesh_.lookupObjectRef<volScalarField>(Tname_);

    // (BMA:Eq. 8)
    return beta_*alphat()*(fvc::grad(T) & g_)();
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceEpsilon
(
    fvMatrix<scalar>& eqn
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const dictionary& turbDict = turbPtr->coeffDict();
    const dimensionedScalar C1(turbDict.getOrDefault<scalar>("C1", 1.44));
    const dimensionedScalar Cmu(turbDict.getOrDefault<scalar>("Cmu", 0.09));

    const volScalarField& epsilon = eqn.psi();
    const volScalarField::Internal& k = turbPtr->k()();
    const volVectorField::Internal& U = turbPtr->U()();
    const dimensionedScalar k0(k.dimensions(), SMALL);

    // (BMA:Eq. 9)
    const vector gHat(g_.value()/mag(g_.value()));

    volScalarField::Internal v(gHat & U);
    volScalarField::Internal u
    (
        mag(U - gHat*v)
      + dimensionedScalar(dimVelocity, SMALL)
    );

    // (BMA:Eq. 6)
    eqn -= fvm::SuSp(C1*tanh(mag(u/v))*B()/(k + k0), epsilon);
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega
(
    fvMatrix<scalar>& eqn
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    const volScalarField::Internal& nut = turbPtr->nut()();
    const auto& gamma =
        mesh_.lookupObjectRef<volScalarField::Internal>
        (
            IOobject::scopedName(turbPtr->type(), "gamma")
        );

    // (Heuristic approximation to BMA:Eq. 6)
    eqn -= gamma/(nut + dimensionedScalar(nut.dimensions(), SMALL))*B();
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK
(
    fvMatrix<scalar>& eqn
) const
{
    const volScalarField& k = eqn.psi();
    const dimensionedScalar k0(k.dimensions(), SMALL);

    // (BMA: Eq. 5)
    eqn -= fvm::Sp(B()/(k() + k0), k);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::buoyancyTurbSource::buoyancyTurbSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    isEpsilon_(false),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    alphatName_(coeffs_.getOrDefault<word>("alphat", "alphat")),
    Tname_(coeffs_.getOrDefault<word>("T", "T")),
    beta_
    (
        dimensionedScalar
        (
            dimless/dimTemperature,
            coeffs_.getCheckOrDefault<scalar>
            (
                "beta",
                3.3e-3,
                [=](const scalar x){ return x > SMALL; }
            )
        )
    ),
    g_
    (
        dimLength/sqr(dimTime),
        meshObjects::gravity::New(mesh_.time()).value()
    )
{
    if (mag(g_.value()) < SMALL)
    {
        FatalErrorInFunction
            << "Gravitational field cannot be equal to or less than zero"
            << exit(FatalError);
    }

    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << exit(FatalError);
    }

    fieldNames_.resize(2);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (tepsilon.is_reference())
    {
        isEpsilon_ = true;
        fieldNames_[0] = tepsilon().name();
    }
    else if (tomega.is_reference())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find an omega or epsilon field." << nl
            << "buoyancyTurbSource needs an omega- or epsilon-based model."
            << exit(FatalError);
    }

    fieldNames_[1] = turbPtr->k()().name();

    fv::option::resetApplied();

    Log << "    Applying buoyancyTurbSource to: "
        << fieldNames_[0] << " and " << fieldNames_[1]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(eqn);
        return;
    }

    if (isEpsilon_)
    {
        buoyancyTurbSourceEpsilon(eqn);
    }
    else
    {
        buoyancyTurbSourceOmega(eqn);
    }
}


void Foam::fv::buoyancyTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(geometricOneField(), rho, eqn, fieldi);
        return;
    }
}


void Foam::fv::buoyancyTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(alpha, rho, eqn, fieldi);
        return;
    }
}


// ************************************************************************* //
