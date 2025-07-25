/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "atmBuoyancyTurbSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmBuoyancyTurbSource, 0);
    addToRunTimeSelectionTable(option, atmBuoyancyTurbSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fv::atmBuoyancyTurbSource::calcB()
{
    //- Temperature field [K]
    const volScalarField& T = mesh_.lookupObjectRef<volScalarField>("T");

    //- Kinematic turbulent thermal conductivity field [m2/s]
    const volScalarField& alphat =
        mesh_.lookupObjectRef<volScalarField>("alphat");

    // (ARAL:Eq. 7)
    B_ = beta_*alphat()*(fvc::grad(T) & g_)();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::atmBuoyancyTurbSource::calcC3
(
    const volScalarField::Internal& k,
    const volScalarField::Internal& epsilon,
    const volScalarField::Internal& G
) const
{
    // Gradient Richardson number (ARAL:p. 4)
    const volScalarField::Internal Rig
    (
        -B_/(G + dimensionedScalar(G.dimensions(), SMALL))
    );

    // Mixing-length scale estimation (P:Eq. 10.37 & p. 374) normalised by Lmax_
    const volScalarField::Internal LbyLmax
    (
        (pow(Cmu_, 0.75)/Lmax_)*pow(k, 1.5)/epsilon
    );

    // (ARAL:Eq. 10), with a typo of (C2_) instead of using (C2_ - 1.0)
    volScalarField::Internal alphaB(1.0 - LbyLmax);

    alphaB =
        neg0(Rig)*(1.0 - (1.0 + (C2_ - 1.0)/(C2_ - C1_))*LbyLmax)
      + pos(Rig)*(1.0 - LbyLmax);

    // (SKL:Eq. 18, rhs-term:3); (ARAL:Eq. 5, rhs-term:3) has a typo
    return (C1_ - C2_)*alphaB + 1.0;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::atmBuoyancyTurbSource::calcC3
(
    const volScalarField::Internal& k,
    const volScalarField::Internal& omega,
    const volScalarField::Internal& G,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    // Gradient Richardson number (ARAL:p. 4)
    const volScalarField::Internal Rig
    (
        -B_/(G + dimensionedScalar(G.dimensions(), SMALL))
    );

    // Mixing-length scale estimation (L:Eq. 3.20) normalised by Lmax_
    const volScalarField::Internal LbyLmax
    (
        (1.0/(pow025(Cmu_)*Lmax_))*sqrt(k)/omega
    );

    // (ARAL:Eq. 10)
    volScalarField::Internal alphaB(1.0 - LbyLmax);

    alphaB =
        neg0(Rig)*(1.0 - (1.0 + beta/(beta - gamma))*LbyLmax)
      + pos(Rig)*(1.0 - LbyLmax);

    // (SKL:Eq. 19, rhs-term:3); (ARAL:Eq. 5, rhs-term:3) has a typo
    return (gamma - beta)*alphaB;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmBuoyancyTurbSource::atmBuoyancyTurbSource
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
    Lmax_
    (
        dimensionedScalar
        (
            dimLength,
            coeffs_.getCheckOrDefault<scalar>
            (
                "Lmax",
                41.575,
                [&](const scalar Lmax){ return Lmax > SMALL; }
            )
        )
    ),
    beta_
    (
        dimensionedScalar
        (
            dimless/dimTemperature,
            coeffs_.getCheckOrDefault<scalar>
            (
                "beta",
                3.3e-3,
                [&](const scalar x){ return x > SMALL; }
            )
        )
    ),
    Cmu_(Zero),
    C1_(Zero),
    C2_(Zero),
    g_
    (
        "g",
        dimLength/sqr(dimTime),
        meshObjects::gravity::New(mesh_.time()).value()
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero)
    )
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.resize(2);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (tepsilon.is_reference())
    {
        isEpsilon_ = true;
        fieldNames_[0] = tepsilon().name();

        const dictionary& turbDict = turbPtr->coeffDict();

        Cmu_.read("Cmu", turbDict);
        C1_.read("C1", turbDict);
        C2_.read("C2", turbDict);
    }
    else if (tomega.is_reference())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();

        const dictionary& turbDict = turbPtr->coeffDict();

        Cmu_.read("betaStar", turbDict);
    }
    else
    {
        FatalErrorInFunction
            << "Needs either epsilon or omega field."
            << abort(FatalError);
    }

    fieldNames_[1] = turbPtr->k()().name();

    fv::option::resetApplied();

    Log << "    Applying atmBuoyancyTurbSource to: "
        << fieldNames_[0] << " and " << fieldNames_[1]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmBuoyancyTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmBuoyancyTurbSourceK
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
        return;
    }

    calcB();

    if (isEpsilon_)
    {
        atmBuoyancyTurbSourceEpsilon
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
    else
    {
        atmBuoyancyTurbSourceOmega
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
}


void Foam::fv::atmBuoyancyTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmBuoyancyTurbSourceK(geometricOneField(), rho, eqn, fieldi);
        return;
    }

    calcB();

    if (isEpsilon_)
    {
        atmBuoyancyTurbSourceEpsilon(geometricOneField(), rho, eqn, fieldi);
    }
    else
    {
        atmBuoyancyTurbSourceOmega(geometricOneField(), rho, eqn, fieldi);
    }
}


void Foam::fv::atmBuoyancyTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmBuoyancyTurbSourceK(alpha, rho, eqn, fieldi);
        return;
    }

    calcB();

    if (isEpsilon_)
    {
        atmBuoyancyTurbSourceEpsilon(alpha, rho, eqn, fieldi);
    }
    else
    {
        atmBuoyancyTurbSourceOmega(alpha, rho, eqn, fieldi);
    }
}


// ************************************************************************* //
