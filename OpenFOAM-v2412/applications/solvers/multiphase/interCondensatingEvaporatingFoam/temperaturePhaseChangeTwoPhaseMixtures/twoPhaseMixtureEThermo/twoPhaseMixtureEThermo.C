/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "twoPhaseMixtureEThermo.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "twoPhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseMixtureEThermo, 0);
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureEThermo::twoPhaseMixtureEThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    basicThermo(U.mesh(), word::null),
    thermoIncompressibleTwoPhaseMixture(U, phi),

    TSat_("TSat", dimTemperature, static_cast<const basicThermo&>(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseMixtureEThermo::correct()
{
    incompressibleTwoPhaseMixture::correct();
}


Foam::word Foam::twoPhaseMixtureEThermo::thermoName() const
{
    NotImplemented;
    return word::null;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::hc() const
{
    const fvMesh& mesh = this->T_.mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "hc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("hc", Hf2() - Hf1())
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::Cp() const
{
    const volScalarField limitedAlpha1
    (
        clamp(alpha1_, zero_one{})
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
            limitedAlpha1*Cp1() + (scalar(1) - limitedAlpha1)*Cp2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    return
    (
        alpha1p*Cp1().value() + (scalar(1) - alpha1p)*Cp2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::rho() const
{
    const volScalarField limitedAlpha1
    (
        clamp(alpha1_, zero_one{})
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*rho1().value()
          + (scalar(1) - limitedAlpha1)*rho2().value()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::rho
(
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    return
    (
        alpha1p*rho1().value() + (scalar(1) - alpha1p)*rho2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::Cv() const
{
    const volScalarField limitedAlpha1
    (
        clamp(alpha1_, zero_one{})
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cv",
            limitedAlpha1*Cv1() + (scalar(1) - limitedAlpha1)*Cv2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    return
    (
        alpha1p*Cv1().value() + (scalar(1) - alpha1p)*Cv2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::gamma() const
{
    return tmp<volScalarField>
    (
        (alpha1_*Cp1() + alpha2_*Cp2())/(alpha1_*Cv1() + alpha2_*Cv2())
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
    (
        gamma()().boundaryField()[patchi]
    );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::Cpv() const
{
     // This is an e thermo (Cpv = Cv)
     return Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    // This is an e thermo (Cpv = Cv)
    return Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::CpByCpv() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::W() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::kappa() const
{
    const volScalarField limitedAlpha1
    (
        clamp(alpha1_, zero_one{})
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "kappa",
            limitedAlpha1*kappa1() + (scalar(1) - limitedAlpha1)*kappa2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::kappa
(
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    return (alpha1p*kappa1().value() + (1 - alpha1p)*kappa2().value());
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::alphahe() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::alphahe
(
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::kappaEff
(
    const scalarField& kappat,
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    return
        (alpha1p*kappa1().value() + (1 - alpha1p)*kappa2().value()) + kappat;

}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    const volScalarField rho
    (
        alpha1_*rho1() + (1.0 - alpha1_)*rho2()
    );

    return (kappa()/Cp()/rho + alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureEThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    const scalarField alpha1p
    (
        clamp(alpha1_.boundaryField()[patchi], zero_one{})
    );

    const scalarField rho
    (
        alpha1p*rho1().value() + (1.0 - alpha1p)*rho2().value()
    );

    const scalarField kappa
    (
        alpha1p*kappa1().value() + (1.0 - alpha1p)*kappa2().value()
    );

    const scalarField Cp
    (
        alpha1p*Cp1().value() + (1.0 - alpha1p)*Cp2().value()
    );

    return kappa/Cp/rho + alphat;
}


bool Foam::twoPhaseMixtureEThermo::read()
{
    if (basicThermo::read() && thermoIncompressibleTwoPhaseMixture::read())
    {
        basicThermo::readEntry("TSat", TSat_);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
