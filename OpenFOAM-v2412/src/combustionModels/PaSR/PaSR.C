/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "PaSR.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::PaSR<ReactionThermo>::PaSR
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    Cmix_(this->coeffs().getScalar("Cmix")),
    kappa_
    (
        IOobject
        (
            thermo.phaseScopedName(typeName, "kappa"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::PaSR<ReactionThermo>::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::PaSR<ReactionThermo>::correct()
{
    if (this->active())
    {
        laminar<ReactionThermo>::correct();

        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const auto& epsilon = tepsilon();

        tmp<volScalarField> tmuEff(this->turbulence().muEff());
        const auto& muEff = tmuEff();

        tmp<volScalarField> ttc(this->tc());
        const auto& tc = ttc();

        tmp<volScalarField> trho(this->rho());
        const auto& rho = trho();

        forAll(epsilon, i)
        {
            const scalar tk =
                Cmix_*sqrt(max(muEff[i]/rho[i]/(epsilon[i] + SMALL), 0));

            if (tk > SMALL)
            {
                kappa_[i] = tc[i]/(tc[i] + tk);
            }
            else
            {
                kappa_[i] = 1.0;
            }
        }

        // Evaluate bcs
        kappa_.correctBoundaryConditions();
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR<ReactionThermo>::R(volScalarField& Y) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<ReactionThermo>::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phaseScopedName(typeName, "Qdot"),
        IOobject::NO_REGISTER,
        kappa_*laminar<ReactionThermo>::Qdot()
    );
}


template<class ReactionThermo>
bool Foam::combustionModels::PaSR<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        this->coeffs().readEntry("Cmix", Cmix_);
        return true;
    }

    return false;
}


// ************************************************************************* //
