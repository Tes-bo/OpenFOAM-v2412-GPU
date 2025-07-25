/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "EDC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::EDC<ReactionThermo>::EDC
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    version_
    (
        EDCversionNames.getOrDefault
        (
            "version",
            this->coeffs(),
            EDCdefaultVersion
        )
    ),
    C1_(this->coeffs().getOrDefault("C1", 0.05774)),
    C2_(this->coeffs().getOrDefault("C2", 0.5)),
    Cgamma_(this->coeffs().getOrDefault("Cgamma", 2.1377)),
    Ctau_(this->coeffs().getOrDefault("Ctau", 0.4083)),
    exp1_(this->coeffs().getOrDefault("exp1", EDCexp1[int(version_)])),
    exp2_(this->coeffs().getOrDefault("exp2", EDCexp2[int(version_)])),
    kappa_
    (
        IOobject
        (
            this->thermo().phaseScopedName(typeName, "kappa"),
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
Foam::combustionModels::EDC<ReactionThermo>::~EDC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::EDC<ReactionThermo>::correct()
{
    if (this->active())
    {
        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const auto& epsilon = tepsilon();

        tmp<volScalarField> tmu(this->turbulence().mu());
        const auto& mu = tmu();

        tmp<volScalarField> tk(this->turbulence().k());
        const auto& k = tk();

        tmp<volScalarField> trho(this->rho());
        const auto& rho = trho();

        scalarField tauStar(epsilon.size(), Zero);

        if (version_ == EDCversions::v2016)
        {
            tmp<volScalarField> ttc(this->chemistryPtr_->tc());
            const auto& tc = ttc();

            forAll(tauStar, i)
            {
                const scalar nu = mu[i]/(rho[i] + SMALL);

                const scalar Da = clamp
                (
                    sqrt(nu/(epsilon[i] + SMALL))/tc[i],
                    scalar(1e-10), scalar(10)
                );

                const scalar ReT = sqr(k[i])/(nu*epsilon[i] + SMALL);
                const scalar CtauI = min(C1_/(Da*sqrt(ReT + 1)), 2.1377);

                const scalar CgammaI =
                    clamp(C2_*sqrt(Da*(ReT + 1)), scalar(0.4082), scalar(5));

                const scalar gammaL =
                    CgammaI*pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));

                tauStar[i] = CtauI*sqrt(nu/(epsilon[i] + SMALL));

                if (gammaL >= 1)
                {
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }

            // Evaluate bcs
            kappa_.correctBoundaryConditions();
        }
        else
        {
            forAll(tauStar, i)
            {
                const scalar nu = mu[i]/(rho[i] + SMALL);
                const scalar gammaL =
                    Cgamma_*pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));

                tauStar[i] = Ctau_*sqrt(nu/(epsilon[i] + SMALL));
                if (gammaL >= 1)
                {
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }

            // Evaluate bcs
            kappa_.correctBoundaryConditions();
        }

        Info<< "Chemistry time solved max/min : "
            << gMax(tauStar) << " / " << gMin(tauStar) << endl;

        this->chemistryPtr_->solve(tauStar);
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::EDC<ReactionThermo>::R(volScalarField& Y) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC<ReactionThermo>::Qdot() const
{
    auto tQdot = volScalarField::New
    (
        this->thermo().phaseScopedName(typeName, "Qdot"),
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    if (this->active())
    {
        tQdot.ref() = kappa_*this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class ReactionThermo>
bool Foam::combustionModels::EDC<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        version_ =
        (
            EDCversionNames.getOrDefault
            (
                "version",
                this->coeffs(),
                EDCdefaultVersion
            )
        );
        C1_ = this->coeffs().getOrDefault("C1", 0.05774);
        C2_ = this->coeffs().getOrDefault("C2", 0.5);
        Cgamma_ = this->coeffs().getOrDefault("Cgamma", 2.1377);
        Ctau_ = this->coeffs().getOrDefault("Ctau", 0.4083);
        exp1_ = this->coeffs().getOrDefault("exp1", EDCexp1[int(version_)]);
        exp2_ = this->coeffs().getOrDefault("exp2", EDCexp2[int(version_)]);

        return true;
    }

    return false;
}


// ************************************************************************* //
