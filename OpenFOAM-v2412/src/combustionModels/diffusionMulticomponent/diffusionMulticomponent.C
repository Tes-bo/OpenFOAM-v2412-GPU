/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "diffusionMulticomponent.H"
#include "fvcGrad.H"
#include "reactingMixture.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::combustionModels::
diffusionMulticomponent<ReactionThermo, ThermoType>::init()
{
    // Load default values
    this->coeffs().readIfPresent("Ci", Ci_);
    this->coeffs().readIfPresent("YoxStream", YoxStream_);
    this->coeffs().readIfPresent("YfStream", YfStream_);
    this->coeffs().readIfPresent("sigma", sigma_);
    this->coeffs().readIfPresent("ftCorr", ftCorr_);
    this->coeffs().readIfPresent("alpha", alpha_);
    this->coeffs().readIfPresent("laminarIgn", laminarIgn_);

    typedef typename Reaction<ThermoType>::specieCoeffs specieCoeffs;

    const speciesTable& species = this->thermo().composition().species();

    scalarList specieStoichCoeffs(species.size());
    const label nReactions = reactions_.size();

    for (label k=0; k < nReactions; k++)
    {
        RijPtr_.set
        (
            k,
            new volScalarField
            (
                this->mesh_.newIOobject("Rijk" + Foam::name(k)),
                this->mesh_,
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero),
                fvPatchFieldBase::zeroGradientType()
            )
        );

        RijPtr_[k].storePrevIter();

        const List<specieCoeffs>& lhs = reactions_[k].lhs();
        const List<specieCoeffs>& rhs = reactions_[k].rhs();

        const label fuelIndex = species.find(fuelNames_[k]);
        const label oxidantIndex = species.find(oxidantNames_[k]);

        const scalar Wu = specieThermo_[fuelIndex].W();
        const scalar Wox = specieThermo_[oxidantIndex].W();

        forAll(lhs, i)
        {
            const label specieI = lhs[i].index;
            specieStoichCoeffs[specieI] = -lhs[i].stoichCoeff;
            qFuel_[k] +=
                specieThermo_[specieI].hc()*lhs[i].stoichCoeff/Wu;
        }

        forAll(rhs, i)
        {
            const label specieI = rhs[i].index;
            specieStoichCoeffs[specieI] = rhs[i].stoichCoeff;
            qFuel_[k] -=
                specieThermo_[specieI].hc()*rhs[i].stoichCoeff/Wu;
        }

        Info << "Fuel heat of combustion : " << qFuel_[k] << endl;

        s_[k] =
            (Wox*mag(specieStoichCoeffs[oxidantIndex]))
          / (Wu*mag(specieStoichCoeffs[fuelIndex]));

        Info << "stoichiometric oxygen-fuel ratio : " << s_[k] << endl;

        stoicRatio_[k] = s_[k]*YfStream_[k]/YoxStream_[k];

        Info << "stoichiometric air-fuel ratio : " << stoicRatio_[k] << endl;

        const scalar fStoich = 1.0/(1.0 + stoicRatio_[k]);

        Info << "stoichiometric mixture fraction : " << fStoich << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::combustionModels::diffusionMulticomponent<ReactionThermo, ThermoType>::
diffusionMulticomponent
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ChemistryCombustion<ReactionThermo>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(thermo)
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(thermo).
            speciesData()
    ),
    RijPtr_(reactions_.size()),
    Ci_(reactions_.size(), 1.0),
    fuelNames_(this->coeffs().lookup("fuels")),
    oxidantNames_(this->coeffs().lookup("oxidants")),
    qFuel_(reactions_.size()),
    stoicRatio_(reactions_.size()),
    s_(reactions_.size()),
    YoxStream_(reactions_.size(), 0.23),
    YfStream_(reactions_.size(), 1.0),
    sigma_(reactions_.size(), 0.02),
    oxidantRes_(this->coeffs().lookup("oxidantRes")),
    ftCorr_(reactions_.size(), Zero),
    alpha_(1),
    laminarIgn_(false)
{
    init();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField> Foam::combustionModels::
diffusionMulticomponent<ReactionThermo, ThermoType>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class ReactionThermo, class ThermoType>
void Foam::combustionModels::
diffusionMulticomponent<ReactionThermo, ThermoType>::correct()
{
    if (this->active())
    {
        typedef typename Reaction<ThermoType>::specieCoeffs specieCoeffs;
        const speciesTable& species = this->thermo().composition().species();

        const label nReactions = reactions_.size();

        PtrList<volScalarField> RijlPtr(nReactions);

        for (label k=0; k < nReactions; k++)
        {
            RijlPtr.set
            (
                k,
                new volScalarField
                (
                    this->mesh_.newIOobject("Rijl" + Foam::name(k)),
                    this->mesh_,
                    dimensionedScalar(dimMass/dimTime/dimVolume, Zero),
                    fvPatchFieldBase::zeroGradientType()
                )
            );

            volScalarField& Rijl = RijlPtr[k];

            // Obtain individual reaction rates for each reaction
            const label fuelIndex = species.find(fuelNames_[k]);

            if (laminarIgn_)
            {
                Rijl.ref() = -this->chemistryPtr_->calculateRR(k, fuelIndex);
            }


             // Look for the fuelStoic
            const List<specieCoeffs>& rhs = reactions_[k].rhs();
            const List<specieCoeffs>& lhs = reactions_[k].lhs();

            // Set to zero RR's
            forAll(lhs, l)
            {
                const label lIndex = lhs[l].index;
                this->chemistryPtr_->RR(lIndex) =
                    dimensionedScalar(dimMass/dimTime/dimVolume, Zero);
            }

            forAll(rhs, l)
            {
                const label rIndex = rhs[l].index;
                this->chemistryPtr_->RR(rIndex) =
                    dimensionedScalar(dimMass/dimTime/dimVolume, Zero);
            }
        }


        for (label k=0; k < nReactions; k++)
        {
            const label fuelIndex = species.find(fuelNames_[k]);
            const label oxidantIndex = species.find(oxidantNames_[k]);

            const volScalarField& Yfuel =
                this->thermo().composition().Y(fuelIndex);

            const volScalarField& Yox =
                this->thermo().composition().Y(oxidantIndex);

            const volScalarField ft
            (
                "ft" + Foam::name(k),
                (
                    s_[k]*Yfuel - (Yox - YoxStream_[k])
                )
               /(
                    s_[k]*YfStream_[k] + YoxStream_[k]
                )
            );

            const scalar sigma = sigma_[k];

            const scalar fStoich = 1.0/(1.0 + stoicRatio_[k]) + ftCorr_[k];

            const volScalarField OAvailScaled
            (
                "OAvailScaled",
                Yox/max(oxidantRes_[k], 1e-3)
            );

            const volScalarField preExp
            (
                "preExp" + Foam::name(k),
                 1.0  + sqr(OAvailScaled)
            );

            const volScalarField filter
            (
                (1.0/(sigma*sqrt(2.0*constant::mathematical::pi)))
              * exp(-sqr(ft - fStoich)/(2*sqr(sigma)))
            );

            const volScalarField topHatFilter(pos(filter - 1e-3));

            const volScalarField prob("prob" + Foam::name(k), preExp*filter);

            const volScalarField RijkDiff
            (
               "RijkDiff",
                Ci_[k]*this->turbulence().muEff()*prob*
                (
                    mag(fvc::grad(Yfuel) & fvc::grad(Yox))
                )
               *pos(Yox)*pos(Yfuel)
            );

            volScalarField& Rijk = RijPtr_[k];

            if (laminarIgn_)
            {
                Rijk =
                    min(RijkDiff, topHatFilter*RijlPtr[k]*pos(Yox)*pos(Yfuel));
            }
            else
            {
                Rijk = RijkDiff;
            }

            Rijk.relax(alpha_);

            if (debug && this->mesh_.time().writeTime())
            {
                Rijk.write();
                ft.write();
            }

            // Look for the fuelStoic
            const List<specieCoeffs>& rhs = reactions_[k].rhs();
            const List<specieCoeffs>& lhs = reactions_[k].lhs();

            scalar fuelStoic = 1.0;
            forAll(lhs, l)
            {
                const label lIndex = lhs[l].index;
                if (lIndex == fuelIndex)
                {
                    fuelStoic = lhs[l].stoichCoeff;
                    break;
                }
            }

            const scalar MwFuel = specieThermo_[fuelIndex].W();

            // Update left hand side species
            forAll(lhs, l)
            {
                const label lIndex = lhs[l].index;

                const scalar stoichCoeff = lhs[l].stoichCoeff;

                this->chemistryPtr_->RR(lIndex) +=
                   -Rijk*stoichCoeff*specieThermo_[lIndex].W()/fuelStoic/MwFuel;

            }

            // Update right hand side species
            forAll(rhs, r)
            {
                const label rIndex = rhs[r].index;

                const scalar stoichCoeff = rhs[r].stoichCoeff;

                this->chemistryPtr_->RR(rIndex) +=
                    Rijk*stoichCoeff*specieThermo_[rIndex].W()/fuelStoic/MwFuel;
            }
        }
    }
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::diffusionMulticomponent<ReactionThermo, ThermoType>::R
(
    volScalarField& Y
) const
{
    auto tSu = tmp<fvScalarMatrix>::New(Y, dimMass/dimTime);
    auto& Su = tSu.ref();

    if (this->active())
    {
        const label specieI =
            this->thermo().composition().species().find(Y.member());

        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::
diffusionMulticomponent<ReactionThermo, ThermoType>::Qdot() const
{
    auto tQdot = volScalarField::New
    (
        "Qdot",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, Zero),
        fvPatchFieldBase::zeroGradientType()
    );

    if (this->active())
    {
        volScalarField& Qdot = tQdot.ref();
        Qdot = this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
bool Foam::combustionModels::
diffusionMulticomponent<ReactionThermo, ThermoType>::read()
{
    if (ChemistryCombustion<ReactionThermo>::read())
    {
        this->coeffs().readIfPresent("Ci", Ci_);
        this->coeffs().readIfPresent("sigma", sigma_);
        this->coeffs().readIfPresent("oxidantRes", oxidantRes_);
        this->coeffs().readIfPresent("ftCorr", ftCorr_);
        this->coeffs().readIfPresent("alpha", alpha_);
        this->coeffs().readIfPresent("laminarIgn", laminarIgn_);
        return true;
    }

    return false;
}


// ************************************************************************* //
