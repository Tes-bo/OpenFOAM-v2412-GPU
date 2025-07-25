/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd
-------------------------------------------------------------------------------
License

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

#include "eddyDissipationModelBase.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModelBase<ReactionThermo, ThermoType>::eddyDissipationModelBase
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion<ReactionThermo, ThermoType>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    CEDC_(this->coeffs().getScalar("CEDC"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModelBase<ReactionThermo, ThermoType>::
~eddyDissipationModelBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModelBase<ReactionThermo, ThermoType>::rtTurb() const
{
    return
        CEDC_*this->turbulence().epsilon()
      / max
        (
            this->turbulence().k(),
            dimensionedScalar("SMALL",sqr(dimVelocity), SMALL)
        );
}


template<class ReactionThermo, class ThermoType>
void eddyDissipationModelBase<ReactionThermo, ThermoType>::correct()
{
    this->wFuel_ == dimensionedScalar(dimMass/dimVolume/dimTime, Zero);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel = this->thermo_.composition().Y()[fuelI];

        const dimensionedScalar s = this->singleMixturePtr_->s();

        if (this->thermo_.composition().contains("O2"))
        {
            const volScalarField& YO2 = this->thermo_.composition().Y("O2");

            this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                * timeScale();
        }
        else
        {
            FatalErrorInFunction
                << "You selected a combustion model that requires O2 mass"
                << " to be present in the mixture"
                << exit(FatalError);
        }
    }
}


template<class ReactionThermo, class ThermoType>
bool eddyDissipationModelBase<ReactionThermo, ThermoType>::read()
{
    if (singleStepCombustion<ReactionThermo, ThermoType>::read())
    {
        this->coeffs().readEntry("CEDC", CEDC_);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
