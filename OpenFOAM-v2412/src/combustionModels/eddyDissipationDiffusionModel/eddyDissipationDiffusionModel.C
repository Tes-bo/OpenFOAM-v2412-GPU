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

#include "eddyDissipationDiffusionModel.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationDiffusionModel<ReactionThermo, ThermoType>::
eddyDissipationDiffusionModel
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    eddyDissipationModelBase<ReactionThermo, ThermoType>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    Cd_(this->coeffs().getScalar("Cd"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationDiffusionModel<ReactionThermo, ThermoType>::
~eddyDissipationDiffusionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionModel<ReactionThermo, ThermoType>::timeScale()
{
    return (max(this->rtTurb(), this->rtDiff()));
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationDiffusionModel<ReactionThermo, ThermoType>::rtDiff() const
{
    auto tdelta = volScalarField::New
    (
        "tdelta",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimLength, Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    auto& delta = tdelta.ref();

    delta.ref() = cbrt(this->mesh().V());
    delta.correctBoundaryConditions();

    // NOTE: Assume Prt = 1
    return Cd_*this->turbulence().nuEff()/sqr(delta);
}


template<class ReactionThermo, class ThermoType>
bool eddyDissipationDiffusionModel<ReactionThermo, ThermoType>::read()
{
    if (eddyDissipationModelBase<ReactionThermo, ThermoType>::read())
    {
        this->coeffs().readEntry("Cd", Cd_);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
