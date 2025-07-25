/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd
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

#include "Kutadeladze.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phasePairKey.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleateFluxModels
{
    defineTypeNameAndDebug(Kutadeladze, 0);
    addToRunTimeSelectionTable
    (
        nucleateFluxModel,
        Kutadeladze,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleateFluxModels::Kutadeladze::Kutadeladze
(
    const dictionary& dict
)
:
    nucleateFluxModel(),
    Cn_(dict.getOrDefault<scalar>("Cn", 5.66e-10)),
    an_(dict.getOrDefault<scalar>("an", 2.5)),
    bn_(dict.getOrDefault<scalar>("bn", 1)),
    n_(dict.getOrDefault<scalar>("n", 1))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::nucleateFluxModels::Kutadeladze::qNucleate
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    const auto& p = liquid.mesh().lookupObject<volScalarField>("p");
    const scalarField& pb = p.boundaryField()[patchi];

    const labelUList& cells = liquid.mesh().boundary()[patchi].faceCells();

    tmp<scalarField> trhoVapor = vapor.thermo().rhoEoS(pb, Tsatw, cells);
    const scalarField& rhoVapor = trhoVapor.ref();

    tmp<scalarField> trhoLiq = liquid.thermo().rhoEoS(pb, Tsatw, cells);
    const scalarField& rhoLiq = trhoLiq.ref();

    const phasePairKey pair(liquid.name(), vapor.name());
    const scalarField sigma
    (
        liquid.fluid().sigma(pair)().boundaryField()[patchi]
    );

    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];

    const scalarField kappaLiquid(liquid.kappa(patchi));

    tmp<scalarField> tCpliq = liquid.thermo().Cp(pb, Tsatw, cells);
    const scalarField& Cpliquid = tCpliq();

    const scalarField muLiquid(liquid.mu(patchi));

    const scalarField deltaTsub
    (
        pow(max((Tw-Tsatw), scalar(0)), an_)
    );

    return
        Cn_*kappaLiquid*pow(Cpliquid, 1.5)*pow(rhoLiq,1.28)*pow(pb,1.75)
       *deltaTsub
       /
       (pow(muLiquid, 0.625)*pow(sigma,0.9)*pow(L,1.5)*pow(rhoVapor,1.5));
}


void Foam::wallBoilingModels::nucleateFluxModels::Kutadeladze::write
(
    Ostream& os
) const
{
    nucleateFluxModel::write(os);
    os.writeEntry("Cn", Cn_);
    os.writeEntry("an", an_);
    os.writeEntry("bn", bn_);
    os.writeEntry("n", n_);
}


// ************************************************************************* //
