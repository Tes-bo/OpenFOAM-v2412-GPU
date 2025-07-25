/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd
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

#include "Bromley.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{
    defineTypeNameAndDebug(Bromley, 0);
    addToRunTimeSelectionTable
    (
        filmBoilingModel,
        Bromley,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::Bromley::Bromley
(
    const dictionary& dict
)
:
    filmBoilingModel(),
    Cn_(dict.getOrDefault<scalar>("Cn", 0.62)),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 1)),
    L_(dict.get<scalar>("L"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::filmBoilingModels::Bromley::htcFilmBoil
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];
    const auto& g =
        liquid.mesh().time().lookupObject<uniformDimensionedVectorField>("g");

    const labelUList& cells = liquid.mesh().boundary()[patchi].faceCells();

    const scalarField& pw = liquid.thermo().p().boundaryField()[patchi];

    tmp<scalarField> trhoVapor = vapor.thermo().rhoEoS(pw, Tsatw, cells);
    const scalarField& rhoVapor = trhoVapor.ref();

    tmp<scalarField> trhoLiq = liquid.thermo().rhoEoS(pw, Tsatw, cells);
    const scalarField& rhoLiq = trhoLiq.ref();


    const scalarField kappaVapor(vapor.kappa(patchi));

    tmp<scalarField> tCp = vapor.thermo().Cp(pw, Tsatw, cells);
    const scalarField& CpVapor = tCp();

    const scalarField muVapor(vapor.mu(patchi));

    const scalarField htcRad
    (
        emissivity_*physicoChemical::sigma.value()*(pow4(Tw) - pow4(Tsatw))
      / max((Tw - Tsatw), scalar(1e-4))
    );

    return
        Cn_*pow025
        (
            pow3(kappaVapor)
           *rhoVapor*(rhoLiq - rhoVapor)*mag(g.value())
           *(L + scalar(0.4)*CpVapor*max((Tw-Tsatw), scalar(0)))
           /(L_*muVapor*max((Tw-Tsatw), scalar(1e-4)))
        ) + scalar(0.75)*htcRad;
}


void Foam::wallBoilingModels::filmBoilingModels::Bromley::write
(
    Ostream& os
) const
{
    filmBoilingModel::write(os);
    os.writeEntry("Cn", Cn_);
    os.writeEntry("L", L_);
    os.writeEntry("emissivity", emissivity_);
}


// ************************************************************************* //
