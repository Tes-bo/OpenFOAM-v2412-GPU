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

#include "BreenWestwater.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phasePairKey.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{
    defineTypeNameAndDebug(BreenWestwater, 0);
    addToRunTimeSelectionTable
    (
        filmBoilingModel,
        BreenWestwater,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::BreenWestwater::BreenWestwater
(
    const dictionary& dict
)
:
    filmBoilingModel(),
    Cn_(dict.getOrDefault<scalar>("Cn", 0.37))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::filmBoilingModels::BreenWestwater::htcFilmBoil
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


    const scalarField kappaLiquid(liquid.kappa(patchi));

    tmp<scalarField> tCp = vapor.thermo().Cp(pw, Tsatw, cells);
    const scalarField& CpVapor = tCp();

    const scalarField nuLiquid(liquid.nu(patchi));

    const scalarField Leff
    (
        L*sqr(1 + 0.34*CpVapor*max((Tw-Tsatw), scalar(0))/L)
    );

    const scalarField rhoDiff(rhoLiq - rhoVapor);

    const phasePairKey pair(liquid.name(), vapor.name());
    const scalarField sigma
    (
        liquid.fluid().sigma(pair)().boundaryField()[patchi]
    );

    return
         Cn_
        /pow(sigma/mag(g.value())/rhoDiff, 1/8)
        /pow
         (
             nuLiquid*max((Tw-Tsatw), scalar(1e-3))
            /(pow3(kappaLiquid)*rhoVapor*Leff*mag(g.value())*rhoDiff)
            , 0.25
         );
}


void Foam::wallBoilingModels::filmBoilingModels::BreenWestwater::write
(
    Ostream& os
) const
{
    filmBoilingModel::write(os);
    os.writeEntry("Cn", Cn_);
}


// ************************************************************************* //
