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

#include "Zuber.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phasePairKey.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{
    defineTypeNameAndDebug(Zuber, 0);
    addToRunTimeSelectionTable
    (
        CHFModel,
        Zuber,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::Zuber::Zuber
(
    const dictionary& dict
)
:
    CHFModel(),
    Cn_(dict.getOrDefault<scalar>("Cn", 0.131))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::CHFModels::Zuber::CHF
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    const auto& g =
        liquid.mesh().time().lookupObject<uniformDimensionedVectorField>("g");

    const labelUList& cells = liquid.mesh().boundary()[patchi].faceCells();

    const scalarField& pw = liquid.thermo().p().boundaryField()[patchi];

    tmp<scalarField> trhoVapor = vapor.thermo().rhoEoS(pw ,Tsatw, cells);
    const scalarField& rhoVapor = trhoVapor.ref();

    tmp<scalarField> trhoLiq = liquid.thermo().rhoEoS(pw, Tsatw, cells);
    const scalarField& rhoLiq = trhoLiq.ref();
    const phasePairKey pair(liquid.name(), vapor.name());
    const scalarField sigma
    (
        liquid.fluid().sigma(pair)().boundaryField()[patchi]
    );

    return
        Cn_*rhoVapor*L*pow
        (
            sigma*mag(g.value())*(max(rhoLiq - rhoVapor, Zero)/sqr(rhoVapor)),
            0.25
        );
}


void Foam::wallBoilingModels::CHFModels::Zuber::write(Ostream& os) const
{
    CHFModel::write(os);
    os.writeEntry("Cn", Cn_);
}


// ************************************************************************* //
