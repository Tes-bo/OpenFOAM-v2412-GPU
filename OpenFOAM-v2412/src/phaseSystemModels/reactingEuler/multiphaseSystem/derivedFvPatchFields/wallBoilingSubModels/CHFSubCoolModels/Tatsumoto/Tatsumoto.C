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

#include "Tatsumoto.H"
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
    defineTypeNameAndDebug(Tatsumoto, 0);
    addToRunTimeSelectionTable
    (
        CHFSubCoolModel,
        Tatsumoto,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::CHFModels::Tatsumoto::Tatsumoto
(
    const dictionary& dict
)
:
    CHFSubCoolModel(),
    K_(dict.getOrDefault<scalar>("K", 0.23))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::CHFModels::Tatsumoto::CHFSubCool
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    const labelUList& cells = liquid.mesh().boundary()[patchi].faceCells();
    const scalarField& pw = liquid.thermo().p().boundaryField()[patchi];

    tmp<scalarField> trhoVapor = vapor.thermo().rhoEoS(pw, Tsatw, cells);
    const scalarField& rhoVapor = trhoVapor.ref();

    tmp<scalarField> trhoLiq = liquid.thermo().rhoEoS(pw, Tsatw, cells);
    const scalarField& rhoLiq = trhoLiq.ref();

    tmp<scalarField> tCp = liquid.thermo().Cp(pw, Tsatw, cells);
    const scalarField& Cp = tCp();

    return
        1 + K_*pow(rhoVapor/rhoLiq, 0.8)*Cp*max(Tsatw - Tl, scalar(0))/L;
}


void Foam::wallBoilingModels::CHFModels::Tatsumoto::write
(
    Ostream& os
) const
{
    CHFSubCoolModel::write(os);
    os.writeEntry("K", K_);
}


// ************************************************************************* //
