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

#include "Shirai.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{
    defineTypeNameAndDebug(Shirai, 0);
    addToRunTimeSelectionTable
    (
        TDNBModel,
        Shirai,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::TDNBModels::Shirai::Shirai
(
    const dictionary& dict
)
:
    TDNBModel(),
    Tc_(dict.get<scalar>("Tc")),
    Pc_(dict.get<scalar>("Pc"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::TDNBModels::Shirai::TDNB
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    tmp<scalarField> tp = liquid.thermo().p().boundaryField()[patchi];

    const scalarField pRatio(clamp(tp/Pc_, zero_one{}));

    return
    (
        (0.8823*pow3(pRatio) - 1.8938*sqr(pRatio) + 1.4322*pRatio + 0.6289)*Tc_
    );
}


void Foam::wallBoilingModels::TDNBModels::Shirai::write
(
    Ostream& os
) const
{
    TDNBModel::write(os);
    os.writeEntry("Tc", Tc_);
    os.writeEntry("Pc", Pc_);
}


// ************************************************************************* //
