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

#include "exponential.H"
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
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable
    (
        nucleateFluxModel,
        exponential,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::nucleateFluxModels::exponential::exponential
(
    const dictionary& dict
)
:
    nucleateFluxModel(),
    a_(dict.getOrDefault<scalar>("a", 6309)),
    b_(dict.getOrDefault<scalar>("b", 2.52))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::nucleateFluxModels::exponential::qNucleate
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

    return a_*pow(max((Tw-Tsatw), scalar(0)), b_);
}


void Foam::wallBoilingModels::nucleateFluxModels::exponential::write
(
    Ostream& os
) const
{
    nucleateFluxModel::write(os);
    os.writeEntry("a", a_);
    os.writeEntry("b", b_);
}


// ************************************************************************* //
