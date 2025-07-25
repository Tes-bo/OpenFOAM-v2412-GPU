/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "Moraga.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(Moraga, 0);
    addToRunTimeSelectionTable(liftModel, Moraga, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::Moraga::Moraga
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::Moraga::~Moraga()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::Moraga::Cl() const
{
    volScalarField Re(pair_.Re());

    volScalarField sqrSr
    (
        sqr(pair_.dispersed().d())
       /pair_.continuous().nu()
       *mag(fvc::grad(pair_.continuous().U()))
    );

    if
    (
        min(Re).value() < 1200.0
     || max(Re).value() > 18800.0
     || min(sqrSr).value() < 0.0016
     || max(sqrSr).value() > 0.04
    )
    {
        WarningInFunction
            << "Re and/or Sr are out of the range of applicability of the "
            << "Moraga model. Clamping to range bounds"
            << endl;
    }

    Re.clamp_range(1200.0, 18800.0);
    sqrSr.clamp_range(0.0016, 0.04);

    return 0.2*exp(- Re*sqrSr/3.6e5 - 0.12)*exp(Re*sqrSr/3.0e7);
}


// ************************************************************************* //
