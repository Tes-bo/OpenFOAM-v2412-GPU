/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "betaMaxStepRamp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(betaMaxStepRamp, 0);
    addToRunTimeSelectionTable
    (
       betaMax,
       betaMaxStepRamp,
       dictionary
    );
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::betaMaxStepRamp::betaMaxStepRamp
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    betaMax(mesh, dict),
    betaMin_(dict.getOrDefault<scalar>("betaMin", Zero)),
    funcPtr_(nullptr)
{
    funcPtr_.reset
        (Function1<scalar>::New("betaMaxStepRamp", dict, "stepRamp"));

    DebugInfo
        << "betaMaxStepRamp:: creating interpolation function of type "
        << funcPtr_().type() << endl;
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

Foam::scalar Foam::betaMaxStepRamp::value() const
{
    const scalar t = mesh_.time().timeOutputValue();
    const scalar value = betaMin_ + (value_ - betaMin_)*funcPtr_().value(t);
    DebugInfo
        << "stepRamp betaMax:: t, betaMax value "
        << t << ", " << value << endl;

    return value;
}

// ************************************************************************* //
