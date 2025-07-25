/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "APIdiffCoefFunc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(APIdiffCoefFunc, 0);
    addToRunTimeSelectionTable
    (
        thermophysicalFunction,
        APIdiffCoefFunc,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::APIdiffCoefFunc::APIdiffCoefFunc
(
    const scalar a,
    const scalar b,
    const scalar wf,
    const scalar wa
)
:
    a_(a),
    b_(b),
    wf_(wf),
    wa_(wa),
    alpha_(sqrt(1/wf_ + 1/wa_)),
    beta_(sqr(cbrt(a_) + cbrt(b_)))
{}


Foam::APIdiffCoefFunc::APIdiffCoefFunc(const dictionary& dict)
:
    a_(dict.get<scalar>("a")),
    b_(dict.get<scalar>("b")),
    wf_(dict.get<scalar>("wf")),
    wa_(dict.get<scalar>("wa")),
    alpha_(sqrt(1/wf_ + 1/wa_)),
    beta_(sqr((cbrt(a_) + cbrt(b_))))
{}


// ************************************************************************* //
