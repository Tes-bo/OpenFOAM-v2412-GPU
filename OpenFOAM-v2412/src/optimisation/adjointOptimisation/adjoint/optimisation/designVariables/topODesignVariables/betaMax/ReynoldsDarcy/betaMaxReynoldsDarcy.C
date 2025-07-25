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

#include "betaMaxReynoldsDarcy.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(betaMaxReynoldsDarcy, 0);
    addToRunTimeSelectionTable(betaMax, betaMaxReynoldsDarcy, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::betaMaxReynoldsDarcy::betaMaxReynoldsDarcy
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    betaMax(mesh, dict),
    ReynoldsDarcyNumber_
    (
        dict.subDict(type() + "Coeffs").getOrDefault<scalar>("ReDa", 1.e-5)
    ),
    length_(computeLength(dict)),
    Uref_(dict.subDict(type() + "Coeffs").get<scalar>("Uref"))
{
    value_ = Uref_/ReynoldsDarcyNumber_/length_;
    Info
        << "Computed a betaMax value of " << value_
        << " based on a length of " << length_ << nl << endl;
}


// ************************************************************************* //
