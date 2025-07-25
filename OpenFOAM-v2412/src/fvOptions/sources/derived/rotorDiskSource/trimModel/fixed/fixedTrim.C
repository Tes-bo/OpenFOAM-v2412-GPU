/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "fixedTrim.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedTrim, 0);
    addToRunTimeSelectionTable(trimModel, fixedTrim, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTrim::fixedTrim
(
    const fv::rotorDiskSource& rotor,
    const dictionary& dict
)
:
    trimModel(rotor, dict, typeName),
    thetag_(rotor.cells().size(), Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedTrim::read(const dictionary& dict)
{
    trimModel::read(dict);

    const scalar theta0 = degToRad(coeffs_.get<scalar>("theta0"));
    const scalar theta1c = degToRad(coeffs_.get<scalar>("theta1c"));
    const scalar theta1s = degToRad(coeffs_.get<scalar>("theta1s"));

    const List<point>& x = rotor_.x();
    forAll(thetag_, i)
    {
        scalar psi = x[i].y();
        thetag_[i] = theta0 + theta1c*cos(psi) + theta1s*sin(psi);
    }
}


Foam::tmp<Foam::scalarField> Foam::fixedTrim::thetag() const
{
    return tmp<scalarField>(thetag_);
}


void Foam::fixedTrim::correct
(
    const vectorField& U,
    vectorField& force
)
{}


void Foam::fixedTrim::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force
)
{}


// ************************************************************************* //
