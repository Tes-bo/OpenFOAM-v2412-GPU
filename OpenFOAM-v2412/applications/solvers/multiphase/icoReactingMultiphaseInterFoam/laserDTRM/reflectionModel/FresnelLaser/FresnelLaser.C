/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "FresnelLaser.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(FresnelLaser, 0);

        addToRunTimeSelectionTable
        (
            reflectionModel,
            FresnelLaser,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::FresnelLaser::FresnelLaser
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    reflectionModel(dict, mesh),
    epsilon_(dict.get<scalar>("epsilon"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::FresnelLaser::rho
(
    const scalar cosTheta
) const
{
    //scalar cosTheta = cos(incidentAngle);

    scalar rho =
        0.5
      * (
            (1 + sqr(1 - epsilon_*cosTheta))/(1 + sqr(1 + epsilon_*cosTheta))
        +
            (sqr(epsilon_) - 2*epsilon_*cosTheta + 2*sqr(cosTheta))
          /
            (sqr(epsilon_) + 2*epsilon_*cosTheta + 2*sqr(cosTheta))
        );

    return rho;
}


Foam::vector Foam::radiation::FresnelLaser::R
(
    const vector& i,
    const vector& n
) const
{
    return  i + 2.0*(-i & n) * n;
}


// ************************************************************************* //
