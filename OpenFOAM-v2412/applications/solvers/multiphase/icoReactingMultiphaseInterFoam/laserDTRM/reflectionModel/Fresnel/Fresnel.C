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

#include "Fresnel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(Fresnel, 0);

        addToRunTimeSelectionTable
        (
            reflectionModel,
            Fresnel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::Fresnel::Fresnel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    reflectionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    nk1_(coeffsDict_.lookup("nk1")),
    nk2_(coeffsDict_.lookup("nk2"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::Fresnel::rho
(
    const scalar incidentAngle
) const
{
    // absorbing madium
    scalar n1 = sqr(nk1_[0]);
    //scalar k1 = sqr(nk1_[1]);

    // dialectric
    scalar n2 = sqr(nk2_[0]);
    scalar k2 = sqr(nk2_[1]);

    scalar sinTheta1 = sin(incidentAngle);

    scalar sqrP =
        0.5*
        (
            sqrt
            (
                sqr(n2-k2-n1*sqr(sinTheta1)) + 4*n2*k2
            )
          + (n2-k2-n1*sqr(sinTheta1))
        );

   scalar sqrQ =
        0.5*
        (
            sqrt
            (
                sqr(n2-k2-n1*sqr(sinTheta1)) + 4*n2*k2
            )
          - (n2-k2-n1*sqr(sinTheta1))
        );

    scalar cosTheta1 = cos(incidentAngle);
    scalar tanTheta1 = tan(incidentAngle);

    scalar rhoP =
        (
            (sqr(sqrt(n1)*cosTheta1 - sqrt(sqrP)) + sqrQ)
           /
            (sqr(sqrt(n1)*cosTheta1 + sqrt(sqrP)) + sqrQ)
        );

    scalar rhoN =
        (
            (sqr(sqrt(sqrP) - sqrt(n1)*sinTheta1*tanTheta1) + sqrQ)
           /
            (sqr(sqrt(sqrP) + sqrt(n1)*sinTheta1*tanTheta1) + sqrQ)
        )*rhoP;

    return 0.5*(rhoP + rhoN);
}


Foam::vector Foam::radiation::Fresnel::R
(
    const vector& i,
    const vector& n
) const
{
    return  i + 2.0*(-i & n) * n;
}


// ************************************************************************* //
