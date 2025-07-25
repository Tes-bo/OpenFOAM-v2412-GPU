/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "BrunDrippingInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BrunDrippingInjection, 0);
addToRunTimeSelectionTable
(
    injectionModel,
    BrunDrippingInjection,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BrunDrippingInjection::BrunDrippingInjection
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    ubarStar_
    (
        coeffDict_.getCheckOrDefault<scalar>
        (
            "ubarStar",
            1.62208,
            scalarMinMax::ge(SMALL)
        )
    ),
    dCoeff_(coeffDict_.getOrDefault<scalar>("dCoeff", 3.3)),
    deltaStable_(coeffDict_.getOrDefault<scalar>("deltaStable", 0)),
    diameter_(film.regionMesh().nFaces(), -1.0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void BrunDrippingInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    // Calculate available dripping mass
    tmp<areaScalarField> tsinAlpha = -film().gn()/mag(film().g());
    const scalarField& sinAlpha = tsinAlpha();

    const areaScalarField& delta = film().h();
    const areaScalarField& rho = film().rho();
    const areaScalarField& sigma = film().sigma();
    const scalar magg = mag(film().g().value());

    forAll(delta, facei)
    {
        bool dripping = false;

        if (sinAlpha[facei] > SMALL && delta[facei] > deltaStable_)
        {
            const scalar rhoc = rho[facei];
            const scalar lc = sqrt(sigma[facei]/(rhoc*magg));
            const scalar deltaStable = max
            (
                3*lc*sqrt(1 - sqr(sinAlpha[facei]))
               /(ubarStar_*sqrt(sinAlpha[facei])*sinAlpha[facei]),
                deltaStable_
            );

            if (delta[facei] > deltaStable)
            {
                const scalar massDrip =
                    availableMass[facei]*(delta[facei] - deltaStable);

                if (massDrip > 0)
                {
                    const scalar diam = dCoeff_*lc;
                    diameter_[facei] = diam;

                    massToInject[facei] += massDrip;
                    availableMass[facei] -= massDrip;

                    diameterToInject[facei] = diam;
                    addToInjectedMass(massDrip);

                    dripping = true;
                }
            }
        }

        if (!dripping)
        {
            diameterToInject[facei] = 0;
            massToInject[facei] = 0;
        }
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
