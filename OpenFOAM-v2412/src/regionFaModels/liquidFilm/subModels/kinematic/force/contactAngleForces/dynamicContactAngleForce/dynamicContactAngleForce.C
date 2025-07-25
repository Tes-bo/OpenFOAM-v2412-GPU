/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "dynamicContactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "Function1.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    dynamicContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicContactAngleForce::dynamicContactAngleForce
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    U_vs_thetaPtr_
    (
        Function1<scalar>::NewIfPresent
        (
            "Utheta",
            coeffDict_,
            word::null,
            &film.primaryMesh()
        )
    ),
    T_vs_thetaPtr_
    (
        Function1<scalar>::NewIfPresent
        (
            "Ttheta",
            coeffDict_,
            word::null,
            &film.primaryMesh()
        )
    ),
    rndGen_(label(0)),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    )
{
    if (U_vs_thetaPtr_ && T_vs_thetaPtr_)
    {
        FatalIOErrorInFunction(dict)
            << "Entries Utheta and Ttheta could not be used together"
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<areaScalarField> dynamicContactAngleForce::theta() const
{
    auto ttheta = tmp<areaScalarField>::New
    (
        IOobject
        (
            IOobject::scopedName(typeName, "theta"),
            film().regionMesh().time().timeName(),
            film().regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        film().regionMesh(),
        dimensionedScalar(dimless, Zero)
    );
    areaScalarField& theta = ttheta.ref();
    scalarField& thetai = theta.ref();


    if (U_vs_thetaPtr_)
    {
        // Initialize with the function of film speed
        const areaVectorField& U = film().Uf();

        thetai = U_vs_thetaPtr_->value(mag(U()));
    }

    if (T_vs_thetaPtr_)
    {
        // Initialize with the function of film temperature
        const areaScalarField& T = film().Tf();

        thetai = T_vs_thetaPtr_->value(T());
    }

    // Add the stochastic perturbation
    forAll(thetai, facei)
    {
        thetai[facei] += distribution_->sample();
    }

    return ttheta;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
