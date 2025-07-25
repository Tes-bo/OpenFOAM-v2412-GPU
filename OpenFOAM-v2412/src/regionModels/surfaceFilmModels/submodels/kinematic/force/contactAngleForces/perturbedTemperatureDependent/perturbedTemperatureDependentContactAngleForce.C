/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "perturbedTemperatureDependentContactAngleForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(perturbedTemperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    perturbedTemperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
perturbedTemperatureDependentContactAngleForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_, &film.regionMesh())),
    rndGen_(label(0)),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
~perturbedTemperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField>
perturbedTemperatureDependentContactAngleForce::theta() const
{
    auto ttheta = volScalarField::New
    (
        IOobject::scopedName(typeName, "theta"),
        IOobject::NO_REGISTER,
        filmModel_.regionMesh(),
        dimensionedScalar(dimless, Zero)
    );
    auto& thetaIf = ttheta.ref().internalFieldRef();
    auto& thetaBf = ttheta.ref().boundaryFieldRef();

    const volScalarField& T = filmModel_.T();

    // Initialize with the function of temperature
    thetaIf.field() = thetaPtr_->value(T());

    // Add the stochastic perturbation
    forAll(thetaIf, celli)
    {
        thetaIf[celli] += distribution_->sample();
    }

    forAll(thetaBf, patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            fvPatchField<scalar>& thetaf = thetaBf[patchi];

            // Initialize with the function of temperature
            thetaf = thetaPtr_->value(T.boundaryField()[patchi]);

            // Add the stochastic perturbation
            forAll(thetaf, facei)
            {
                thetaf[facei] += distribution_->sample();
            }
        }
    }

    return ttheta;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
