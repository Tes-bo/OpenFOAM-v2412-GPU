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

#include "temperatureDependentContactAngleForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    temperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::temperatureDependentContactAngleForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_, &film.regionMesh()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::~temperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> temperatureDependentContactAngleForce::theta() const
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

    thetaIf.field() = thetaPtr_->value(T());

    forAll(thetaBf, patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            thetaBf[patchi] = thetaPtr_->value(T.boundaryField()[patchi]);
        }
    }

    return ttheta;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
