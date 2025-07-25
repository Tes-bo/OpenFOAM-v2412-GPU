/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "filmSeparation.H"
#include "filmSeparationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmSeparation, 0);
addToRunTimeSelectionTable
(
    injectionModel,
    filmSeparation,
    dictionary
);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::filmSeparation::filmSeparation
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    filmSeparationModelPtr_(filmSeparationModel::New(film, coeffDict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::filmSeparation::~filmSeparation()
{}  // filmSeparationModel was forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::areaSurfaceFilmModels::filmSeparation::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const faMesh& mesh = film().regionMesh();

    // Calculate the mass ratio of film separation
    tmp<scalarField> tmassRatio = filmSeparationModelPtr_->separatedMassRatio();
    const auto& massRatio = tmassRatio.cref();

    // Update various film properties based on the mass ratio
    massToInject = massRatio*availableMass;
    diameterToInject = massRatio*film().h();
    availableMass -= massRatio*availableMass;

    addToInjectedMass(sum(massToInject));

    injectionModel::correct();


    if (debug && mesh.time().writeTime())
    {
        {
            areaScalarField areaSeparated
            (
                mesh.newIOobject("separated"),
                mesh,
                dimensionedScalar(dimMass, Zero)
            );
            areaSeparated.primitiveFieldRef() = massRatio;
            areaSeparated.write();
        }

        {
            areaScalarField areaMassToInject
            (
                mesh.newIOobject("massToInject"),
                mesh,
                dimensionedScalar(dimMass, Zero)
            );
            areaMassToInject.primitiveFieldRef() = massToInject;
            areaMassToInject.write();
        }
    }
}


// ************************************************************************* //
