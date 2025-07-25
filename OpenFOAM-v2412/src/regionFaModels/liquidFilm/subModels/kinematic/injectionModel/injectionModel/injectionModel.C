/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(injectionModel, 0);
defineRunTimeSelectionTable(injectionModel, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void injectionModel::addToInjectedMass(const scalar dMass)
{
    injectedMass_ += dMass;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

injectionModel::injectionModel(liquidFilmBase& film)
:
    filmSubModelBase(film),
    injectedMass_(0.0)
{}


injectionModel::injectionModel
(
    const word& modelType,
    liquidFilmBase& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType),
    injectedMass_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

injectionModel::~injectionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void injectionModel::correct()
{
    if (writeTime())
    {
        scalar injectedMass0 = getModelProperty<scalar>("injectedMass");
        injectedMass0 += returnReduce(injectedMass_, sumOp<scalar>());
        setModelProperty<scalar>("injectedMass", injectedMass0);
        injectedMass_ = 0.0;
    }
}


scalar injectionModel::injectedMassTotal() const
{
    const scalar injectedMass0 = getModelProperty<scalar>("injectedMass");
    return injectedMass0 + returnReduce(injectedMass_, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
