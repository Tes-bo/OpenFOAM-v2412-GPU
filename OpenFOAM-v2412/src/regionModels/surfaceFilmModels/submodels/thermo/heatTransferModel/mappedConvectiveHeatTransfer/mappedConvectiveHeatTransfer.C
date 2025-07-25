/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "mappedConvectiveHeatTransfer.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mappedConvectiveHeatTransfer, 0);

addToRunTimeSelectionTable
(
    heatTransferModel,
    mappedConvectiveHeatTransfer,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::mappedConvectiveHeatTransfer
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    heatTransferModel(film),
    htcConvPrimary_
    (
        IOobject
        (
            "htcConv",
            film.time().timeName(),
            film.primaryMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.primaryMesh()
    ),
    htcConvFilm_
    (
        IOobject
        (
            htcConvPrimary_.name(), // must have same name as above for mapping
            film.time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime)/dimTemperature, Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    )
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Pull the data from the primary region via direct mapped BCs
    htcConvFilm_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::~mappedConvectiveHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void mappedConvectiveHeatTransfer::correct()
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Pull the data from the primary region via direct mapped BCs
    htcConvFilm_.correctBoundaryConditions();
}


tmp<volScalarField> mappedConvectiveHeatTransfer::h() const
{
    return htcConvFilm_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
