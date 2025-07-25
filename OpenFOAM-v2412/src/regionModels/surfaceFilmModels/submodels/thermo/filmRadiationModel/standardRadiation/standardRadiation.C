/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "standardRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    standardRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardRadiation::standardRadiation
(
     surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, film, dict),
    qinPrimary_
    (
        IOobject
        (
            "qin", // same name as qin on primary region to enable mapping
            film.time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    ),
    qrNet_
    (
        IOobject
        (
            "qrNet",
            film.time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    beta_(coeffDict_.get<scalar>("beta")),
    kappaBar_(coeffDict_.get<scalar>("kappaBar"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardRadiation::correct()
{
    // Transfer qr from primary region
    qinPrimary_.correctBoundaryConditions();
}


tmp<volScalarField> standardRadiation::Shs()
{
    auto tShs = volScalarField::New
    (
        IOobject::scopedName(typeName, "Shs"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    );
    scalarField& Shs = tShs.ref();

    const scalarField& qinP = qinPrimary_;
    const scalarField& delta = filmModel_.delta();
    const scalarField& alpha = filmModel_.alpha();

    Shs = beta_*qinP*alpha*(1.0 - exp(-kappaBar_*delta));

    // Update net qr on local region
    qrNet_.primitiveFieldRef() = qinP - Shs;
    qrNet_.correctBoundaryConditions();

    return tShs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
