/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "constantRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    constantRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantRadiation::constantRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, film, dict),
    qrConst_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "qrConst"),
            film.time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    ),
    mask_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "mask"),
            film.time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    absorptivity_(coeffDict_.get<scalar>("absorptivity")),
    timeStart_(coeffDict_.get<scalar>("timeStart")),
    duration_(coeffDict_.get<scalar>("duration"))
{
    mask_ = pos0(mask_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantRadiation::~constantRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void constantRadiation::correct()
{}


tmp<volScalarField> constantRadiation::Shs()
{
    auto tShs = volScalarField::New
    (
        IOobject::scopedName(typeName, "Shs"),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    );
    auto& Shs = tShs.ref().primitiveFieldRef();

    const scalar time = film().time().value();

    if ((time >= timeStart_) && (time <= timeStart_ + duration_))
    {
        const scalarField& qr = qrConst_;
        const scalarField& alpha = filmModel_.alpha();

        Shs = mask_*qr*alpha*absorptivity_;
    }

    return tShs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
