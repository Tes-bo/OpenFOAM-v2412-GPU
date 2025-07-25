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

#include "laminar.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminar, 0);
addToRunTimeSelectionTable(filmTurbulenceModel, laminar, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminar::laminar
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    filmTurbulenceModel(type(), film, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<areaScalarField> laminar::mut() const
{
    auto tmut = tmp<areaScalarField>::New
    (
        IOobject
        (
            "mut",
            film().regionMesh().time().timeName(),
            film().regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film().regionMesh(),
        dimensionedScalar(dimMass/dimLength/dimTime)
    );

    return tmut;
}


void laminar::correct()
{}


tmp<faVectorMatrix> laminar::Su(areaVectorField& U) const
{
    return primaryRegionFriction(U) + wallFriction(U);
}


tmp<faVectorMatrix> laminar::wallFriction(areaVectorField& U) const
{
    // local references to film fields
    tmp<areaVectorField> Uw = film_.Uw();
    tmp<areaScalarField> wf = Cw();

    return
    (
       - fam::Sp(wf(), U) + wf()*Uw() // wall contribution
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
