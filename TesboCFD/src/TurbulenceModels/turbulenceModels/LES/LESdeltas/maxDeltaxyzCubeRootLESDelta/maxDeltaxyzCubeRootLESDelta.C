/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "maxDeltaxyzCubeRootLESDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(maxDeltaxyzCubeRootLESDelta, 0);
    addToRunTimeSelectionTable
    (
        LESdelta,
        maxDeltaxyzCubeRootLESDelta,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::calcDelta()
{
    maxDeltaxyz_.calcDelta();
    cubeRootVolDelta_.calcDelta();

    delta_ =
        max
        (
            static_cast<const volScalarField&>(maxDeltaxyz_),
            static_cast<const volScalarField&>(cubeRootVolDelta_)
        );

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::maxDeltaxyzCubeRootLESDelta::maxDeltaxyzCubeRootLESDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    maxDeltaxyz_
    (
        name + "maxDeltaxyz",
        turbulence,
        dict.subDict(typeName + "Coeffs")
    ),
    cubeRootVolDelta_
    (
        name + "cubeRootVolDelta",
        turbulence,
        dict.subDict(typeName + "Coeffs")
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::read(const dictionary& dict)
{
   maxDeltaxyz_.read(dict.subDict(typeName + "Coeffs"));
   cubeRootVolDelta_.read(dict.subDict(typeName + "Coeffs"));

   calcDelta();
}


void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
