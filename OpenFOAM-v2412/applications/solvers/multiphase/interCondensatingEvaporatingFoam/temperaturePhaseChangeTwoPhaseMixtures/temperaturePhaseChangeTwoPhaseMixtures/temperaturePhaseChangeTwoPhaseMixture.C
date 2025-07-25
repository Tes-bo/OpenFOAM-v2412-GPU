/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "temperaturePhaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(temperaturePhaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable
    (
        temperaturePhaseChangeTwoPhaseMixture,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeTwoPhaseMixture::
temperaturePhaseChangeTwoPhaseMixture
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseChangeProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mixture_(mixture),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff
    (
        1.0/mixture_.rho1() - mixture_.alpha1()
       *(1.0/mixture_.rho1() - 1.0/mixture_.rho2())
    );

    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixture::vDot() const
{
    dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());
    Pair<tmp<volScalarField>> mDot = this->mDot();

    return Pair<tmp<volScalarField>>(pCoeff*mDot[0], pCoeff*mDot[1]);
}


bool Foam::temperaturePhaseChangeTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
