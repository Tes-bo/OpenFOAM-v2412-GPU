/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "incompressibleRhoTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleRhoTurbulenceModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::incompressibleRhoTurbulenceModel::incompressibleRhoTurbulenceModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    turbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    ),
    rho_(rho)
{}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleRhoTurbulenceModel::mu() const
{
    return rho_*nu();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleRhoTurbulenceModel::mu(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleRhoTurbulenceModel::mut() const
{
    return rho_*nut();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleRhoTurbulenceModel::mut(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleRhoTurbulenceModel::muEff() const
{
    return rho_*nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleRhoTurbulenceModel::muEff(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nuEff(patchi);
}


// ************************************************************************* //
