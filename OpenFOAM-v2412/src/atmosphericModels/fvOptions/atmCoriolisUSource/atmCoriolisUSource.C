/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 CENER
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "atmCoriolisUSource.H"
#include "fvMatrices.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmCoriolisUSource, 0);
    addToRunTimeSelectionTable(option, atmCoriolisUSource, dictionary);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::fv::atmCoriolisUSource::planetaryRotationVector() const
{
    return vector
    (
        Zero,
        twoPi/(planetaryRotationPeriod_*3600.0)*cos(degToRad(latitude_)),
        twoPi/(planetaryRotationPeriod_*3600.0)*sin(degToRad(latitude_))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmCoriolisUSource::atmCoriolisUSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    latitude_
    (
        coeffs_.getCheckOrDefault<scalar>
        (
            "latitude",
            0.0,
            [&](const scalar x){ return (90 >= mag(x)) && (mag(x) >= 0); }
        )
    ),
    planetaryRotationPeriod_
    (
        coeffs_.getCheckOrDefault<scalar>
        (
            "planetaryRotationPeriod",
            23.9344694,
            [&](const scalar x){ return x > SMALL; }
        )
    ),
    Omega_
    (
        dimensionedVector
        (
            dimless/dimTime,
            coeffs_.getOrDefault<vector>
            (
                "Omega",
                planetaryRotationVector()
            )
        )
    )
{
    if (mag(Omega_.value()) < SMALL)
    {
        WarningInFunction
            << "The magnitude of the rotation vector in atmCoriolisUSource is "
            << "effectively zero, mag(Omega) = " << mag(Omega_.value()) << nl
            << "Please check input values in atmCoriolisUSource settings."
            << endl;
    }

    fieldNames_.resize(1, "U");

    fv::option::resetApplied();

    Log << "    Applying atmCoriolisUSource to: " << fieldNames_[0] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmCoriolisUSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        eqn -= (2.0*Omega_)^U;
    }
}


void Foam::fv::atmCoriolisUSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        eqn -= rho*((2.0*Omega_)^U);
    }
}


void Foam::fv::atmCoriolisUSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    if (V_ > VSMALL)
    {
        eqn -= alpha*rho*((2.0*Omega_)^U);
    }
}


// ************************************************************************* //
