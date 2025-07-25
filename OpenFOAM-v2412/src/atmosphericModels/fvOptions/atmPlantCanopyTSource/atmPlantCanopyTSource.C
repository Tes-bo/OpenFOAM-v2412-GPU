/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
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

#include "atmPlantCanopyTSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmPlantCanopyTSource, 0);
    addToRunTimeSelectionTable(option, atmPlantCanopyTSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmPlantCanopyTSource::atmPlantCanopyTSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    Cp0_
    (
        dimensionedScalar
        (
            sqr(dimLength)/sqr(dimTime)/dimTemperature,
            coeffs_.getCheckOrDefault<scalar>
            (
                "Cp0",
                1005.0,
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    qPlant_
    (
        IOobject
        (
            "qPlant",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    fieldNames_.resize(1, "T");

    fv::option::resetApplied();

    Log << "    Applying atmPlantCanopyTSource to: " << fieldNames_[0] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmPlantCanopyTSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (V_ > VSMALL)
    {
        eqn -= -qPlant_/Cp0_;
    }
}


void Foam::fv::atmPlantCanopyTSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (V_ > VSMALL)
    {
        eqn -= -rho*qPlant_/Cp0_;
    }
}


void Foam::fv::atmPlantCanopyTSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (V_ > VSMALL)
    {
        eqn -= -alpha*rho*qPlant_/Cp0_;
    }
}


// ************************************************************************* //
