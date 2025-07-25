/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "shallowWaterAbsorption.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(shallowWaterAbsorption, 0);
    addToRunTimeSelectionTable(waveModel, shallowWaterAbsorption, patch);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::waveModels::shallowWaterAbsorption::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    level = waterDepthRef_;
}


void Foam::waveModels::shallowWaterAbsorption::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    // Apply zero-gradient condition to z-component of velocity only
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    U_ = U.boundaryField()[patch_.index()].patchInternalField();

    for (vector& vel : U_)
    {
        vel.x() = 0;
        vel.y() = 0;
    }
}


void Foam::waveModels::shallowWaterAbsorption::setAlpha
(
    const scalarField& level
)
{
    // Set alpha as zero-gradient
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    alpha_ = alpha.boundaryField()[patch_.index()].patchInternalField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::shallowWaterAbsorption::shallowWaterAbsorption
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveAbsorptionModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::shallowWaterAbsorption::readDict
(
    const dictionary& overrideDict
)
{
    return waveAbsorptionModel::readDict(overrideDict);
}


// ************************************************************************* //
