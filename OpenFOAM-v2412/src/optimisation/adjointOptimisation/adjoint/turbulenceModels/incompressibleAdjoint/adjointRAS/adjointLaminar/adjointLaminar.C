/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "adjointLaminar.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{
namespace adjointRASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointLaminar, 0);
addToRunTimeSelectionTable(adjointRASModel, adjointLaminar, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointLaminar::adjointLaminar
(
    incompressibleVars& primalVars,
    incompressibleAdjointMeanFlowVars& adjointVars,
    objectiveManager& objManager,
    const word& adjointTurbulenceModelName,
    const word& modelName
)
:
    adjointRASModel
    (
        modelName,
        primalVars,
        adjointVars,
        objManager,
        adjointTurbulenceModelName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> adjointLaminar::devReff() const
{
    const volVectorField& Ua = adjointVars_.Ua();
    return devReff(Ua);
}


tmp<volSymmTensorField> adjointLaminar::devReff
(
    const volVectorField& U
) const
{
    return volSymmTensorField::New
    (
        "devRhoReff",
        IOobject::NO_REGISTER,
        (
           -nu()*devTwoSymm(fvc::grad(U))
        )
    );
}


tmp<fvVectorMatrix> adjointLaminar::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<volVectorField> adjointLaminar::adjointMeanFlowSource()
{
    return volVectorField::New
    (
        "adjointMeanFlowSource",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(dimensionSet(0, 1, -2, 0, 0), Zero)
    );
}


const boundaryVectorField& adjointLaminar::adjointMomentumBCSource() const
{
    // zero contribution
    return adjMomentumBCSourcePtr_();
}


const boundaryVectorField& adjointLaminar::wallShapeSensitivities()
{
    return wallShapeSensitivitiesPtr_();
}


const boundaryVectorField& adjointLaminar::wallFloCoSensitivities()
{
    return wallFloCoSensitivitiesPtr_();
}


tmp<volScalarField> adjointLaminar::distanceSensitivities()
{
    return volScalarField::New
    (
        "adjointEikonalSource" + type(),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimLength/pow3(dimTime), Zero)
    );
}


tmp<volTensorField> adjointLaminar::FISensitivityTerm()
{
    return volTensorField::New
    (
        "volumeSensTerm" + type(),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedTensor(dimensionSet(0, 2, -3, 0, 0), Zero)
    );
}


tmp<scalarField> adjointLaminar::topologySensitivities
(
    const word& designVarsName
) const
{
    return tmp<scalarField>::New(mesh_.nCells(), Zero);
}


void adjointLaminar::nullify()
{
    // Does nothing. No fields to nullify
}


bool adjointLaminar::read()
{
    return adjointRASModel::read();
}


void adjointLaminar::correct()
{
    adjointTurbulenceModel::correct();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adjointRASModels
} // End namespace incompressibleAdjoint
} // End namespace Foam

// ************************************************************************* //
