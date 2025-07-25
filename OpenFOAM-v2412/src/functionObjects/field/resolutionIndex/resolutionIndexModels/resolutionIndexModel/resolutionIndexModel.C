/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "resolutionIndexModel.H"
#include "fvMesh.H"
#include "ListOps.H"
#include "turbulenceFields.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(resolutionIndexModel, 0);
    defineRunTimeSelectionTable(resolutionIndexModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::resolutionIndexModel::V() const
{
    auto tV = volScalarField::New
    (
        "V",
        IOobject::NO_REGISTER,
        mesh_,
        dimVolume,
        fvPatchFieldBase::zeroGradientType()
    );

    tV.ref().primitiveFieldRef() = mesh_.V();
    tV.ref().correctBoundaryConditions();

    return tV;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resolutionIndexModel::resolutionIndexModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    resultName_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::resolutionIndexModel::read(const dictionary& dict)
{
    resultName_ = dict.getOrDefault("result", type());

    auto* indexPtr = mesh_.getObjectPtr<volScalarField>(resultName_);

    if (!indexPtr)
    {
        indexPtr = new volScalarField
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::LAZY_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        regIOobject::store(indexPtr);
    }

    return true;
}


// ************************************************************************* //
