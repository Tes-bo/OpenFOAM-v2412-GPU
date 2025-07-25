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

#include "multiphaseInterHtcModel.H"
#include "heatTransferCoeffModel.H"
#include "multiphaseInterSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(multiphaseInterHtcModel, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        multiphaseInterHtcModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::multiphaseInterHtcModel::calc()
{
    auto& htc =
        htcModelPtr_->mesh().lookupObjectRef<volScalarField>(resultName_);

    htcModelPtr_->calc(htc, q());

    return true;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::functionObjects::multiphaseInterHtcModel::q() const
{
    const fvMesh& mesh = htcModelPtr_->mesh();

    const auto& T = mesh.lookupObject<volScalarField>(htcModelPtr_->TName());

    const volScalarField::Boundary& Tbf = T.boundaryField();

    auto tq = tmp<FieldField<Field, scalar>>::New(Tbf.size());
    auto& q = tq.ref();

    forAll(q, patchi)
    {
        q.set(patchi, new Field<scalar>(Tbf[patchi].size(), Zero));
    }

    const auto* fluidPtr =
        mesh.cfindObject<multiphaseInterSystem>("phaseProperties");

    if (!fluidPtr)
    {
        FatalErrorInFunction
            << "Unable to find a valid phaseSystem to evaluate q" << nl
            << exit(FatalError);
    }

    const multiphaseInterSystem& fluid = *fluidPtr;

    for (const label patchi : htcModelPtr_->patchIDs())
    {
        q[patchi] += fluid.kappaEff(patchi)()*Tbf[patchi].snGrad();
    }

    // Add radiative heat flux contribution if present

    const auto* qrPtr =
        mesh.cfindObject<volScalarField>(htcModelPtr_->qrName());

    if (qrPtr)
    {
        const volScalarField::Boundary& qrbf = qrPtr->boundaryField();

        for (const label patchi : htcModelPtr_->patchIDs())
        {
            q[patchi] += qrbf[patchi];
        }
    }

    return tq;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::multiphaseInterHtcModel::multiphaseInterHtcModel
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    htcModelPtr_(heatTransferCoeffModel::New(dict, mesh_, fieldName_))
{
    read(dict);

    // scopedName?
    setResultName(typeName, "htc:" + htcModelPtr_->type());

    auto* htcPtr =
        new volScalarField
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimPower/dimArea/dimTemperature, Zero)
        );

    regIOobject::store(htcPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::multiphaseInterHtcModel::read
(
    const dictionary& dict
)
{
    if (!fieldExpression::read(dict) || !htcModelPtr_->read(dict))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
