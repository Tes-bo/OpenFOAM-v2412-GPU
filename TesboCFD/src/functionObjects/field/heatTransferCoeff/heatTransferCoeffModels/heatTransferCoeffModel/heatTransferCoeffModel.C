/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "heatTransferCoeffModel.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferCoeffModel, 0);
    defineRunTimeSelectionTable(heatTransferCoeffModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModel::heatTransferCoeffModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    mesh_(mesh),
    TName_(TName),
    qrName_("qr")
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::heatTransferCoeffModel::q() const
{
    const auto& T = mesh_.lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();

    auto tq = tmp<FieldField<Field, scalar>>::New(Tbf.size());
    auto& q = tq.ref();

    forAll(q, patchi)
    {
        q.set(patchi, new Field<scalar>(Tbf[patchi].size(), Zero));
    }

    typedef compressible::turbulenceModel cmpTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        // Note: calling he(p,T) instead of he()
        const volScalarField he(turb.transport().he(turb.transport().p(), T));
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField alphaEff(turb.alphaEff());
        const volScalarField::Boundary& alphaEffbf = alphaEff.boundaryField();

        for (const label patchi : patchIDs_)
        {
            q[patchi] = alphaEffbf[patchi]*hebf[patchi].snGrad();
        }
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        // Note: calling he(p,T) instead of he()
        const volScalarField he(thermo.he(thermo.p(), T));
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField& alpha(thermo.alpha());
        const volScalarField::Boundary& alphabf = alpha.boundaryField();

        for (const label patchi : patchIDs_)
        {
            q[patchi] = alphabf[patchi]*hebf[patchi].snGrad();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find a valid thermo model to evaluate q. " << nl
            << "Database contents are: " << mesh_.objectRegistry::sortedToc()
            << exit(FatalError);
    }

    // Add radiative heat flux contribution if present

    const auto* qrPtr = mesh_.cfindObject<volScalarField>(qrName_);

    if (qrPtr)
    {
        const volScalarField::Boundary& qrbf = qrPtr->boundaryField();

        for (const label patchi : patchIDs_)
        {
            q[patchi] += qrbf[patchi];
        }
    }

    return tq;
}


bool Foam::heatTransferCoeffModel::calc
(
    volScalarField& result,
    const FieldField<Field, scalar>& q
)
{
    htc(result, q);

    return true;
}


bool Foam::heatTransferCoeffModel::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    dict.readIfPresent("qr", qrName_);

    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();

    return true;
}


// ************************************************************************* //
