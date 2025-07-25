/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "ATCstandard.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ATCstandard, 0);
addToRunTimeSelectionTable
(
    ATCModel,
    ATCstandard,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ATCstandard::ATCstandard
(
    const fvMesh& mesh,
    const incompressibleVars& primalVars,
    const incompressibleAdjointVars& adjointVars,
    const dictionary& dict
)
:
    ATCModel(mesh, primalVars, adjointVars, dict),
    gradU_
    (
        IOobject
        (
            "gradUATC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimless/dimTime, Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ATCstandard::addATC(fvVectorMatrix& UaEqn)
{
    addProfiling(ATCstandard, "ATCstandard::addATC");
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.UaInst();
    const surfaceScalarField& phi = primalVars_.phi();


    // Main ATC term
    ATC_ = gradU_ & Ua;

    if (extraConvection_ > 0)
    {
        // Implicit part added to increase diagonal dominance
        UaEqn += ATClimiter_*extraConvection_*fvm::div(-phi, Ua);

        // Correct rhs due to implicitly augmenting the adjoint convection
        ATC_ += extraConvection_*(fvc::grad(Ua, "gradUaATC")().T() & U);
    }

    //zero ATC on cells next to given patch types
    smoothATC();

    // actual ATC term
    UaEqn += ATC_.internalField();
}


tmp<volTensorField> ATCstandard::getFISensitivityTerm() const
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.Ua();

    // We only need to modify the boundaryField of gradU locally.
    // If grad(U) is cached then
    // a. The .ref() call fails since the tmp is initialised from a
    //    const ref
    // b. we would be changing grad(U) for all other places in the code
    //    that need it
    // So, always allocate new memory and avoid registering the new field
    tmp<volTensorField> tgradU(volTensorField::New("gradULocal", fvc::grad(U)));
    volTensorField::Boundary& gradUbf = tgradU.ref().boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            gradUbf[patchI] = tnf*U.boundaryField()[patchI].snGrad();
        }
    }

    const volScalarField& mask = getLimiter();

    return
        tmp<volTensorField>::New
        (
            "ATCFISensitivityTerm" + type(),
          - (tgradU & Ua)*U*mask
        );
}


void ATCstandard::updatePrimalBasedQuantities()
{
    const volVectorField& U = primalVars_.U();
    const surfaceScalarField& phi = primalVars_.phi();
    // Build U to go into the ATC term, based on whether to smooth field or not
    autoPtr<volVectorField> UForATC(nullptr);
    if (reconstructGradients_)
    {
        gradU_ = fvc::grad(fvc::reconstruct(phi), "gradUATC");
    }
    else
    {
        gradU_ = fvc::grad(U, "gradUATC");
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
