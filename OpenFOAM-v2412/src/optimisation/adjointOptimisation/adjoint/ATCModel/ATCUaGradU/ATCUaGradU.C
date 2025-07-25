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

#include "ATCUaGradU.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ATCUaGradU, 0);
addToRunTimeSelectionTable
(
    ATCModel,
    ATCUaGradU,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ATCUaGradU::ATCUaGradU
(
    const fvMesh& mesh,
    const incompressibleVars& primalVars,
    const incompressibleAdjointVars& adjointVars,
    const dictionary& dict
)
:
    ATCModel(mesh, primalVars, adjointVars, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ATCUaGradU::addATC(fvVectorMatrix& UaEqn)
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.UaInst();
    const surfaceScalarField& phi = primalVars_.phi();
    const surfaceScalarField& phia = adjointVars_.phiaInst();

    // Build Ua to go into the ATC term, based on whether to smooth field or not
    autoPtr<volVectorField> UaForATC(nullptr);
    if (reconstructGradients_)
    {
        UaForATC.reset(new volVectorField(fvc::reconstruct(phia)));
    }
    else
    {
        UaForATC.reset(new volVectorField(Ua));
    }

    // Main ATC term
    ATC_ = -fvc::grad(UaForATC(), "gradUaATC") & U;

    if (extraConvection_ > 0)
    {
        // Implicit part added to increase diagonal dominance
        UaEqn += ATClimiter_*extraConvection_*fvm::div(-phi, Ua);

        // Correct rhs due to implicitly augmenting the adjoint convection
        ATC_ += extraConvection_*(fvc::grad(UaForATC(), "gradUaATC")().T() & U);
    }

    // Zero ATC on cells next to given patch types
    smoothATC();

    // Actual ATC term
    UaEqn += ATC_.internalField();
}


tmp<volTensorField> ATCUaGradU::getFISensitivityTerm() const
{
    auto tvolSDTerm = volTensorField::New
    (
        "ATCFISensitivityTerm" + type(),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedTensor(sqr(dimLength)/pow(dimTime, 3), Zero)
    );
    auto& volSDTerm = tvolSDTerm.ref();

    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.Ua();

    //const volScalarField& mask = getLimiter();

    volSDTerm -=
        Ua.component(0) * fvc::grad(U.component(0) * U)
      + Ua.component(1) * fvc::grad(U.component(1) * U)
      + Ua.component(2) * fvc::grad(U.component(2) * U);

    return tvolSDTerm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
