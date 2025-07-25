/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "adjointOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFvPatchVectorField::assignBoundaryValue()
{
    const scalarField& magSf = patch().magSf();
    tmp<vectorField> tnf(patch().nf());
    const vectorField& nf = tnf();

    // Primal normal velocity
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    const scalarField phiOverSurf(phip/magSf);

    // Ua patch adjacent
    vectorField Uac(this->patchInternalField());

    // Tangent component of internalField
    vectorField Uac_t(Uac - nf*(Uac & nf));

    // Adjoint normal velocity
    const fvsPatchField<scalar>& phiab = boundaryContrPtr_->phiab();

    // Inverse distance
    const scalarField& delta = patch().deltaCoeffs();

    // Objective function and other explicit contributions
    tmp<vectorField> tsource(boundaryContrPtr_->tangentVelocitySource());
    const vectorField& source = tsource();

    // Momentum diffusion coefficient
    tmp<scalarField> tmomentumDiffusion
    (
        boundaryContrPtr_->momentumDiffusion()
    );
    const scalarField& momentumDiffusion = tmomentumDiffusion();

    // Part of the diffusive flux related to div(nuEff*dev(grad(Ua).T()))
    const word& fieldName = internalField().name();
    tmp<tensorField> tgradUaf(computePatchGrad<vector>(fieldName));
    const tensorField& gradUaf = tgradUaf();
    const vectorField explDiffusiveFlux
    (
        momentumDiffusion
       *(gradUaf - sphericalTensor::oneThirdI*tr(gradUaf)) & nf
    );
    const vectorField explDiffusiveFlux_t
    (
        explDiffusiveFlux - (explDiffusiveFlux & nf)*nf
    );

    // Auxiliary quantities
    scalarField nd(momentumDiffusion*delta);
    // Denominator. Susceptible to zero values in case of back flow
    // Should use adjointOutletVelocityFlux in such cases.
    scalarField denom(phiOverSurf + nd);

    vectorField Uat((nd*Uac_t - explDiffusiveFlux_t - source)/denom);

    operator==((phiab/magSf)*nf + Uat);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, word::null)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointVectorBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointVectorBoundaryCondition(pivpvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    assignBoundaryValue();
    fvPatchVectorField::evaluate();
}


void Foam::adjointOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::adjointOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
