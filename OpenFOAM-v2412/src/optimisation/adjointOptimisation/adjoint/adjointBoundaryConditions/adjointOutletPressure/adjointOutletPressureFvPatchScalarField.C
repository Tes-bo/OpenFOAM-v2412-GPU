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

#include "adjointOutletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "emptyFvPatch.H"
#include "ATCUaGradU.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletPressureFvPatchScalarField::
adjointOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointScalarBoundaryCondition(p, iF, word::null)
{}


Foam::adjointOutletPressureFvPatchScalarField::
adjointOutletPressureFvPatchScalarField
(
    const adjointOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointOutletPressureFvPatchScalarField::
adjointOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointScalarBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::adjointOutletPressureFvPatchScalarField::
adjointOutletPressureFvPatchScalarField
(
    const adjointOutletPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    adjointScalarBoundaryCondition(tppsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch normal and surface
    const scalarField& magSf = patch().magSf();
    const vectorField nf(patch().nf());

    // Primal flux
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    const fvPatchField<vector>& Up = boundaryContrPtr_->Ub();

    // Adjoint velocity
    const fvPatchField<vector>& Uap = boundaryContrPtr_->Uab();
    scalarField snGradUan(Uap.snGrad() & nf);

    // Patch normal adjoint velocity
    scalarField Uap_n(Uap & nf);

    // Patch normal velocity
    scalarField phiOverSurf(phip/magSf);

    // Momentum diffusion coefficient
    tmp<scalarField> tmomentumDiffusion =
        boundaryContrPtr_->momentumDiffusion();
    const scalarField& momentumDiffusion = tmomentumDiffusion();

    // Part of the diffusive flux related to div(nuEff*dev(grad(Ua).T()))
    const word& UaName = boundaryContrPtr_->Uab().internalField().name();
    tmp<tensorField> tgradUab = computePatchGrad<vector>(UaName);
    const tensorField& gradUab = tgradUab();
    vectorField explDiffusiveFlux
    (
        momentumDiffusion*(gradUab - sphericalTensor::oneThirdI*tr(gradUab))
      & nf
    );
    scalarField normalExplDifFlux(explDiffusiveFlux & nf);

    // Objective function and other explicit contributions
    tmp<scalarField> tsource = boundaryContrPtr_->pressureSource();
    scalarField& source = tsource.ref();

    // Contribution from the ATC part (if UaGradU)
    if (addATCUaGradUTerm())
    {
        source += Uap & Up;
    }

    operator==
    (
        (Uap_n*phiOverSurf)
      + momentumDiffusion*snGradUan
      + normalExplDifFlux
      + source
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
