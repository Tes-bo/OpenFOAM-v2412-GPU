/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "adjointOutletWaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointOutletWaFvPatchScalarField::adjointOutletWaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointScalarBoundaryCondition(p, iF, word::null)
{}


adjointOutletWaFvPatchScalarField::adjointOutletWaFvPatchScalarField
(
    const adjointOutletWaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


adjointOutletWaFvPatchScalarField::adjointOutletWaFvPatchScalarField
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


adjointOutletWaFvPatchScalarField::adjointOutletWaFvPatchScalarField
(
    const adjointOutletWaFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    adjointScalarBoundaryCondition(tppsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjointOutletWaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    vectorField nf(patch().nf());

    const fvPatchField<vector>& Ub = boundaryContrPtr_->Ub();
    tmp<scalarField> tnuEff(boundaryContrPtr_->TMVariable2Diffusion());
    const scalarField& nuEff = tnuEff();

    // Patch-adjacent wa
    tmp<scalarField> twaNei = patchInternalField();
    const scalarField& waNei = twaNei();

    const scalarField& delta = patch().deltaCoeffs();

    // Source from the objective
    tmp<scalarField> source = boundaryContrPtr_->adjointTMVariable2Source();

    operator==
    (
        (nuEff*delta*waNei - source)
       /((Ub & nf) + nuEff*delta)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void adjointOutletWaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, adjointOutletWaFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
