/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022 PCOpt/NTUA
    Copyright (C) 2014-2022 FOSS GP
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

#include "adjointFarFieldTMVar2FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointFarFieldTMVar2FvPatchScalarField::
adjointFarFieldTMVar2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointScalarBoundaryCondition(p, iF, word::null)
{}


adjointFarFieldTMVar2FvPatchScalarField::
adjointFarFieldTMVar2FvPatchScalarField
(
    const adjointFarFieldTMVar2FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


adjointFarFieldTMVar2FvPatchScalarField::
adjointFarFieldTMVar2FvPatchScalarField
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


adjointFarFieldTMVar2FvPatchScalarField::
adjointFarFieldTMVar2FvPatchScalarField
(
    const adjointFarFieldTMVar2FvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    adjointScalarBoundaryCondition(tppsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjointFarFieldTMVar2FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    vectorField nf(patch().nf());

    tmp<scalarField> tnuEff(boundaryContrPtr_->TMVariable2Diffusion());
    const scalarField& nuEff = tnuEff();
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    const scalarField& magSf = patch().magSf();

    // Patch-adjacent values
    tmp<scalarField> intf(patchInternalField());

    // Patch deltas
    const scalarField& delta = patch().deltaCoeffs();

    operator==
    (
        pos(phip)
       *(
            (nuEff*delta*intf)
           /(phip/magSf + nuEff*delta)
        )
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Field<scalar>> adjointFarFieldTMVar2FvPatchScalarField::
valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // For fixedValue patches
    return tmp<Field<scalar>>::New(neg(phip)*pTraits<scalar>::one);
}


tmp<Field<scalar>> adjointFarFieldTMVar2FvPatchScalarField::
valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // For zeroGradient patches
    return tmp<Field<scalar>>::New(pos(phip)*(*this));
}


void adjointFarFieldTMVar2FvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    adjointFarFieldTMVar2FvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
