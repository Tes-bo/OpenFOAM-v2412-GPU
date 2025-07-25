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

#include "adjointFarFieldVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointFarFieldVelocityFvPatchVectorField::
adjointFarFieldVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, word::null)
{}


Foam::adjointFarFieldVelocityFvPatchVectorField::
adjointFarFieldVelocityFvPatchVectorField
(
    const adjointFarFieldVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointVectorBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointFarFieldVelocityFvPatchVectorField::
adjointFarFieldVelocityFvPatchVectorField
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


Foam::adjointFarFieldVelocityFvPatchVectorField::
adjointFarFieldVelocityFvPatchVectorField
(
    const adjointFarFieldVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointVectorBoundaryCondition(pivpvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointFarFieldVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& faceMag = patch().magSf();
    const vectorField nf(patch().nf());

    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    scalarField phiOverSurf(phip/faceMag);

    // Ua patch adjacent
    vectorField Uac(this->patchInternalField());

    // Tangent component of internalField
    vectorField Uac_t(Uac - nf*(Uac & nf));

    // Inverse distance
    const scalarField& delta = patch().deltaCoeffs();

    // Objective function and other explicit contributions for
    // zeroGradient boundaries
    tmp<vectorField> tsourceVelocity =
        boundaryContrPtr_->tangentVelocitySource();
    vectorField& sourceVelocity = tsourceVelocity.ref();

    // Objective function contribution for fixedValue boundaries
    tmp<vectorField> tsourcePressure =
        boundaryContrPtr_->normalVelocitySource();
    vectorField& sourcePressure = tsourcePressure.ref();

    // Momentum diffusion coefficient
    tmp<scalarField> tmomentumDiffusion =
        boundaryContrPtr_->momentumDiffusion();
    scalarField& momentumDiffusion = tmomentumDiffusion.ref();

    scalarField denom(phiOverSurf + momentumDiffusion*delta);

    operator==
    (
        // Inlet
      - neg(phip)*sourcePressure

        // Outlet
      + pos(phip)
       *((Uac&nf)*nf + (Uac_t*(momentumDiffusion*delta) - sourceVelocity)/denom)
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointFarFieldVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // For fixedValue U patches
    return tmp<Field<vector>>
    (
        new Field<vector>
        (
            neg(phip)*pTraits<vector>::one
        )
    );
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointFarFieldVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // For zeroGradient U patches
    return tmp<Field<vector>>
    (
        new Field<vector>
        (
            pos(phip)*(*this)
        )
    );
}


void Foam::adjointFarFieldVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointFarFieldVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
