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

#include "adjointFarFieldPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ATCUaGradU.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointFarFieldPressureFvPatchScalarField::
adjointFarFieldPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointScalarBoundaryCondition(p, iF, word::null)
{}


Foam::adjointFarFieldPressureFvPatchScalarField::
adjointFarFieldPressureFvPatchScalarField
(
    const adjointFarFieldPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointFarFieldPressureFvPatchScalarField::
adjointFarFieldPressureFvPatchScalarField
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


Foam::adjointFarFieldPressureFvPatchScalarField::
adjointFarFieldPressureFvPatchScalarField
(
    const adjointFarFieldPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    adjointScalarBoundaryCondition(tppsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointFarFieldPressureFvPatchScalarField::updateCoeffs()
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

    // Adjoint flux
    //const fvsPatchField<scalar>& phiap =
    //    patch().lookupPatchField<surfaceScalarField>("phia");

    // Primal velocity
    const fvPatchField<vector>& Up = boundaryContrPtr_->Ub();

    // Adjoint velocity
    const fvPatchField<vector>& Uap = boundaryContrPtr_->Uab();

    // Patch-adjacent normal adjoint velocity
    scalarField Uac_n(Uap.patchInternalField()() & nf);

    // Patch normal adjoint velocity
    scalarField Uap_n(Uap & nf);
    //scalarField Uap_n = phiap/magSf;

    // Patch normal velocity Uap_n
    scalarField phiOverSurf(phip/magSf);

    // Patch deltas
    const scalarField& delta = patch().deltaCoeffs();

    // snGrad Ua_n
    scalarField snGradUan(delta*(Uap_n - Uac_n));

    // Momentum diffusion coefficient
    tmp<scalarField> tmomentumDiffusion =
        boundaryContrPtr_->momentumDiffusion();
    scalarField& momentumDiffusion = tmomentumDiffusion.ref();

    // Objective function and other explicit contributions
    tmp<scalarField> tsource = boundaryContrPtr_->pressureSource();
    scalarField source = tsource.ref();

    // Contribution from the ATC part (if UaGradU)
    if (addATCUaGradUTerm())
    {
        source += Uap & Up;
    }

    operator==
    (
        // Inlet
        neg(phip)*(patchInternalField())

        // Outlet
      + pos(phip)*
        (
            Uap_n*phiOverSurf
          + 2*momentumDiffusion*snGradUan
          + source
        )
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::adjointFarFieldPressureFvPatchScalarField::snGrad() const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    return tmp<Field<scalar>>
    (
        new Field<scalar>
        (
            pos(phip)*patch().deltaCoeffs()*(*this - patchInternalField())
        )
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::adjointFarFieldPressureFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    return tmp<Field<scalar>>
    (
        new Field<scalar>
        (
            neg(phip)*pTraits<scalar>::one
        )
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::adjointFarFieldPressureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    return tmp<Field<scalar>>
    (
        new Field<scalar>
        (
            pos(phip)*(*this)
        )
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::adjointFarFieldPressureFvPatchScalarField::gradientInternalCoeffs() const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // Act as a zeroGradient pa bc
    return tmp<Field<scalar>>
    (
        new Field<scalar>
        (
            -pos(phip)*pTraits<scalar>::one*(this->patch().deltaCoeffs())
        )
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::adjointFarFieldPressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();

    // Act as a zeroGradient pa bc
    return tmp<Field<scalar>>
    (
        new Field<scalar>
        (
            pos(phip)*(this->patch().deltaCoeffs()*(*this))
        )
    );
}


void Foam::adjointFarFieldPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::adjointFarFieldPressureFvPatchScalarField::operator=
(
    const UList<scalar>& ul
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*ul + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    check(ptf);
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*ptf + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator+=
(
    const fvPatchField<scalar>& ptf
)
{
    check(ptf);
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this) + ptf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator-=
(
    const fvPatchField<scalar>& ptf
)
{
    check(ptf);
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this) - ptf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch() != &ptf.patch())
    {
        FatalErrorInFunction
            << "Incompatible patches for patch fields"
            << abort(FatalError);
    }

    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)*ptf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    if (&patch() != &ptf.patch())
    {
        FatalErrorInFunction
            << "Incompatible patches for patch fields"
            << abort(FatalError);
    }

    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)/ptf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator+=
(
    const Field<scalar>& tf
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this) + tf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator-=
(
    const Field<scalar>& tf
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)-tf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator*=
(
    const scalarField& tf
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)*tf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator/=
(
    const scalarField& tf
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)/tf) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator=
(
    const scalar t
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*t + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator+=
(
    const scalar t
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this) + t) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator-=
(
    const scalar t
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value
    (
        neg(phip)*((*this)-t)
      + pos(phip)*(*this)
    );

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator*=
(
    const scalar s
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)*s) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


void Foam::adjointFarFieldPressureFvPatchScalarField::operator/=
(
    const scalar s
)
{
    const fvsPatchField<scalar>& phip = boundaryContrPtr_->phib();
    scalarField value(neg(phip)*((*this)/s) + pos(phip)*(*this));

    Field<scalar>::operator=(value);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointFarFieldPressureFvPatchScalarField
    );
}

// ************************************************************************* //
