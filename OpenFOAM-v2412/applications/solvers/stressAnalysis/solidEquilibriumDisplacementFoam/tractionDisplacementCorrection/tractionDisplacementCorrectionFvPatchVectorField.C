/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "tractionDisplacementCorrectionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), Zero),
    pressure_(p.size(), Zero)
{
    extrapolateInternal();
    gradient() = Zero;
}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const tractionDisplacementCorrectionFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    extrapolateInternal();
    gradient() = Zero;
}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const tractionDisplacementCorrectionFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionDisplacementCorrectionFvPatchVectorField::
tractionDisplacementCorrectionFvPatchVectorField
(
    const tractionDisplacementCorrectionFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementCorrectionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementCorrectionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementCorrectionFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementCorrectionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementCorrectionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& mechanicalProperties = db().lookupObject<IOdictionary>
    (
        "mechanicalProperties"
    );

    const auto& rho = patch().lookupPatchField<volScalarField>("rho");
    const auto& rhoE = patch().lookupPatchField<volScalarField>("E");
    const auto& nu = patch().lookupPatchField<volScalarField>("nu");

    scalarField E(rhoE/rho);
    scalarField mu(E/(2.0*(1.0 + nu)));
    scalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));

    if (mechanicalProperties.get<bool>("planeStress"))
    {
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
    }

    vectorField n(patch().nf());

    const auto& sigmaD = patch().lookupPatchField<volSymmTensorField>("sigmaD");
    const auto& sigmaExp = patch().lookupPatchField<volTensorField>("sigmaExp");

    gradient() =
    (
        (traction_ + pressure_*n)/rho - (n & (sigmaD + sigmaExp))
    )/(2.0*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void tractionDisplacementCorrectionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementCorrectionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
