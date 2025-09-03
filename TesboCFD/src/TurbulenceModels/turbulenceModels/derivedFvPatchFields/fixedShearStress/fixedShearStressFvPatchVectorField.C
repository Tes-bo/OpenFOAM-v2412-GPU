/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "fixedShearStressFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(Zero)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    tau0_(dict.getOrDefault<vector>("tau", Zero))
{
    this->extrapolateInternal();  // Zero-gradient patch values
}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tau0_(ptf.tau0_)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    tau0_(ptf.tau0_)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    tau0_(ptf.tau0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedShearStressFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    tmp<scalarField> nuEff(turbModel.nuEff(patch().index()));

    const vectorField Uc(patchInternalField());

    vector tauHat = tau0_/(mag(tau0_) + ROOTVSMALL);

    const scalarField& ry = patch().deltaCoeffs();

    operator==(tauHat*(tauHat & (tau0_*(1.0/(ry*nuEff)) + Uc)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedShearStressFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("tau", tau0_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedShearStressFvPatchVectorField
    );
}

// ************************************************************************* //
