/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "cylindricalInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(Zero),
    axis_(Zero),
    axialVelocity_(),
    radialVelocity_(),
    rpm_()
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_.clone()),
    radialVelocity_(ptf.radialVelocity_.clone()),
    rpm_(ptf.rpm_.clone())
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    origin_(dict.getCompat<vector>("origin", {{"centre", 1712}})),
    axis_(dict.lookup("axis")),
    axialVelocity_(Function1<scalar>::New("axialVelocity", dict, &db())),
    radialVelocity_
    (
        Function1<scalar>::New("radialVelocity", dict, &db())
    ),
    rpm_(Function1<scalar>::New("rpm", dict, &db()))
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_.clone()),
    radialVelocity_(ptf.radialVelocity_.clone()),
    rpm_(ptf.rpm_.clone())
{}


Foam::cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_.clone()),
    radialVelocity_(ptf.radialVelocity_.clone()),
    rpm_(ptf.rpm_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylindricalInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const scalar axialVelocity = axialVelocity_->value(t);
    const scalar radialVelocity = radialVelocity_->value(t);
    const scalar omega = rpmToRads(rpm_->value(t));

    const vector axisHat = axis_/mag(axis_);

    const vectorField r(patch().Cf() - origin_);
    const vectorField d(r - (axisHat & r)*axisHat);

    tmp<vectorField> tangVel
    (
        (omega * axisHat) ^ d
    );

    operator==(tangVel + axisHat*axialVelocity + radialVelocity*d/mag(d));

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::cylindricalInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    axialVelocity_->writeData(os);
    radialVelocity_->writeData(os);
    rpm_->writeData(os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       cylindricalInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
