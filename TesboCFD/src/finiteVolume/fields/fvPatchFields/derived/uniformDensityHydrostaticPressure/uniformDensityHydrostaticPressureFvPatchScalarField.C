/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "uniformDensityHydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    rho_(0.0),
    pRefValue_(0.0),
    pRefPoint_(Zero)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    rho_(dict.get<scalar>("rho")),
    pRefValue_(dict.get<scalar>("pRefValue")),
    pRefPoint_(dict.lookup("pRefPoint"))
{
    if (!this->readValueEntry(dict))
    {
        evaluate();
    }
}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rho_(ptf.rho_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const uniformDensityHydrostaticPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    rho_(ptf.rho_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
    const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    rho_(ptf.rho_),
    pRefValue_(ptf.pRefValue_),
    pRefPoint_(ptf.pRefPoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(db().time());

    operator==
    (
        pRefValue_
      + rho_*((g.value() & patch().Cf()) - (g.value() & pRefPoint_))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("rho", rho_);
    os.writeEntry("pRefValue", pRefValue_);
    os.writeEntry("pRefPoint", pRefPoint_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformDensityHydrostaticPressureFvPatchScalarField
    );
}

// ************************************************************************* //
