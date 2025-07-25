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

#include "alphaFixedPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaFixedPressureFvPatchScalarField::
alphaFixedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p_(p.size(), Zero)
{}


Foam::alphaFixedPressureFvPatchScalarField::
alphaFixedPressureFvPatchScalarField
(
    const alphaFixedPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p_(ptf.p_, mapper)
{}


Foam::alphaFixedPressureFvPatchScalarField::
alphaFixedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    // The 'value' is optional (handled below)
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    p_("p", dict, p.size())
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(p_);
    }
}


Foam::alphaFixedPressureFvPatchScalarField::
alphaFixedPressureFvPatchScalarField
(
    const alphaFixedPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    p_(tppsf.p_)
{}


Foam::alphaFixedPressureFvPatchScalarField::
alphaFixedPressureFvPatchScalarField
(
    const alphaFixedPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    p_(tppsf.p_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaFixedPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    p_.autoMap(m);
}


void Foam::alphaFixedPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const alphaFixedPressureFvPatchScalarField& tiptf =
        refCast<const alphaFixedPressureFvPatchScalarField>(ptf);

    p_.rmap(tiptf.p_, addr);
}


void Foam::alphaFixedPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(db().time());

    const auto& rho = patch().lookupPatchField<volScalarField>("rho");

    operator==(p_ - rho*(g.value() & patch().Cf()));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::alphaFixedPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    p_.writeEntry("p", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        alphaFixedPressureFvPatchScalarField
    );
}

// ************************************************************************* //
