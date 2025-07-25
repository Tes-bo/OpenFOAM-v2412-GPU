/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

#include "freeSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "gravityMeshObject.H"
#include "turbulentTransportModel.H"
#include "interfaceTrackingFvMesh.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pa_(p.size(), Zero)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    pa_("pa", dict, p.size())
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(pa_);
    }
}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pa_(ptf.pa_, mapper)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    pa_(ptf.pa_)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    pa_(ptf.pa_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfacePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    pa_.autoMap(m);
}


void Foam::freeSurfacePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const freeSurfacePressureFvPatchScalarField& tiptf =
        refCast<const freeSurfacePressureFvPatchScalarField>(ptf);

    pa_.rmap(tiptf.pa_, addr);
}


void Foam::freeSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // refCast<interfaceTrackingFvMesh>
    auto* itm =
        isA_constCast<interfaceTrackingFvMesh>
        (
            patch().boundaryMesh().mesh()
        );

    if (itm)
    {
        operator==
        (
            pa_ + itm->freeSurfacePressureJump()
        );
    }
    else
    {
        // FatalError
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::freeSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    pa_.writeEntry("pa", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freeSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
