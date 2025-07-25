/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "energyJumpFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::buffered);
    }
}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf
)
:
    fixedJumpFvPatchField<scalar>(ptf)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        const basicThermo& thermo = basicThermo::lookupThermo(*this);
        label patchID = patch().index();

        const scalarField& pp = thermo.p().boundaryField()[patchID];
        const fixedJumpFvPatchScalarField& TbPatch =
            refCast<const fixedJumpFvPatchScalarField>
            (
                thermo.T().boundaryField()[patchID]
            );

        fixedJumpFvPatchScalarField& Tbp =
            const_cast<fixedJumpFvPatchScalarField&>(TbPatch);

        // force update of jump
        Tbp.evaluate(Pstream::commsTypes::buffered);

        const labelUList& faceCells = this->patch().faceCells();

        setJump
        (
            thermo.he(pp, Tbp+Tbp.jump(), faceCells)
          - thermo.he(pp, Tbp, faceCells)
        );
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::energyJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       energyJumpFvPatchScalarField
   );
}

// ************************************************************************* //
