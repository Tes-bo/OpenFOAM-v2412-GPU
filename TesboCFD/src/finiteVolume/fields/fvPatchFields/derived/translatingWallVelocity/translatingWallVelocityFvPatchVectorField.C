/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "translatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(nullptr)
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    U_(Function1<vector>::New("U", dict, &db()))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U_(ptf.U_.clone())
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf
)
:
    fixedValueFvPatchField<vector>(twvpvf),
    U_(twvpvf.U_.clone())
{}


Foam::translatingWallVelocityFvPatchVectorField::
translatingWallVelocityFvPatchVectorField
(
    const translatingWallVelocityFvPatchVectorField& twvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(twvpvf, iF),
    U_(twvpvf.U_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const vector U = U_->value(t);

    // Remove the component of U normal to the wall in case the wall is not flat
    const vectorField n(patch().nf());
    vectorField::operator=(U - n*(n & U));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::translatingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    U_->writeData(os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        translatingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
