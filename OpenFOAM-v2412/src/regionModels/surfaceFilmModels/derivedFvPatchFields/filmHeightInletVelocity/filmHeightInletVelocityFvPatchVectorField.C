/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "filmHeightInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltafName_("deltaf")
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltafName_(ptf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    deltafName_(dict.getOrDefault<word>("deltaf", "deltaf"))
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf
)
:
    fixedValueFvPatchVectorField(fhivpvf),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


Foam::filmHeightInletVelocityFvPatchVectorField::
filmHeightInletVelocityFvPatchVectorField
(
    const filmHeightInletVelocityFvPatchVectorField& fhivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fhivpvf, iF),
    phiName_(fhivpvf.phiName_),
    rhoName_(fhivpvf.rhoName_),
    deltafName_(fhivpvf.deltafName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);
    const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);
    const auto& deltafp = patch().lookupPatchField<volScalarField>(deltafName_);

    vectorField n(patch().nf());
    const scalarField& magSf = patch().magSf();

    operator==(n*phip/(rhop*magSf*deltafp + ROOTVSMALL));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmHeightInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("deltaf", "deltaf", deltafName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::filmHeightInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmHeightInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
