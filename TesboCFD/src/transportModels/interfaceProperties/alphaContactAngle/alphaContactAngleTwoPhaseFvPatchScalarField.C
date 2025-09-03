/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::alphaContactAngleTwoPhaseFvPatchScalarField::limitControls
>
Foam::alphaContactAngleTwoPhaseFvPatchScalarField::limitControlNames_
({
    { limitControls::lcNone, "none" },
    { limitControls::lcGradient, "gradient" },
    { limitControls::lcZeroGradient, "zeroGradient" },
    { limitControls::lcAlpha, "alpha" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaContactAngleTwoPhaseFvPatchScalarField::
alphaContactAngleTwoPhaseFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    limit_(limitControls::lcZeroGradient)
{}


Foam::alphaContactAngleTwoPhaseFvPatchScalarField::
alphaContactAngleTwoPhaseFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),  // Bypass dictionary constructor
    limit_(limitControlNames_.get("limit", dict))
{
    if (this->readGradientEntry(dict))
    {
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        // Fallback: set to zero-gradient
        extrapolateInternal();
        gradient() = Zero;
    }
}


Foam::alphaContactAngleTwoPhaseFvPatchScalarField::
alphaContactAngleTwoPhaseFvPatchScalarField
(
    const alphaContactAngleTwoPhaseFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(acpsf, p, iF, mapper),
    limit_(acpsf.limit_)
{}


Foam::alphaContactAngleTwoPhaseFvPatchScalarField::
alphaContactAngleTwoPhaseFvPatchScalarField
(
    const alphaContactAngleTwoPhaseFvPatchScalarField& acpsf
)
:
    fixedGradientFvPatchScalarField(acpsf),
    limit_(acpsf.limit_)
{}


Foam::alphaContactAngleTwoPhaseFvPatchScalarField::
alphaContactAngleTwoPhaseFvPatchScalarField
(
    const alphaContactAngleTwoPhaseFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(acpsf, iF),
    limit_(acpsf.limit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaContactAngleTwoPhaseFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (limit_ == limitControls::lcGradient)
    {
        gradient() =
        patch().deltaCoeffs()
       *(
           clamp(*this + gradient()/patch().deltaCoeffs(), zero_one{})
         - *this
       );
    }
    else if (limit_ == limitControls::lcZeroGradient)
    {
        gradient() = Zero;
    }

    fixedGradientFvPatchScalarField::evaluate();

    if (limit_ == limitControls::lcAlpha)
    {
        scalarField::operator=(clamp(*this, zero_one{}));
    }
}


void Foam::alphaContactAngleTwoPhaseFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    os.writeEntry("limit", limitControlNames_[limit_]);
}


// ************************************************************************* //
