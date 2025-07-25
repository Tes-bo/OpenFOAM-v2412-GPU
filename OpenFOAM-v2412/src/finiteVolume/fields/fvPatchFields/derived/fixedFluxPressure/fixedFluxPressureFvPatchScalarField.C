/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),  // Bypass dictionary constructor
    curTimeIndex_(-1)
{
    fvPatchFieldBase::readDict(dict);

    if (!this->readGradientEntry(dict) || !this->readValueEntry(dict))
    {
        extrapolateInternal();
        gradient() = Zero;
    }
}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{
    patchType() = ptf.patchType();

    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = Zero;
    gradient().map(ptf.gradient(), mapper);

    // Evaluate the value field from the gradient if the internal field is valid
    if (notNull(iF))
    {
        if (iF.size())
        {
            // Note: cannot ask for nf() if zero faces

            scalarField::operator=
            (
                //patchInternalField() + gradient()/patch().deltaCoeffs()
                // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
                // which fails for AMI patches for some mapping operations
                patchInternalField()
              + gradient()*(patch().nf() & patch().delta())
            );
        }
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value
        this->map(ptf, mapper);
    }
}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    curTimeIndex_(-1)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxPressureFvPatchScalarField::updateSnGrad
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }

    curTimeIndex_ = this->db().time().timeIndex();

    gradient() = snGradp;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }
}


void Foam::fixedFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxPressureFvPatchScalarField
    );
}


// ************************************************************************* //
