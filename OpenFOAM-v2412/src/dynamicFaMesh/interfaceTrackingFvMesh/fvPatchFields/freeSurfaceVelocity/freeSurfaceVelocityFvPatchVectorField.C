/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "freeSurfaceVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interfaceTrackingFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    fvPatchFieldBase::readDict(dict);
    // Apply zero-gradient condition on start-up
    this->extrapolateInternal();
}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(fcvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfaceVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // refCast<interfaceTrackingFvMesh>
    auto* itm =
        const_cast<interfaceTrackingFvMesh*>
        (
            isA<interfaceTrackingFvMesh>(patch().boundaryMesh().mesh())
        );

    if (!itm)
    {
        // FatalError
    }
    else
    {
        gradient() = itm->freeSurfaceSnGradU();
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::freeSurfaceVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchField<vector>::write(os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freeSurfaceVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
