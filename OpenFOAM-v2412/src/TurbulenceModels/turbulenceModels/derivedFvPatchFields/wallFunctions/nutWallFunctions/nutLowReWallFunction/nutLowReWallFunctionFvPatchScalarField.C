/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "nutLowReWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::nutLowReWallFunctionFvPatchScalarField::
calcNut() const
{
    return tmp<scalarField>::New(patch().size(), Zero);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutLowReWallFunctionFvPatchScalarField::
nutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


Foam::nutLowReWallFunctionFvPatchScalarField::
nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::nutLowReWallFunctionFvPatchScalarField::
nutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


Foam::nutLowReWallFunctionFvPatchScalarField::
nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& nlrwfpsf
)
:
    nutWallFunctionFvPatchScalarField(nlrwfpsf)
{}


Foam::nutLowReWallFunctionFvPatchScalarField::
nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& nlrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(nlrwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::nutLowReWallFunctionFvPatchScalarField::
yPlus() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];

    tmp<scalarField> tnuEff = turbModel.nuEff(patchi);
    const scalarField& nuEff = tnuEff();

    return y*sqrt(nuEff*mag(Uw.snGrad()))/nuw;
}


void Foam::nutLowReWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutLowReWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
