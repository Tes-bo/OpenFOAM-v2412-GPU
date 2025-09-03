/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "timeVaryingAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    t0_(0.0),
    thetaT0_(0.0),
    te_(0.0),
    thetaTe_(0.0)
{}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const timeVaryingAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.te_)
{}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    t0_(dict.get<scalar>("t0")),
    thetaT0_(dict.get<scalar>("thetaT0")),
    te_(dict.get<scalar>("te")),
    thetaTe_(dict.get<scalar>("thetaTe"))
{
    evaluate();
}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const timeVaryingAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.thetaTe_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::timeVaryingAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    scalar t = patch().boundaryMesh().mesh().time().value();
    scalar theta0 = thetaT0_;

    if (t < t0_)
    {
        theta0 = thetaT0_;
    }
    else if (t > te_)
    {
        theta0 = thetaTe_;
    }
    else
    {
        theta0 = thetaT0_ + (t - t0_)*(thetaTe_ - thetaT0_)/(te_ - t0_);
    }

    return tmp<scalarField>::New(size(), theta0);
}


void Foam::timeVaryingAlphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntry("t0", t0_);
    os.writeEntry("thetaT0", thetaT0_);
    os.writeEntry("te", te_);
    os.writeEntry("thetaTe", thetaTe_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingAlphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
