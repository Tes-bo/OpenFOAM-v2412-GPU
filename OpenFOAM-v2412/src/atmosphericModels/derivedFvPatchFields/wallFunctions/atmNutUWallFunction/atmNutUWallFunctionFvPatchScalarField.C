/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "atmNutUWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> atmNutUWallFunctionFvPatchScalarField::calcNut() const
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

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const vectorField Up(Uw.patchInternalField() - Uw);
    const scalarField magUpn(mag(Up - (Up & patch().nf())*patch().nf()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar kappa = wallCoeffs_.kappa();

    const scalar t = db().time().timeOutputValue();
    const scalarField z0(z0_->value(t));

    #ifdef FULLDEBUG
    for (const scalar z : z0)
    {
        if (z < VSMALL)
        {
            FatalErrorInFunction
                << "z0 field can only contain positive values. "
                << "Please check input field z0."
                << exit(FatalError);
        }
    }
    #endif

    auto tnutw = tmp<scalarField>::New(*this);
    auto& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        const scalar Edash = (y[facei] + z0[facei])/(z0[facei] + 1e-4);
        const scalar uStar = magUpn[facei]*kappa/log(max(Edash, 1.0 + 1e-4));

        nutw[facei] = sqr(uStar)/max(magUpn[facei], 1e-6)*y[facei] - nuw[facei];
    }

    if (boundNut_)
    {
        nutw = max(nutw, scalar(0));
    }

    return tnutw;
}


void Foam::atmNutUWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<bool>("boundNut", true, boundNut_);

    if (z0_)
    {
        z0_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmNutUWallFunctionFvPatchScalarField::atmNutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutUWallFunctionFvPatchScalarField(p, iF),
    boundNut_(true),
    z0_(nullptr)
{}


atmNutUWallFunctionFvPatchScalarField::atmNutUWallFunctionFvPatchScalarField
(
    const atmNutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutUWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    boundNut_(ptf.boundNut_),
    z0_(ptf.z0_.clone(p.patch()))
{}


atmNutUWallFunctionFvPatchScalarField::atmNutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutUWallFunctionFvPatchScalarField(p, iF, dict),
    boundNut_(dict.getOrDefault<bool>("boundNut", true)),
    z0_(PatchFunction1<scalar>::New(p.patch(), "z0", dict))
{}


atmNutUWallFunctionFvPatchScalarField::atmNutUWallFunctionFvPatchScalarField
(
    const atmNutUWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutUWallFunctionFvPatchScalarField(rwfpsf),
    boundNut_(rwfpsf.boundNut_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


atmNutUWallFunctionFvPatchScalarField::atmNutUWallFunctionFvPatchScalarField
(
    const atmNutUWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutUWallFunctionFvPatchScalarField(rwfpsf, iF),
    boundNut_(rwfpsf.boundNut_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmNutUWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutUWallFunctionFvPatchScalarField::autoMap(m);

    if (z0_)
    {
        z0_->autoMap(m);
    }
}


void atmNutUWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutUWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const atmNutUWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const atmNutUWallFunctionFvPatchScalarField>(ptf);

    if (z0_)
    {
        z0_->rmap(nrwfpsf.z0_(), addr);
    }
}


void atmNutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmNutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
