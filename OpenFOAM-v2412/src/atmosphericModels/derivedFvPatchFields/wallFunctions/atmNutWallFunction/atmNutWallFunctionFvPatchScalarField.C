/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 CENER
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

#include "atmNutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> atmNutWallFunctionFvPatchScalarField::calcNut() const
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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    auto tnutw = tmp<scalarField>::New(*this);
    auto& nutw = tnutw.ref();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

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

    const labelList& faceCells = patch().faceCells();

    forAll(nutw, facei)
    {
        const label celli = faceCells[facei];

        // (RH:Eq. 6)
        const scalar Edash = (y[facei] + z0[facei])/(z0[facei] + z0Min_);

        // (RH:Eq. 6)
        const scalar uStarU = magUp[facei]*kappa/log(max(Edash, 1 + SMALL));

        // (RH:Eq. 7)
        const scalar uStarK = Cmu25*Foam::sqrt(k[celli]);

        // (SBJM:Eq. 7; SM:Eq. 25)
        const scalar tauw = uStarU*uStarK;

        nutw[facei] =
            max(tauw*y[facei]/(max(magUp[facei], SMALL)) - nuw[facei], 0.0);
    }

    return tnutw;
}


void Foam::atmNutWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<scalar>("z0Min", SMALL, z0Min_);

    if (z0_)
    {
        z0_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmNutWallFunctionFvPatchScalarField::atmNutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0Min_(SMALL),
    z0_(nullptr)
{}


atmNutWallFunctionFvPatchScalarField::atmNutWallFunctionFvPatchScalarField
(
    const atmNutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0Min_(ptf.z0Min_),
    z0_(ptf.z0_.clone(p.patch()))
{}


atmNutWallFunctionFvPatchScalarField::atmNutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0Min_
    (
        dict.getCheckOrDefault<scalar>
        (
            "z0Min",
            SMALL,
            scalarMinMax::ge(0)
        )
    ),
    z0_(PatchFunction1<scalar>::New(p.patch(), "z0", dict))
{}


atmNutWallFunctionFvPatchScalarField::atmNutWallFunctionFvPatchScalarField
(
    const atmNutWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0Min_(rwfpsf.z0Min_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


atmNutWallFunctionFvPatchScalarField::atmNutWallFunctionFvPatchScalarField
(
    const atmNutWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0Min_(rwfpsf.z0Min_),
    z0_(rwfpsf.z0_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmNutWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);

    if (z0_)
    {
        z0_->autoMap(m);
    }
}


void atmNutWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const auto& nrwfpsf =
        refCast<const atmNutWallFunctionFvPatchScalarField>(ptf);

    if (z0_)
    {
        z0_->rmap(nrwfpsf.z0_(), addr);
    }
}


void atmNutWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmNutWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
