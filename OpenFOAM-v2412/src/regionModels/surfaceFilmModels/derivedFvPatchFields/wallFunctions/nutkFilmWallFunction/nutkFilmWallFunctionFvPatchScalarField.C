/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "nutkFilmWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFilmRegionModel.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    auto tuTau = tmp<scalarField>::New(patch().size(), Zero);
    auto& uTau = tuTau.ref();

    const auto* filmModelPtr = db().time().findObject
        <regionModels::surfaceFilmModels::surfaceFilmRegionModel>
        (filmRegionName_);

    if (!filmModelPtr)
    {
        // Do nothing on construction - film model doesn't exist yet
        return tuTau;
    }

    const label patchi = patch().index();

    // Retrieve phase change mass from surface film model
    const auto& filmModel = *filmModelPtr;

    const label filmPatchi = filmModel.regionPatchID(patchi);

    tmp<volScalarField> mDotFilm = filmModel.primaryMassTrans();
    scalarField mDotFilmp = mDotFilm().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, mDotFilmp);


    // Retrieve RAS turbulence model
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

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();

    forAll(uTau, facei)
    {
        const label faceCelli = patch().faceCells()[facei];

        const scalar ut = Cmu25*sqrt(k[faceCelli]);

        const scalar yPlus = y[facei]*ut/nuw[facei];

        const scalar mStar = mDotFilmp[facei]/(y[facei]*ut);

        scalar factor = 0;
        if (yPlus > yPlusCrit_)
        {
            const scalar expTerm = exp(min(scalar(50), B_*mStar));
            const scalar powTerm = pow(yPlus, mStar/kappa);
            factor = mStar/(expTerm*powTerm - 1.0 + ROOTVSMALL);
        }
        else
        {
            const scalar expTerm = exp(min(scalar(50), mStar));

            factor = mStar/(expTerm*yPlus - 1.0 + ROOTVSMALL);
        }

        uTau[facei] = sqrt(max(scalar(0), magGradU[facei]*ut*factor));
    }

    return tuTau;
}


tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::calcNut() const
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

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


void nutkFilmWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntryIfDifferent<scalar>("B", 5.5, B_);
    os.writeEntryIfDifferent<scalar>("yPlusCrit", 11.05, yPlusCrit_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    B_(5.5),
    yPlusCrit_(11.05)
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const nutkFilmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    B_(5.5),
    yPlusCrit_(11.05)
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    filmRegionName_
    (
        dict.getOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    B_(dict.getOrDefault("B", 5.5)),
    yPlusCrit_(dict.getOrDefault("yPlusCrit", 11.05))
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const nutkFilmWallFunctionFvPatchScalarField& wfpsf
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf),
    filmRegionName_(wfpsf.filmRegionName_),
    B_(wfpsf.B_),
    yPlusCrit_(wfpsf.yPlusCrit_)
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const nutkFilmWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf, iF),
    filmRegionName_(wfpsf.filmRegionName_),
    B_(wfpsf.B_),
    yPlusCrit_(wfpsf.yPlusCrit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::yPlus() const
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

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutkFilmWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutkFilmWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
