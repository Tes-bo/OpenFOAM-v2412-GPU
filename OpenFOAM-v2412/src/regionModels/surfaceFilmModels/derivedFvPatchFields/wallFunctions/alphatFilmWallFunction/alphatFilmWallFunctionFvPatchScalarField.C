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

#include "alphatFilmWallFunctionFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"
#include "surfaceFilmRegionModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFilmWallFunctionFvPatchScalarField::
alphatFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    B_(5.5),
    yPlusCrit_(11.05),
    Cmu_(0.09),
    kappa_(0.41),
    Prt_(0.85)
{}


alphatFilmWallFunctionFvPatchScalarField::
alphatFilmWallFunctionFvPatchScalarField
(
    const alphatFilmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    B_(ptf.B_),
    yPlusCrit_(ptf.yPlusCrit_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    Prt_(ptf.Prt_)
{}


alphatFilmWallFunctionFvPatchScalarField::
alphatFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    filmRegionName_
    (
        dict.getOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    B_(dict.getOrDefault("B", 5.5)),
    yPlusCrit_(dict.getOrDefault("yPlusCrit", 11.05)),
    Cmu_(dict.getOrDefault("Cmu", 0.09)),
    kappa_(dict.getOrDefault("kappa", 0.41)),
    Prt_(dict.getOrDefault("Prt", 0.85))
{}


alphatFilmWallFunctionFvPatchScalarField::
alphatFilmWallFunctionFvPatchScalarField
(
    const alphatFilmWallFunctionFvPatchScalarField& fwfpsf
)
:
    fixedValueFvPatchScalarField(fwfpsf),
    filmRegionName_(fwfpsf.filmRegionName_),
    B_(fwfpsf.B_),
    yPlusCrit_(fwfpsf.yPlusCrit_),
    Cmu_(fwfpsf.Cmu_),
    kappa_(fwfpsf.kappa_),
    Prt_(fwfpsf.Prt_)
{}


alphatFilmWallFunctionFvPatchScalarField::
alphatFilmWallFunctionFvPatchScalarField
(
    const alphatFilmWallFunctionFvPatchScalarField& fwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fwfpsf, iF),
    filmRegionName_(fwfpsf.filmRegionName_),
    B_(fwfpsf.B_),
    yPlusCrit_(fwfpsf.yPlusCrit_),
    Cmu_(fwfpsf.Cmu_),
    kappa_(fwfpsf.kappa_),
    Prt_(fwfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFilmWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto* filmModelPtr = db().time().findObject
        <regionModels::surfaceFilmModels::surfaceFilmRegionModel>
        (filmRegionName_);

    if (!filmModelPtr)
    {
        // Do nothing on construction - film model doesn't exist yet
        return;
    }

    const auto& filmModel = *filmModelPtr;


    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    const label patchi = patch().index();

    // Retrieve phase change mass from surface film model
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

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> talpha = turbModel.alpha(patchi);
    const scalarField& alphaw = talpha();

    const scalar Cmu25 = pow025(Cmu_);

    // Populate alphat field values
    scalarField& alphat = *this;
    forAll(alphat, facei)
    {
        const label faceCelli = patch().faceCells()[facei];

        const scalar uTau = Cmu25*sqrt(k[faceCelli]);

        const scalar yPlus = y[facei]*uTau/(muw[facei]/rhow[facei]);

        const scalar Pr = muw[facei]/alphaw[facei];

        scalar factor = 0;
        const scalar mStar = mDotFilmp[facei]/(y[facei]*uTau);
        if (yPlus > yPlusCrit_)
        {
            const scalar expTerm = exp(min(scalar(50), yPlusCrit_*mStar*Pr));
            const scalar yPlusRatio = yPlus/yPlusCrit_;
            const scalar powTerm = mStar*Prt_/kappa_;

            factor =
                mStar/(expTerm*(pow(yPlusRatio, powTerm)) - 1.0 + ROOTVSMALL);
        }
        else
        {
            const scalar expTerm = exp(min(scalar(50), yPlus*mStar*Pr));

            factor = mStar/(expTerm - 1.0 + ROOTVSMALL);
        }

        const scalar dx = patch().deltaCoeffs()[facei];

        const scalar alphaEff = dx*rhow[facei]*uTau*factor;

        alphat[facei] = max(alphaEff - alphaw[facei], scalar(0));
    }

    UPstream::msgType(oldTag);  // Restore tag

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFilmWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntryIfDifferent<scalar>("B", 5.5, B_);
    os.writeEntryIfDifferent<scalar>("yPlusCrit", 11.05, yPlusCrit_);
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("Prt", 0.85, Prt_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFilmWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
