/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "filmPyrolysisTemperatureCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "pyrolysisModel.H"
#include "surfaceFilmRegionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::
filmPyrolysisTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::
filmPyrolysisTemperatureCoupledFvPatchScalarField
(
    const filmPyrolysisTemperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    pyrolysisRegionName_(ptf.pyrolysisRegionName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::
filmPyrolysisTemperatureCoupledFvPatchScalarField
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
    pyrolysisRegionName_
    (
        dict.getOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho"))
{}


Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::
filmPyrolysisTemperatureCoupledFvPatchScalarField
(
    const filmPyrolysisTemperatureCoupledFvPatchScalarField& fptpsf
)
:
    fixedValueFvPatchScalarField(fptpsf),
    filmRegionName_(fptpsf.filmRegionName_),
    pyrolysisRegionName_(fptpsf.pyrolysisRegionName_),
    phiName_(fptpsf.phiName_),
    rhoName_(fptpsf.rhoName_)
{}


Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::
filmPyrolysisTemperatureCoupledFvPatchScalarField
(
    const filmPyrolysisTemperatureCoupledFvPatchScalarField& fptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fptpsf, iF),
    filmRegionName_(fptpsf.filmRegionName_),
    pyrolysisRegionName_(fptpsf.pyrolysisRegionName_),
    phiName_(fptpsf.phiName_),
    rhoName_(fptpsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Film model
    const auto* filmModelPtr = db().time().findObject
        <regionModels::surfaceFilmModels::surfaceFilmRegionModel>
        (filmRegionName_);

    // Pyrolysis model
    const auto* pyrModelPtr = db().time().findObject
        <regionModels::pyrolysisModels::pyrolysisModel>
        (pyrolysisRegionName_);

    if (!filmModelPtr || !pyrModelPtr)
    {
        // Do nothing on construction - film model doesn't exist yet
        return;
    }

    const auto& filmModel = *filmModelPtr;
    const auto& pyrModel = *pyrModelPtr;


    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    scalarField& Tp = *this;

    const label patchi = patch().index();

    // The film model
    const label filmPatchi = filmModel.regionPatchID(patchi);

    scalarField alphaFilm = filmModel.alpha().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, alphaFilm);

    scalarField TFilm = filmModel.Ts().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, TFilm);

    // The pyrolysis model
    const label pyrPatchi = pyrModel.regionPatchID(patchi);

    scalarField TPyr = pyrModel.T().boundaryField()[pyrPatchi];
    pyrModel.toPrimary(pyrPatchi, TPyr);


    // Evaluate temperature
    Tp = alphaFilm*TFilm + (1.0 - alphaFilm)*TPyr;

    UPstream::msgType(oldTag);  // Restore tag

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::filmPyrolysisTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntryIfDifferent<word>
    (
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        filmPyrolysisTemperatureCoupledFvPatchScalarField
    );
}


// ************************************************************************* //
