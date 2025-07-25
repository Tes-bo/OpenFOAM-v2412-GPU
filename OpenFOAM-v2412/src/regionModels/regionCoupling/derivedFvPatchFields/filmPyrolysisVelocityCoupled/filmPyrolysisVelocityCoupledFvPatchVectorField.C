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

#include "filmPyrolysisVelocityCoupledFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "pyrolysisModel.H"
#include "surfaceFilmRegionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    pyrolysisRegionName_(ptf.pyrolysisRegionName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
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


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& fpvpvf
)
:
    fixedValueFvPatchVectorField(fpvpvf),
    filmRegionName_(fpvpvf.filmRegionName_),
    pyrolysisRegionName_(fpvpvf.pyrolysisRegionName_),
    phiName_(fpvpvf.phiName_),
    rhoName_(fpvpvf.rhoName_)
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& fpvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fpvpvf, iF),
    filmRegionName_(fpvpvf.filmRegionName_),
    pyrolysisRegionName_(fpvpvf.pyrolysisRegionName_),
    phiName_(fpvpvf.phiName_),
    rhoName_(fpvpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::updateCoeffs()
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

    vectorField& Up = *this;

    const label patchi = patch().index();

    // The film model
    const label filmPatchi = filmModel.regionPatchID(patchi);

    scalarField alphaFilm = filmModel.alpha().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, alphaFilm);

    vectorField UFilm = filmModel.Us().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, UFilm);

    // The pyrolysis model
    const label pyrPatchi = pyrModel.regionPatchID(patchi);

    scalarField phiPyr = pyrModel.phiGas().boundaryField()[pyrPatchi];
    pyrModel.toPrimary(pyrPatchi, phiPyr);


    const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVolume/dimTime)
    {}
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);
        phiPyr /= rhop;
    }
    else
    {
        FatalErrorInFunction
            << "Unable to process flux field phi with dimensions "
            << phi.dimensions() << nl
            << "    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    const scalarField UAvePyr(-phiPyr/patch().magSf());
    tmp<vectorField> nf(patch().nf());


    // Evaluate velocity
    Up = alphaFilm*UFilm + (1.0 - alphaFilm)*UAvePyr*nf;

    UPstream::msgType(oldTag);  // Restore tag

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
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
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmPyrolysisVelocityCoupledFvPatchVectorField
    );
}


// ************************************************************************* //
