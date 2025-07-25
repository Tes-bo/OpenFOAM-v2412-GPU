/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "inclinedFilmNusseltInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    filmRegionName_("surfaceFilmProperties"),
    GammaMean_(),
    a_(),
    omega_()
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    filmRegionName_(ptf.filmRegionName_),
    GammaMean_(ptf.GammaMean_.clone()),
    a_(ptf.a_.clone()),
    omega_(ptf.omega_.clone())
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
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
    GammaMean_(Function1<scalar>::New("GammaMean", dict, &db())),
    a_(Function1<scalar>::New("a", dict, &db())),
    omega_(Function1<scalar>::New("omega", dict, &db()))
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& fmfrpvf
)
:
    fixedValueFvPatchVectorField(fmfrpvf),
    filmRegionName_(fmfrpvf.filmRegionName_),
    GammaMean_(fmfrpvf.GammaMean_.clone()),
    a_(fmfrpvf.a_.clone()),
    omega_(fmfrpvf.omega_.clone())
{}


Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::
inclinedFilmNusseltInletVelocityFvPatchVectorField
(
    const inclinedFilmNusseltInletVelocityFvPatchVectorField& fmfrpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fmfrpvf, iF),
    filmRegionName_(fmfrpvf.filmRegionName_),
    GammaMean_(fmfrpvf.GammaMean_.clone()),
    a_(fmfrpvf.a_.clone()),
    omega_(fmfrpvf.omega_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve the film region from the database

    const regionModels::regionModel& region =
        db().time().lookupObject<regionModels::regionModel>(filmRegionName_);

    const regionModels::surfaceFilmModels::kinematicSingleLayer& film =
        dynamic_cast
        <
            const regionModels::surfaceFilmModels::kinematicSingleLayer&
        >(region);


    // Calculate the vector tangential to the patch
    // note: normal pointing into the domain
    const vectorField n(-patch().nf());

    const scalarField gTan(film.gTan(patchi) & n);

    if (patch().size() && (max(mag(gTan)) < SMALL))
    {
        WarningInFunction
            << "is designed to operate on patches inclined with respect to "
            << "gravity"
            << endl;
    }

    const volVectorField& nHat = film.nHat();

    const vectorField nHatp(nHat.boundaryField()[patchi].patchInternalField());

    vectorField nTan(nHatp ^ n);
    nTan /= mag(nTan) + ROOTVSMALL;


    // Calculate distance in patch tangential direction

    const vectorField& Cf = patch().Cf();
    scalarField d(nTan & Cf);


    // Calculate the wavy film height

    const scalar t = db().time().timeOutputValue();

    const scalar GMean = GammaMean_->value(t);
    const scalar a = a_->value(t);
    const scalar omega = omega_->value(t);

    const scalarField G(GMean + a*sin(omega*constant::mathematical::twoPi*d));

    const volScalarField& mu = film.mu();
    const scalarField mup(mu.boundaryField()[patchi].patchInternalField());

    const volScalarField& rho = film.rho();
    const scalarField rhop(rho.boundaryField()[patchi].patchInternalField());

    const scalarField Re(max(G, scalar(0))/mup);

    operator==(n*cbrt(gTan*mup/(3*rhop))*pow(Re, 2.0/3.0));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::inclinedFilmNusseltInletVelocityFvPatchVectorField::write
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
    GammaMean_->writeData(os);
    a_->writeData(os);
    omega_->writeData(os);

    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        inclinedFilmNusseltInletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
