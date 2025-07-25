/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd
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

#include "fanPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "TableFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fanPressureFvPatchScalarField::fanFlowDirection
>
Foam::fanPressureFvPatchScalarField::fanFlowDirectionNames_
({
    { fanFlowDirection::ffdIn, "in" },
    { fanFlowDirection::ffdOut, "out" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    fanCurve_(nullptr),
    direction_(ffdOut),
    nonDimensional_(false),
    rpm_(nullptr),
    dm_(nullptr)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(rhs, p, iF, mapper),
    fanCurve_(rhs.fanCurve_.clone()),
    direction_(rhs.direction_),
    nonDimensional_(rhs.nonDimensional_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    fanCurve_(nullptr),
    direction_(fanFlowDirectionNames_.get("direction", dict)),
    nonDimensional_(dict.getOrDefault("nonDimensional", false)),
    rpm_(nullptr),
    dm_(nullptr)
{
    // Backwards compatibility
    if (dict.found("file"))
    {
        fanCurve_.reset
        (
            new Function1Types::TableFile<scalar>("fanCurve", dict, &this->db())
        );
    }
    else
    {
        fanCurve_.reset(Function1<scalar>::New("fanCurve", dict, &this->db()));
    }

    if (nonDimensional_)
    {
        rpm_.reset(Function1<scalar>::New("rpm", dict, &this->db()));
        dm_.reset(Function1<scalar>::New("dm", dict, &this->db()));
    }
}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& rhs
)
:
    totalPressureFvPatchScalarField(rhs),
    fanCurve_(rhs.fanCurve_.clone()),
    direction_(rhs.direction_),
    nonDimensional_(rhs.nonDimensional_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(rhs, iF),
    fanCurve_(rhs.fanCurve_.clone()),
    direction_(rhs.direction_),
    nonDimensional_(rhs.nonDimensional_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName());

    const int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phip.internalField().dimensions() == dimVolume/dimTime)
    {
        volFlowRate = dir*gSum(phip);
    }
    else if (phip.internalField().dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName());
        volFlowRate = dir*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct\n"
            << "    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    // The non-dimensional parameters

    scalar rpm(0);
    scalar meanDiam(0);

    if (nonDimensional_)
    {
        rpm = rpm_->value(this->db().time().timeOutputValue());
        meanDiam = dm_->value(this->db().time().timeOutputValue());

        // Create an non-dimensional flow rate
        volFlowRate =
            120.0*volFlowRate
          / stabilise
            (
                pow3(constant::mathematical::pi * meanDiam) * rpm,
                VSMALL
            );
    }

    // Pressure drop for this flow rate
    scalar pdFan = fanCurve_->value(max(volFlowRate, scalar(0)));

    if (nonDimensional_)
    {
        // Convert the non-dimensional deltap from curve into deltaP
        pdFan =
        (
            pdFan*pow4(constant::mathematical::pi)
          * sqr(rpm * meanDiam) / 1800.0
        );
    }

    totalPressureFvPatchScalarField::updateCoeffs
    (
        p0() - dir*pdFan,
        patch().lookupPatchField<volVectorField>(UName())
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    fanCurve_->writeData(os);
    os.writeEntry("direction", fanFlowDirectionNames_[direction_]);

    if (nonDimensional_)
    {
        os.writeEntry("nonDimensional", "true");
        rpm_->writeData(os);
        dm_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureFvPatchScalarField
    );
};


// ************************************************************************* //
