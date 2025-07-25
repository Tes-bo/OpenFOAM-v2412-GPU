/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "fanFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchFieldType(scalar, fan);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::fanFvPatchField<Foam::scalar>::calcFanJump()
{
    if (!this->cyclicPatch().owner())
    {
        return;
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);

    scalarField volFlowRate(max(phip, scalar(0)));

    if (phip.internalField().dimensions() == dimVolume/dimTime)
    {
        // No conversion of volFlowRate required
    }
    else if (phip.internalField().dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);
        volFlowRate /= rhop;
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

    scalarField pdFan(patch().size(), Zero);

    switch (operatingMode_)
    {
        case operatingMode::VELOCITY:
        {
            // Note: volFlowRate now becomes face normal velocity
            volFlowRate /= patch().magSf();

            // Per-face values
            pdFan = this->jumpTable_->value(volFlowRate);

            break;
        }
        case operatingMode::UNIFORM_VELOCITY:
        {
            // Note: volFlowRate now becomes face normal velocity
            volFlowRate /= patch().magSf();

            // Set face values to patch area-averaged value
            const scalar area = gSum(patch().magSf());
            const scalar UnAve = gSum(volFlowRate*patch().magSf())/area;

            // Assign uniform value
            pdFan = this->jumpTable_->value(UnAve);

            break;
        }
        case operatingMode::VOL_FLOW_RATE:
        {
            // Face-based volFlowRate converted to patch-based volFlowRate
            // for pd curve lookup
            const scalar sumVolFlowRate = gSum(volFlowRate);

            // Assign uniform value
            pdFan = this->jumpTable_->value(sumVolFlowRate);

            break;
        }
        case operatingMode::NON_DIMENSIONAL:
        {
            // Face-based volFlowRate converted to patch-based volFlowRate
            // for pd curve lookup
            scalar sumVolFlowRate = gSum(volFlowRate);

            rpm = rpm_->value(this->db().time().timeOutputValue());
            meanDiam = dm_->value(this->db().time().timeOutputValue());

            // Create a non-dimensional flow rate
            sumVolFlowRate *=
            (
                120.0
               /stabilise
                (
                    pow3(constant::mathematical::pi*meanDiam)*rpm,
                    VSMALL
                )
            );

            const scalar pdNonDim = this->jumpTable_->value(sumVolFlowRate);

            // Convert uniform non-dimensional pdFan from curve into deltaP
            pdFan =
                pdNonDim
               *pow4(constant::mathematical::pi)*sqr(meanDiam*rpm)/1800.0;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration "
                << operatingModeNames_[operatingMode_]
                << abort(FatalError);
        }
    }


    this->setJump(pdFan);

    this->relax();
}


// ************************************************************************* //
