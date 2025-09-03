/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "uniformTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformTotalPressureFvPatchScalarField::
uniformTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    psiName_("none"),
    gamma_(0.0),
    p0_()
{}


Foam::uniformTotalPressureFvPatchScalarField::
uniformTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    UName_(dict.getOrDefault<word>("U", "U")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? dict.get<scalar>("gamma") : 1),
    p0_(Function1<scalar>::New("p0", dict, &db()))
{
    if (!this->readValueEntry(dict))
    {
        const scalar t = this->db().time().timeOutputValue();
        fvPatchScalarField::operator==(p0_->value(t));
    }
}


Foam::uniformTotalPressureFvPatchScalarField::
uniformTotalPressureFvPatchScalarField
(
    const uniformTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(p, iF),  // Don't map
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_.clone())
{
    patchType() = ptf.patchType();

    // Set the patch pressure to the current total pressure
    // This is not ideal but avoids problems with the creation of patch faces
    const scalar t = this->db().time().timeOutputValue();
    fvPatchScalarField::operator==(p0_->value(t));
}


Foam::uniformTotalPressureFvPatchScalarField::
uniformTotalPressureFvPatchScalarField
(
    const uniformTotalPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_.clone())
{}


Foam::uniformTotalPressureFvPatchScalarField::
uniformTotalPressureFvPatchScalarField
(
    const uniformTotalPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformTotalPressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    scalar p0 = p0_->value(this->db().time().timeOutputValue());

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const auto& rho =
                patch().lookupPatchField<volScalarField>(rhoName_);

            operator==(p0 - 0.5*rho*(neg(phip))*magSqr(Up));
        }
        else
        {
            // High-speed compressible flow

            const auto& psip =
                patch().lookupPatchField<volScalarField>(psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                operator==
                (
                    p0
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(neg(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                operator==(p0/(1.0 + 0.5*psip*(neg(phip))*magSqr(Up)));
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        operator==(p0 - 0.5*(neg(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(patch().lookupPatchField<volVectorField>(UName_));
}


void Foam::uniformTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("rho", rhoName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    p0_->writeData(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
