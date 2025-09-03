/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "smoluchowskiJumpTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    psiName_("thermo:psi"),
    muName_("thermo:mu"),
    accommodationCoeff_(1.0),
    Twall_(p.size(), Zero),
    gamma_(1.4)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    muName_(ptf.muName_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    Twall_(ptf.Twall_),
    gamma_(ptf.gamma_)
{}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.getOrDefault<word>("U", "U")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    muName_(dict.getOrDefault<word>("mu", "thermo:mu")),
    accommodationCoeff_(dict.get<scalar>("accommodationCoeff")),
    Twall_("Twall", dict, p.size()),
    gamma_(dict.getOrDefault<scalar>("gamma", 1.4))
{
    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction(dict)
            << "unphysical accommodationCoeff specified"
            << "(0 < accommodationCoeff <= 1)" << endl
            << exit(FatalIOError);
    }

    if (!this->readValueEntry(dict))
    {
        // Fallback: set to the internal field
        fvPatchField<scalar>::patchInternalField(*this);
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTFvPatchScalarField::smoluchowskiJumpTFvPatchScalarField
(
    const smoluchowskiJumpTFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    accommodationCoeff_(ptpsf.accommodationCoeff_),
    Twall_(ptpsf.Twall_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::smoluchowskiJumpTFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::smoluchowskiJumpTFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void Foam::smoluchowskiJumpTFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& pmu = patch().lookupPatchField<volScalarField>(muName_);
    const auto& prho = patch().lookupPatchField<volScalarField>(rhoName_);
    const auto& ppsi = patch().lookupPatchField<volScalarField>(psiName_);
    const auto& pU = patch().lookupPatchField<volVectorField>(UName_);

    // Prandtl number reading consistent with rhoCentralFoam
    const dictionary& thermophysicalProperties =
        db().lookupObject<IOdictionary>(basicThermo::dictName);

    dimensionedScalar Pr
    (
        "Pr",
        dimless,
        thermophysicalProperties.subDict("mixture").subDict("transport")
    );

    Field<scalar> C2
    (
        pmu/prho
        *sqrt(ppsi*constant::mathematical::piByTwo)
        *2.0*gamma_/Pr.value()/(gamma_ + 1.0)
        *(2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    Field<scalar> aCoeff(prho.snGrad() - prho/C2);
    Field<scalar> KEbyRho(0.5*magSqr(pU));

    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C2));
    refValue() = Twall_;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void Foam::smoluchowskiJumpTFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    Twall_.writeEntry("Twall", os);
    os.writeEntry("gamma", gamma_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        smoluchowskiJumpTFvPatchScalarField
    );
}


// ************************************************************************* //
