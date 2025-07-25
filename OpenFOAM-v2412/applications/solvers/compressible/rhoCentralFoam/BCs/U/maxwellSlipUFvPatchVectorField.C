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

#include "maxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_("T"),
    rhoName_("rho"),
    psiName_("thermo:psi"),
    muName_("thermo:mu"),
    tauMCName_("tauMC"),
    accommodationCoeff_(1.0),
    Uwall_(p.size(), Zero),
    thermalCreep_(true),
    curvature_(true)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(mspvf, p, iF, mapper),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_(dict.getOrDefault<word>("T", "T")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    muName_(dict.getOrDefault<word>("mu", "thermo:mu")),
    tauMCName_(dict.getOrDefault<word>("tauMC", "tauMC")),
    accommodationCoeff_(dict.get<scalar>("accommodationCoeff")),
    Uwall_("Uwall", dict, p.size()),
    thermalCreep_(dict.getOrDefault("thermalCreep", true)),
    curvature_(dict.getOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction(dict)
            << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (this->readValueEntry(dict))
    {
        const auto* hasRefValue = dict.findEntry("refValue", keyType::LITERAL);
        const auto* hasFrac = dict.findEntry("valueFraction", keyType::LITERAL);

        if (hasRefValue && hasFrac)
        {
            this->refValue().assign(*hasRefValue, p.size());
            this->valueFraction().assign(*hasFrac, p.size());
        }
        else
        {
            this->refValue() = *this;
            this->valueFraction() = scalar(1);
        }
    }
}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(mspvf, iF),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::maxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& pmu = patch().lookupPatchField<volScalarField>(muName_);
    const auto& prho = patch().lookupPatchField<volScalarField>(rhoName_);
    const auto& ppsi = patch().lookupPatchField<volScalarField>(psiName_);

    Field<scalar> C1
    (
        sqrt(ppsi*constant::mathematical::piByTwo)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    Field<scalar> pnu(pmu/prho);
    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C1*pnu));

    refValue() = Uwall_;

    if (thermalCreep_)
    {
        const volScalarField& vsfT =
            this->db().objectRegistry::lookupObject<volScalarField>(TName_);
        label patchi = this->patch().index();
        const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];
        Field<vector> gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);
        vectorField n(patch().nf());

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    if (curvature_)
    {
        const auto& ptauMC =
            patch().lookupPatchField<volTensorField>(tauMCName_);
        vectorField n(patch().nf());

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::maxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);
    os.writeEntryIfDifferent<word>("tauMC", "tauMC", tauMCName_);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    Uwall_.writeEntry("Uwall", os);
    os.writeEntry("thermalCreep", thermalCreep_);
    os.writeEntry("curvature", curvature_);

    refValue().writeEntry("refValue", os);
    valueFraction().writeEntry("valueFraction", os);

    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        maxwellSlipUFvPatchVectorField
    );
}

// ************************************************************************* //
