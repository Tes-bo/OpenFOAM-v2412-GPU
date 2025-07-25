/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "enthalpySorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "rhoReactionThermo.H"
#include "speciesSorptionFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::enthalpySorptionFvPatchScalarField::enthalpyModelType
>
Foam::enthalpySorptionFvPatchScalarField::enthalpyModelTypeNames
({
    { enthalpyModelType::estimated, "estimated" },
    { enthalpyModelType::calculated, "calculated" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    enthalpyModel_(enthalpyModelType::estimated),
    includeHs_(false),
    enthalpyMassLoadPtr_(nullptr),
    C_(0),
    Hvap_(0),
    speciesName_("none"),
    pName_("p"),
    TName_("T"),
    dhdt_(p.size(), 0)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    enthalpyModel_(enthalpyModelTypeNames.get("enthalpyModel", dict)),
    includeHs_(dict.getOrDefault<bool>("includeHs", true)),
    enthalpyMassLoadPtr_(nullptr),
    C_(dict.getCheckOrDefault<scalar>("C", 0, scalarMinMax::ge(0))),
    Hvap_(dict.getCheckOrDefault<scalar>("Hvap", 0, scalarMinMax::ge(0))),
    speciesName_(dict.get<word>("species")),
    pName_(dict.getOrDefault<word>("p", "p")),
    TName_(dict.getOrDefault<word>("T", "T")),
    dhdt_("dhdt", dict, p.size(), IOobjectOption::LAZY_READ)
{
    switch (enthalpyModel_)
    {
        case enthalpyModelType::calculated:
        {
            enthalpyMassLoadPtr_ =
                Function1<scalar>::New("enthalpyTable", dict, &iF.db());
            break;
        }
        case enthalpyModelType::estimated:
        {
            break;
        }
    }

    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(Zero);
    }
}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    enthalpyModel_(ptf.enthalpyModel_),
    includeHs_(ptf.includeHs_),
    enthalpyMassLoadPtr_(ptf.enthalpyMassLoadPtr_.clone()),
    C_(ptf.C_),
    Hvap_(ptf.Hvap_),
    speciesName_(ptf.speciesName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    dhdt_(ptf.dhdt_, mapper)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
    enthalpyModel_(ptf.enthalpyModel_),
    includeHs_(ptf.includeHs_),
    enthalpyMassLoadPtr_(ptf.enthalpyMassLoadPtr_.clone()),
    C_(ptf.C_),
    Hvap_(ptf.Hvap_),
    speciesName_(ptf.speciesName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    dhdt_(ptf.dhdt_)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF),
    enthalpyModel_(ptf.enthalpyModel_),
    includeHs_(ptf.includeHs_),
    enthalpyMassLoadPtr_(ptf.enthalpyMassLoadPtr_.clone()),
    C_(ptf.C_),
    Hvap_(ptf.Hvap_),
    speciesName_(ptf.speciesName_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    dhdt_(ptf.dhdt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enthalpySorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchScalarField::autoMap(m);

    dhdt_.autoMap(m);
}


void Foam::enthalpySorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    const auto& tiptf = refCast<const enthalpySorptionFvPatchScalarField>(ptf);

    dhdt_.rmap(tiptf.dhdt_, addr);
}


Foam::tmp<Foam::scalarField> Foam::enthalpySorptionFvPatchScalarField::
patchSource() const
{
    const auto& Yp =
        refCast<const speciesSorptionFvPatchScalarField>
        (
            patch().lookupPatchField<volScalarField>(speciesName_)
        );

    // mass rate [kg/sec/m3]
    tmp<scalarField> tmassb = Yp.patchSource();
    const scalarField& massb = tmassb();

    // The moles absorbed by the solid
    // dhdt[J/kg] * kg/sec/m3 = [J/m3/s]
    scalarField dhdt(dhdt_*massb);

    if (includeHs_)
    {
        const auto& pp = patch().lookupPatchField<volScalarField>(pName_);
        const auto& Tp = patch().lookupPatchField<volScalarField>(TName_);

        const auto& thermo = db().lookupObject<rhoReactionThermo>
        (
            basicThermo::dictName
        );

        const basicSpecieMixture& composition = thermo.composition();

        const label speciesId =
            thermo.composition().species()[speciesName_];

        scalarField hsp(this->patch().size(), 0);

        forAll(pp, facei)
        {
            hsp[facei] = composition.Hs(speciesId, pp[facei], Tp[facei]);
        }

        dhdt += hsp*massb;
    }

    if (debug)
    {
        Info<< " Patch enthalpy rate min/max [J/m3/sec]: "
            << gMin(dhdt) << " - " << gMax(dhdt) << endl;
    }

    return tmp<scalarField>::New(dhdt);
}


void Foam::enthalpySorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& Yp =
        refCast<const speciesSorptionFvPatchScalarField>
        (
            patch().lookupPatchField<volScalarField>(speciesName_)
        );

    switch (enthalpyModel_)
    {
        case enthalpyModelType::estimated:
        {
            dhdt_ = -C_*Hvap_;
            break;
        }
        case enthalpyModelType::calculated:
        {
            // mass [mol/kg]
            tmp<scalarField> tmassb = Yp.mass();
            const scalarField& massb = tmassb.ref();

            forAll(massb, facei)
            {
                const scalar mfacei = massb[facei];

                dhdt_[facei] = enthalpyMassLoadPtr_->value(mfacei);
            }
            break;
        }
        default:
            break;
    }

    if (debug)
    {
        Info<< "  Enthalpy change min/max [J/kg]: "
            << gMin(dhdt_) << " - " << gMax(dhdt_) << endl;
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::enthalpySorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("enthalpyModel", enthalpyModelTypeNames[enthalpyModel_]);

    if (enthalpyMassLoadPtr_)
    {
        enthalpyMassLoadPtr_->writeData(os);
    }

    os.writeEntry("species", speciesName_);

    os.writeEntryIfDifferent<bool>("includeHs", true, includeHs_);
    os.writeEntryIfDifferent<scalar>("C", scalar(0), C_);
    os.writeEntryIfDifferent<scalar>("Hvap", scalar(0), Hvap_);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);

    dhdt_.writeEntry("dhdt", os);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        enthalpySorptionFvPatchScalarField
    );
}

// ************************************************************************* //
