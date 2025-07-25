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

#include "speciesSorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::speciesSorptionFvPatchScalarField::equilibriumModelType
>
Foam::speciesSorptionFvPatchScalarField::equilibriumModelTypeNames
({
    { equilibriumModelType::LANGMUIR, "Langmuir" }
});


const Foam::Enum
<
    Foam::speciesSorptionFvPatchScalarField::kineticModelType
>
Foam::speciesSorptionFvPatchScalarField::kinematicModelTypeNames
({
    { kineticModelType::PseudoFirstOrder, "PseudoFirstOrder" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::speciesSorptionFvPatchScalarField::calcMoleFractions() const
{
    auto tMole = tmp<scalarField>::New(patch().size(), Zero);
    auto& Mole = tMole.ref();

    if (db().foundObject<rhoReactionThermo>(basicThermo::dictName))
    {
        const auto& thermo = db().lookupObject<rhoReactionThermo>
        (
            basicThermo::dictName
        );

        const PtrList<volScalarField>& Y = thermo.composition().Y();

        const volScalarField W(thermo.W());

        const labelUList& faceCells = patch().faceCells();

        const label speciesId =
            thermo.composition().species()[this->internalField().name()];

        const dimensionedScalar Wi
        (
            dimMass/dimMoles,
            thermo.composition().W(speciesId)
        );

        const volScalarField X(W*Y[speciesId]/Wi);

        forAll(faceCells, i)
        {
            const label cellId = faceCells[i];
            Mole[i] = X[cellId];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Thermo type is not 'rhoReactionThermo'. " << nl
            << "This BC is designed to operate with a rho based thermo."
            << exit(FatalError);
    }

    return tMole;
}


Foam::volScalarField&
Foam::speciesSorptionFvPatchScalarField::field
(
    const word& fieldName,
    const dimensionSet& dim
) const
{
    const fvMesh& mesh = this->internalField().mesh();
    auto* ptr = mesh.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            dimensionedScalar(dim, Zero)
        );

        ptr->store();
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    equilibriumModel_(equilibriumModelType::LANGMUIR),
    kinematicModel_(kineticModelType::PseudoFirstOrder),
    thicknessPtr_(nullptr),
    kabs_(1),
    kl_(0),
    max_(1),
    rhoS_(0),
    pName_("p"),
    dfldp_(p.size(), 0),
    mass_(p.size(), 0)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    equilibriumModel_(equilibriumModelTypeNames.get("equilibriumModel", dict)),
    kinematicModel_(kinematicModelTypeNames.get("kinematicModel", dict)),
    thicknessPtr_(PatchFunction1<scalar>::New(p.patch(), "thickness", dict)),
    kabs_(dict.getCheck<scalar>("kabs", scalarMinMax::ge(0))),
    kl_(dict.getCheck<scalar>("kl", scalarMinMax::ge(0))),
    max_(dict.getCheck<scalar>("max", scalarMinMax::ge(0))),
    rhoS_(dict.get<scalar>("rhoS")),
    pName_(dict.getOrDefault<word>("p", "p")),
    dfldp_("dfldp", dict, p.size(), IOobjectOption::LAZY_READ),
    mass_("mass", dict, p.size(), IOobjectOption::LAZY_READ)
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(Zero);
    }
}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    equilibriumModel_(ptf.equilibriumModel_),
    kinematicModel_(ptf.kinematicModel_),
    thicknessPtr_(ptf.thicknessPtr_.clone(patch().patch())),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    rhoS_(ptf.rhoS_),
    pName_(ptf.pName_),
    dfldp_(ptf.dfldp_, mapper),
    mass_(ptf.mass_, mapper)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
    equilibriumModel_(ptf.equilibriumModel_),
    kinematicModel_(ptf.kinematicModel_),
    thicknessPtr_(ptf.thicknessPtr_.clone(patch().patch())),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    rhoS_(ptf.rhoS_),
    pName_(ptf.pName_),
    dfldp_(ptf.dfldp_),
    mass_(ptf.mass_)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF),
    equilibriumModel_(ptf.equilibriumModel_),
    kinematicModel_(ptf.kinematicModel_),
    thicknessPtr_(ptf.thicknessPtr_.clone(patch().patch())),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    rhoS_(ptf.rhoS_),
    pName_(ptf.pName_),
    dfldp_(ptf.dfldp_),
    mass_(ptf.mass_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesSorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchScalarField::autoMap(m);

    dfldp_.autoMap(m);
    mass_.autoMap(m);

    if (thicknessPtr_)
    {
        thicknessPtr_->autoMap(m);
    }
}


void Foam::speciesSorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    const auto& tiptf = refCast<const speciesSorptionFvPatchScalarField>(ptf);

    dfldp_.rmap(tiptf.dfldp_, addr);
    mass_.rmap(tiptf.mass_, addr);

    if (thicknessPtr_)
    {
        thicknessPtr_->rmap(tiptf.thicknessPtr_(), addr);
    }
}


Foam::tmp<Foam::scalarField> Foam::speciesSorptionFvPatchScalarField::
patchSource() const
{
    const auto& thermo = db().lookupObject<rhoReactionThermo>
    (
        basicThermo::dictName
    );

    const label speciesId =
        thermo.composition().species()[this->internalField().name()];

    const scalar Wi(thermo.composition().W(speciesId));

    const scalar t = db().time().timeOutputValue();

    const scalarField h(thicknessPtr_->value(t));

    const scalarField AbyV(this->patch().magSf());

    // Solid mass [kg]
    const scalarField mass(h*AbyV*rhoS_);

    scalarField Vol(this->patch().size());

    forAll(AbyV, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        Vol[facei] = this->internalField().mesh().V()[faceCelli];
    }

    // The moles absorbed by the solid
    // dfldp[mol/kg/sec]* mass[kg]* Wi[kg/mol] / Vol[m3]= [kg/sec/m3]
    const scalarField dfldp(-dfldp_*mass*Wi*1e-3/Vol);

    if (debug)
    {
        Info<< " Patch mass rate min/max [kg/m3/sec]: "
            << gMin(dfldp) << " - " << gMax(dfldp) << endl;
    }

    return tmp<scalarField>::New(dfldp);
}


Foam::tmp<Foam::scalarField> Foam::speciesSorptionFvPatchScalarField::
mass() const
{
    return tmp<scalarField>::New(mass_);
}


void Foam::speciesSorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // equilibrium in mol/kg
    scalarField cEq(patch().size(), 0);

    switch (equilibriumModel_)
    {
        case equilibriumModelType::LANGMUIR:
        {
            // mole fraction
            tmp<scalarField> tco = calcMoleFractions();

            const auto& pp = patch().lookupPatchField<volScalarField>(pName_);

            cEq = max_*(kl_*tco()*pp/(1 + kl_*tco()*pp));
            break;
        }
        default:
            break;
    }

    // source [mol/kg/sec]
    dfldp_ = Zero;

    switch (kinematicModel_)
    {
        case kineticModelType::PseudoFirstOrder:
        {
            dfldp_ = kabs_*(cEq - mass_);
        }
        default:
            break;
    }

    // mass [mol/kg]
    const scalar dt = db().time().deltaTValue();
    mass_ += dfldp_*dt;
    mass_ = max(mass_, scalar(0));

    scalarField& pMass =
        field
        (
            "absorbedMass" + this->internalField().name(),
            dimensionSet(dimMoles/dimMass)
        ).boundaryFieldRef()[patch().index()];

    pMass = mass_;

    if (debug)
    {
        Info<< "  Absorption rate min/max [mol/kg/sec]: "
            << gMin(dfldp_) << " - " << gMax(dfldp_) << endl;
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::speciesSorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry
    (
        "equilibriumModel", equilibriumModelTypeNames[equilibriumModel_]
    );
    os.writeEntry
    (
        "kinematicModel", kinematicModelTypeNames[kinematicModel_]
    );
    if (thicknessPtr_)
    {
        thicknessPtr_->writeData(os);
    }
    os.writeEntry("kabs", kabs_);
    os.writeEntry("kl", kl_);
    os.writeEntry("max", max_);
    os.writeEntry("rhoS", rhoS_);

    dfldp_.writeEntry("dfldp", os);
    mass_.writeEntry("mass", os);
    os.writeEntryIfDifferent<word>("p", "p", pName_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        speciesSorptionFvPatchScalarField
    );
}

// ************************************************************************* //
