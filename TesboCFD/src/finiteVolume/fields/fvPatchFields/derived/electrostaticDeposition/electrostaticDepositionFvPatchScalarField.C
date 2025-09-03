/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "electrostaticDepositionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::electrostaticDepositionFvPatchScalarField&
Foam::electrostaticDepositionFvPatchScalarField::eVPatch
(
    const label patchi
) const
{
    const auto& eV =
        db().lookupObject<volScalarField>(this->internalField().name());

    const volScalarField::Boundary& bf = eV.boundaryField();

    const auto& eVpf =
        refCast<const electrostaticDepositionFvPatchScalarField>(bf[patchi]);

    return const_cast<electrostaticDepositionFvPatchScalarField&>(eVpf);
}


void Foam::electrostaticDepositionFvPatchScalarField::setMaster() const
{
    if (master_ != -1)
    {
        return;
    }

    const auto& eV =
        db().lookupObject<volScalarField>(this->internalField().name());

    const volScalarField::Boundary& bf = eV.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<electrostaticDepositionFvPatchScalarField>(bf[patchi]))
        {
            electrostaticDepositionFvPatchScalarField& eVpf = eVPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            eVpf.master() = master;
        }
    }
}


void Foam::electrostaticDepositionFvPatchScalarField::round
(
    scalarField& fld,
    const scalar dcml
) const
{
    for (auto& f : fld)
    {
        f = std::round(f*dcml)/dcml;
    }
}


void Foam::electrostaticDepositionFvPatchScalarField::writeFilmFields() const
{
    const auto& eV =
        db().lookupObject<volScalarField>(this->internalField().name());

    const volScalarField::Boundary& bf = eV.boundaryField();

    const fvMesh& mesh = eV.mesh();

    volScalarField h
    (
        IOobject
        (
            IOobject::scopedName("electrostaticDeposition", "h"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimLength)
    );

    forAll(bf, patchi)
    {
        if (isA<electrostaticDepositionFvPatchScalarField>(bf[patchi]))
        {
            electrostaticDepositionFvPatchScalarField& eVpf = eVPatch(patchi);

            auto& hp = h.boundaryFieldRef()[patchi];

            hp = eVpf.h();
        }
    }

    h.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrostaticDepositionFvPatchScalarField::
electrostaticDepositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    h_(p.size(), 0),
    qcum_(p.size(), 0),
    Vfilm_(p.size(), 0),
    Ceffptr_(nullptr),
    rptr_(nullptr),
    jMin_(0),
    qMin_(0),
    Rbody_(0),
    Vi_(0),
    Vanode_(GREAT),
    phasesDict_(),
    phaseNames_(),
    phases_(),
    sigmas_(),
    sigma_(sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)), scalar(1)),
    timei_(-1),
    master_(-1)
{}


Foam::electrostaticDepositionFvPatchScalarField::
electrostaticDepositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    h_("h", dict, p.size()),
    qcum_("qCumulative", dict, p.size(), IOobjectOption::LAZY_READ),
    Vfilm_("Vfilm", dict, p.size(), IOobjectOption::LAZY_READ),
    Ceffptr_
    (
        PatchFunction1<scalar>::New(p.patch(), "CoulombicEfficiency", dict)
    ),
    rptr_(PatchFunction1<scalar>::New(p.patch(), "resistivity", dict)),
    jMin_(dict.getCheckOrDefault<scalar>("jMin", 0, scalarMinMax::ge(0))),
    qMin_(dict.getCheckOrDefault<scalar>("qMin", 0, scalarMinMax::ge(0))),
    Rbody_(dict.getCheckOrDefault<scalar>("Rbody", 0, scalarMinMax::ge(0))),
    Vi_(dict.getOrDefault<scalar>("Vi", 0)),
    Vanode_(dict.getOrDefault<scalar>("Vanode", GREAT)),
    phasesDict_(dict.subOrEmptyDict("phases")),
    phaseNames_(),
    phases_(),
    sigmas_(),
    sigma_
    (
        dimensionedScalar
        (
            sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
            dict.getCheckOrDefault<scalar>
            (
                "sigma",
                scalar(1),
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    timei_(-1),
    master_(-1)
{
    if (!this->readValueEntry(dict))
    {
        this->extrapolateInternal();  // Zero-gradient patch values
    }

    // If flow is multiphase
    if (!phasesDict_.empty())
    {
        phaseNames_.setSize(phasesDict_.size());
        phases_.setSize(phasesDict_.size());
        sigmas_.setSize(phasesDict_.size());

        label phasei = 0;
        for (const entry& dEntry : phasesDict_)
        {
            const word& key = dEntry.keyword();

            if (!dEntry.isDict())
            {
                FatalIOErrorInFunction(phasesDict_)
                    << "Entry " << key << " is not a dictionary" << nl
                    << exit(FatalIOError);
            }

            const dictionary& subDict = dEntry.dict();

            phaseNames_[phasei] = key;

            sigmas_.set
            (
                phasei,
                new dimensionedScalar
                (
                    sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
                    subDict.getCheck<scalar>
                    (
                        "sigma",
                        scalarMinMax::ge(SMALL)
                    )
                )
            );

            ++phasei;
        }

        forAll(phaseNames_, i)
        {
            phases_.set
            (
                i,
                db().getObjectPtr<volScalarField>(phaseNames_[i])
            );
        }
    }
}


Foam::electrostaticDepositionFvPatchScalarField::
electrostaticDepositionFvPatchScalarField
(
    const electrostaticDepositionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    h_(ptf.h_, mapper),
    qcum_(ptf.qcum_, mapper),
    Vfilm_(ptf.Vfilm_, mapper),
    Ceffptr_(ptf.Ceffptr_.clone(p.patch())),
    rptr_(ptf.rptr_.clone(p.patch())),
    jMin_(ptf.jMin_),
    qMin_(ptf.qMin_),
    Rbody_(ptf.Rbody_),
    Vi_(ptf.Vi_),
    Vanode_(ptf.Vanode_),
    phasesDict_(ptf.phasesDict_),
    phaseNames_(ptf.phaseNames_),
    phases_(ptf.phases_),
    sigmas_(),
    sigma_(ptf.sigma_),
    timei_(ptf.timei_),
    master_(-1)
{}


Foam::electrostaticDepositionFvPatchScalarField::
electrostaticDepositionFvPatchScalarField
(
    const electrostaticDepositionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    h_(ptf.h_),
    qcum_(ptf.qcum_),
    Vfilm_(ptf.Vfilm_),
    Ceffptr_(ptf.Ceffptr_.clone(patch().patch())),
    rptr_(ptf.rptr_.clone(patch().patch())),
    jMin_(ptf.jMin_),
    qMin_(ptf.qMin_),
    Rbody_(ptf.Rbody_),
    Vi_(ptf.Vi_),
    Vanode_(ptf.Vanode_),
    phasesDict_(ptf.phasesDict_),
    phaseNames_(ptf.phaseNames_),
    phases_(ptf.phases_),
    sigmas_(),
    sigma_(ptf.sigma_),
    timei_(ptf.timei_),
    master_(-1)
{}


Foam::electrostaticDepositionFvPatchScalarField::
electrostaticDepositionFvPatchScalarField
(
    const electrostaticDepositionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    h_(ptf.h_),
    qcum_(ptf.qcum_),
    Vfilm_(ptf.Vfilm_),
    Ceffptr_(ptf.Ceffptr_.clone(patch().patch())),
    rptr_(ptf.rptr_.clone(patch().patch())),
    jMin_(ptf.jMin_),
    qMin_(ptf.qMin_),
    Rbody_(ptf.Rbody_),
    Vi_(ptf.Vi_),
    Vanode_(ptf.Vanode_),
    phasesDict_(ptf.phasesDict_),
    phaseNames_(ptf.phaseNames_),
    phases_(ptf.phases_),
    sigmas_(),
    sigma_(ptf.sigma_),
    timei_(ptf.timei_),
    master_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrostaticDepositionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);

    h_.autoMap(m);
    qcum_.autoMap(m);
    Vfilm_.autoMap(m);

    if (Ceffptr_)
    {
        Ceffptr_->autoMap(m);
    }

    if (rptr_)
    {
        rptr_->autoMap(m);
    }
}


void Foam::electrostaticDepositionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const auto& tiptf =
        refCast<const electrostaticDepositionFvPatchScalarField>(ptf);

    h_.rmap(tiptf.h_, addr);
    qcum_.rmap(tiptf.qcum_, addr);
    Vfilm_.rmap(tiptf.Vfilm_, addr);

    if (Ceffptr_)
    {
        Ceffptr_->rmap(tiptf.Ceffptr_(), addr);
    }

    if (rptr_)
    {
        rptr_->rmap(tiptf.rptr_(), addr);
    }
}


Foam::tmp<Foam::scalarField>
Foam::electrostaticDepositionFvPatchScalarField::sigma() const
{
    const label patchi = patch().index();

    if (phases_.size())
    {
        tmp<scalarField> tsigma =
            phases_[0].boundaryField()[patchi]*sigmas_[0].value();

        for (label i = 1; i < phases_.size(); ++i)
        {
            tsigma.ref() +=
                phases_[i].boundaryField()[patchi]*sigmas_[i].value();
        }

        return tsigma;
    }

    return tmp<scalarField>::New(patch().size(), sigma_.value());
}


void Foam::electrostaticDepositionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (timei_ == db().time().timeIndex())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();
    const scalar dt = db().time().deltaTValue();
    const label patchi = patch().index();

    const auto& eV =
        db().lookupObject<volScalarField>(this->internalField().name());

    // Current density on film interface
    tmp<scalarField> tjnp = -this->sigma()*eV.boundaryField()[patchi].snGrad();
    scalarField& jnp = tjnp.ref();
    jnp = max(jnp, scalar(0)); // experimental - do not allow any negative jnp
    // experimental - avoid micro/nano currents/volts
    // to reduce snowballing effects of lateral gradients on the patch
    round(jnp);


    // Calculate film-thickness finite increments
    tmp<scalarField> tCoulombicEfficiency = Ceffptr_->value(t);
    tmp<scalarField> tdh = tCoulombicEfficiency*(jnp - jMin_)*dt;
    scalarField& dh = tdh.ref();

    // Do not allow any depletion or abrasion of deposition
    dh = max(dh, scalar(0));

    // Do not allow any deposition when accumulative specific
    // charge is less than minimum accumulative specific charge
    qcum_ += jnp*dt;

    forAll(dh, i)
    {
        if (qcum_[i] < qMin_)
        {
            dh[i] = 0;
        }
    }

    // Add finite increments of film thickness to total film thickness
    h_ += dh;


    // Calculate incremental electric potential due to film resistance
    tmp<scalarField> tresistivity = rptr_->value(t);
    tmp<scalarField> tRfilm = tresistivity*tdh;
    tmp<scalarField> tdV = jnp*tRfilm;
    Vfilm_ += tdV;
    Vfilm_ = min(Vfilm_, Vanode_);


    // Calculate electric potential due to body resistance
    tmp<scalarField> tVbody = tjnp*Rbody_;


    // Add all electric potential contributions
    operator==(min(Vi_ + Vfilm_ + tVbody, Vanode_));


    fixedValueFvPatchScalarField::updateCoeffs();

    timei_ = db().time().timeIndex();

    {
        const scalar hMin = gMin(h_);
        const scalar hMax = gMax(h_);
        const scalar hAvg = gAverage(h_);

        if (Pstream::master())
        {
            Info<< "    patch: " << patch().name()
                << ", h: min = " << hMin
                << ", max = " << hMax
                << ", average = " << hAvg << nl
                << endl;
        }
    }

    // Write here to avoid any upset to redistributePar-decompose
    if (db().time().writeTime())
    {
        // Write film thickness fields as patch fields of a volScalarField
        setMaster();

        if (patch().index() == master_)
        {
            writeFilmFields();
        }
    }
}


void Foam::electrostaticDepositionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    h_.writeEntry("h", os);

    if (Ceffptr_)
    {
        Ceffptr_->writeData(os);
    }

    if (rptr_)
    {
        rptr_->writeData(os);
    }

    if (!phasesDict_.empty())
    {
        phasesDict_.writeEntry(phasesDict_.dictName(), os);
    }
    else
    {
        sigma_.writeEntry("sigma", os);
    }

    os.writeEntryIfDifferent<scalar>("jMin", 0, jMin_);
    os.writeEntryIfDifferent<scalar>("qMin", 0, qMin_);
    os.writeEntryIfDifferent<scalar>("Rbody", 0, Rbody_);
    os.writeEntryIfDifferent<scalar>("Vi", 0, Vi_);
    os.writeEntryIfDifferent<scalar>("Vanode", GREAT, Vanode_);
    qcum_.writeEntry("qCumulative", os);
    Vfilm_.writeEntry("Vfilm", os);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electrostaticDepositionFvPatchScalarField
    );
}

// ************************************************************************* //
