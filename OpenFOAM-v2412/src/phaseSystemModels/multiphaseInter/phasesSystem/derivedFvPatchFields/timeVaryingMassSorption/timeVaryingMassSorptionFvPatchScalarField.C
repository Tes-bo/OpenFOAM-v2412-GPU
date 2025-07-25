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

#include "timeVaryingMassSorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::timeVaryingMassSorptionFvPatchScalarField::ddtSchemeType
>
Foam::timeVaryingMassSorptionFvPatchScalarField::ddtSchemeTypeNames_
({
    {
        ddtSchemeType::tsEuler,
        fv::EulerDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsCrankNicolson,
        fv::CrankNicolsonDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsBackward,
        fv::backwardDdtScheme<scalar>::typeName_()
    },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMassSorptionFvPatchScalarField::
timeVaryingMassSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    kabs_(scalar(1)),
    max_(scalar(1)),
    kdes_(scalar(1))
{}


Foam::timeVaryingMassSorptionFvPatchScalarField::
timeVaryingMassSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    kabs_(dict.getCheck<scalar>("kabs", scalarMinMax::ge(0))),
    max_(dict.getCheck<scalar>("max", scalarMinMax::ge(0))),
    kdes_(dict.getCheckOrDefault<scalar>("kdes", 0, scalarMinMax::ge(0)))
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(Zero);
    }
}


Foam::timeVaryingMassSorptionFvPatchScalarField::
timeVaryingMassSorptionFvPatchScalarField
(
    const timeVaryingMassSorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    kabs_(ptf.kabs_),
    max_(ptf.max_),
    kdes_(ptf.kdes_)
{}


Foam::timeVaryingMassSorptionFvPatchScalarField::
timeVaryingMassSorptionFvPatchScalarField
(
    const timeVaryingMassSorptionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    kabs_(ptf.kabs_),
    max_(ptf.max_),
    kdes_(ptf.kdes_)
{}


Foam::timeVaryingMassSorptionFvPatchScalarField::
timeVaryingMassSorptionFvPatchScalarField
(
    const timeVaryingMassSorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    kabs_(ptf.kabs_),
    max_(ptf.max_),
    kdes_(ptf.kdes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingMassSorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::timeVaryingMassSorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


Foam::tmp<Foam::scalarField>
Foam::timeVaryingMassSorptionFvPatchScalarField::source() const
{
    auto tsource = tmp<scalarField>::New(patch().size(), Zero);
    auto& source = tsource.ref();

    const scalarField cp(*this);
    const scalarField w(max(1 - cp/max_, scalar(0)));

    source = -kabs_*w*max(patchInternalField() - cp, scalar(0));

    source += kdes_*max(cp - patchInternalField(), scalar(0));

    return tsource;
}


void Foam::timeVaryingMassSorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const scalar dt = db().time().deltaTValue();

    const auto& fld =
        db().lookupObject<volScalarField>(this->internalField().name());
    const volScalarField& fld0 = fld.oldTime();

    // Lookup d/dt scheme from database
    const word ddtSchemeName(fld.mesh().ddtScheme(fld.name()));
    const ddtSchemeType& ddtScheme = ddtSchemeTypeNames_[ddtSchemeName];

    const scalarField cp(*this);
    const scalarField w(max(1 - cp/max_, scalar(0)));

    tmp<scalarField> dfldp =
         kabs_
        *w
        *max(patchInternalField() - cp, scalar(0))
        *dt;

    dfldp.ref() -=
        kdes_
       *max(cp - patchInternalField(), scalar(0))
       *dt;

    switch (ddtScheme)
    {
        case tsEuler:
        case tsCrankNicolson:
        {
            operator==(fld0.boundaryField()[patchi] + dfldp);

            break;
        }
        case tsBackward:
        {
            const scalar dt0 = db().time().deltaT0Value();

            const scalar c = scalar(1) + dt/(dt + dt0);
            const scalar c00 = dt*dt/(dt0*(dt + dt0));
            const scalar c0 = c + c00;

            operator==
                (
                    (
                        c0*fld0.boundaryField()[patchi]
                     - c00*fld0.oldTime().boundaryField()[patchi]
                     + dfldp
                    )/c
                );

            break;
        }
        default:
        {
            FatalErrorInFunction
                << ddtSchemeName << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingMassSorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("kabs", kabs_);
    os.writeEntry("max", max_);
    os.writeEntryIfDifferent<scalar>("kdes", scalar(0), kdes_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingMassSorptionFvPatchScalarField
    );
}

// ************************************************************************* //
