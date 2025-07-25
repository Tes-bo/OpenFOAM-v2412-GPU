/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "externalCoupledTemperatureMixedFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Enum.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::externalCoupledTemperatureMixedFvPatchScalarField::
    outputTemperatureType
>
Foam::externalCoupledTemperatureMixedFvPatchScalarField::outputTemperatureNames
({
    { outputTemperatureType::FLUID, "fluid" },
    { outputTemperatureType::WALL, "wall" },
});


const Foam::Enum
<
    Foam::externalCoupledTemperatureMixedFvPatchScalarField::
    refTemperatureType
>
Foam::externalCoupledTemperatureMixedFvPatchScalarField::refTemperatureNames
({
    { refTemperatureType::CELL, "cell" },
    { refTemperatureType::USER, "user" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeHeader
(
    Ostream& os
) const
{
    if (outTempType_ == outputTemperatureType::WALL)
    {
        os  << "# Values: area Twall qDot htc" << endl;
    }
    else
    {
        os  << "# Values: area Tfluid qDot htc" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(p, iF),
    outTempType_(outputTemperatureType::WALL),
    refTempType_(refTemperatureType::CELL),
    Tref_(nullptr)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    externalCoupledMixedFvPatchField<scalar>(rhs, p, iF, mapper),
    outTempType_(rhs.outTempType_),
    refTempType_(rhs.refTempType_),
    Tref_(rhs.Tref_.clone())
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //externalCoupledMixedFvPatchField<scalar>(p, iF, dict)
    externalCoupledMixedFvPatchField<scalar>(p, iF),
    outTempType_(outputTemperatureType::WALL),
    refTempType_
    (
        refTemperatureNames.getOrDefault
        (
            "htcRefTemperature",
            dict,
            refTemperatureType::CELL
        )
    ),
    Tref_(nullptr)
{
    if (dict.found("outputTemperature"))
    {
        outTempType_ = outputTemperatureNames.get("outputTemperature", dict);
    }
    else
    {
        WarningInFunction
            << "outputTemperature not specified "
            << flatOutput(outputTemperatureNames) << nl
            << "using 'wall' as compatibility default" << nl
            << endl;
    }

    if (refTempType_ == refTemperatureType::USER)
    {
        Tref_ = Function1<scalar>::New("Tref", dict, &db());
    }

    if (this->readMixedEntries(dict))
    {
        // Initialise same way as mixed
        evaluate();
    }
    else
    {
        // For convenience: initialise as fixedValue with either read value
        // or extrapolated value
        if (!this->readValueEntry(dict))
        {
            fvPatchField<scalar>::extrapolateInternal();
        }

        // Initialise as a fixed value
        this->refValue() = *this;
        this->refGrad() = Zero;
        this->valueFraction() = 1.0;
    }
}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& rhs
)
:
    externalCoupledMixedFvPatchField<scalar>(rhs),
    outTempType_(rhs.outTempType_),
    refTempType_(rhs.refTempType_),
    Tref_(rhs.Tref_.clone())
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(rhs, iF),
    outTempType_(rhs.outTempType_),
    refTempType_(rhs.refTempType_),
    Tref_(rhs.Tref_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeData
(
    Ostream& os
) const
{
    const label patchi = patch().index();

    // Heat flux [W/m2]
    scalarField qDot(this->patch().size(), Zero);

    typedef compressible::turbulenceModel cmpTurbModelType;

    static word turbName
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    static word thermoName("thermophysicalProperties");

    if (db().foundObject<cmpTurbModelType>(turbName))
    {
        const cmpTurbModelType& turbModel =
            db().lookupObject<cmpTurbModelType>(turbName);

        const basicThermo& thermo = turbModel.transport();

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = turbModel.alphaEff(patchi)*hep.snGrad();
    }
    else if (db().foundObject<basicThermo>(thermoName))
    {
        const basicThermo& thermo = db().lookupObject<basicThermo>(thermoName);

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = thermo.alpha().boundaryField()[patchi]*hep.snGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Condition requires either compressible turbulence and/or "
            << "thermo model to be available" << exit(FatalError);
    }


    // Wall temperature [K]
    const scalarField& Twall = *this;

    // Fluid temperature [K]
    scalarField Tfluid(size());

    if (refTempType_ == refTemperatureType::USER)
    {
        // User-specified reference temperature
        const scalar currTref =
            Tref_->value(this->db().time().timeOutputValue());

        Tfluid = currTref;
    }
    else
    {
        // Near wall cell temperature
        this->patchInternalField(Tfluid);
    }

    // Heat transfer coefficient [W/m2/K]
    // This htc needs to be always larger or equal to zero
    //const scalarField htc(qDot/max(Twall - Tfluid, 1e-3));
    scalarField htc(qDot.size(), Zero);
    forAll(qDot, i)
    {
        const scalar deltaT = mag(Twall[i] - Tfluid[i]);
        if (deltaT > 1e-3)
        {
            htc[i] = sign(qDot[i])*qDot[i]/deltaT;
        }
    }

    const Field<scalar>& magSf = this->patch().magSf();

    const UList<scalar>& Tout =
    (
        outTempType_ == outputTemperatureType::FLUID
      ? Tfluid
      : Twall
    );

    forAll(patch(), facei)
    {
        os  << magSf[facei] << token::SPACE
            << Tout[facei] << token::SPACE
            << qDot[facei] << token::SPACE
            << htc[facei] << nl;
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::readData
(
    Istream& is
)
{
    // Assume generic input stream so we can do line-based format and skip
    // unused columns
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line;

    forAll(*this, facei)
    {
        iss.getLine(line);
        IStringStream lineStr(line);

        lineStr
            >> this->refValue()[facei]
            >> this->refGrad()[facei]
            >> this->valueFraction()[facei];
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    externalCoupledMixedFvPatchField::write(os);
    os.writeEntry
    (
        "outputTemperature",
        outputTemperatureNames[outTempType_]
    );
    os.writeEntry
    (
        "htcRefTemperature",
        refTemperatureNames[refTempType_]
    );

    if (Tref_)
    {
        Tref_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalCoupledTemperatureMixedFvPatchScalarField
    );
}


// ************************************************************************* //
