/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "phaseSystem.H"
#include "mappedPatchBase.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::compressible::
    turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::regionType
>
Foam::compressible::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::regionTypeNames_
{
    { regionType::solid, "solid" },
    { regionType::fluid, "fluid" },
};


const Foam::Enum
<
    Foam::compressible::
    turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::KMethodType
>
Foam::compressible::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::KMethodTypeNames_
{
    { KMethodType::mtSolidThermo, "solidThermo" },
    { KMethodType::mtLookup, "lookup" },
    { KMethodType::mtPhaseSystem, "phaseSystem" }
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::tmp<Foam::scalarField> Foam::compressible::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
kappa
(
    const scalarField& Tp
) const
{
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const label patchi = patch().index();

    switch (method_)
    {
        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            return thermo.kappa(patchi);
            break;
        }

        case mtLookup:
        {
            {
                const auto* ptr =
                    mesh.cfindObject<volScalarField>(kappaName_);

                if (ptr)
                {
                    return patch().patchField(*ptr);
                }
            }
            {
                const auto* ptr =
                    mesh.cfindObject<volSymmTensorField>(kappaName_);

                if (ptr)
                {
                    const symmTensorField& KWall = patch().patchField(*ptr);

                    const vectorField n(patch().nf());

                    return n & KWall & n;
                }
            }

            {
                FatalErrorInFunction
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch().name()
                    << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                    << " or volSymmTensorField."
                    << exit(FatalError);
            }
            break;
        }

        case mtPhaseSystem:
        {
            // Lookup the fluid model
            const phaseSystem& fluid =
            (
                mesh.lookupObject<phaseSystem>("phaseProperties")
            );

            auto tkappaEff = tmp<scalarField>::New(patch().size(), Zero);
            auto& kappaEff = tkappaEff.ref();

            forAll(fluid.phases(), phasei)
            {
                const phaseModel& phase = fluid.phases()[phasei];

                const fvPatchScalarField& alpha = phase.boundaryField()[patchi];

                kappaEff += alpha*phase.kappaEff(patchi)();
            }

            return tkappaEff;

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappaMethod' to one of "
                << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                << "and 'kappa' to the name of the volScalar"
                << exit(FatalError);
        }
    }

    // Return zero-sized (not nullptr)
    return tmp<scalarField>::New();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    regionType_(fluid),
    method_(mtLookup),
    kappaName_("none"),
    otherPhaseName_("vapor"),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    regionType_(psf.regionType_),
    method_(psf.method_),
    kappaName_(psf.kappaName_),
    otherPhaseName_(psf.otherPhaseName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_)
{}


turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    regionType_(regionTypeNames_.get("region", dict)),
    method_(KMethodTypeNames_.get("kappaMethod", dict)),
    kappaName_(dict.getOrDefault<word>("kappa", "none")),
    otherPhaseName_(dict.get<word>("otherPhase")),
    TnbrName_(dict.getOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.getOrDefault<word>("qrNbr", "none")),
    qrName_(dict.getOrDefault<word>("qr", "none"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }


    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    regionType_(psf.regionType_),
    method_(psf.method_),
    kappaName_(psf.kappaName_),
    otherPhaseName_(psf.otherPhaseName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const label patchi = patch().index();
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    scalarField& Tp = *this;

    const auto& nbrField =
        refCast
        <
            const turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
        >(nbrPatch.lookupPatchField<volScalarField>(TnbrName_));

    // Swap to obtain full local values of neighbour internal field
    scalarField TcNbr(nbrField.patchInternalField());
    mpp.distribute(TcNbr);


    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    KDeltaNbr = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    mpp.distribute(KDeltaNbr);

    scalarField KDelta(kappa(Tp)*patch().deltaCoeffs());

    scalarField qr(Tp.size(), 0.0);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField>(qrName_);
    }

    scalarField qrNbr(Tp.size(), 0.0);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField>(qrNbrName_);
        mpp.distribute(qrNbr);
    }


    if (regionType_ == solid)
    {
        // Lookup the fluid model in the nbrFvMesh
        const phaseSystem& fluid =
        (
            nbrMesh.lookupObject<phaseSystem>("phaseProperties")
        );

        // The BC is applied to the liquid phase of the fluid
        const phaseModel& liquid
        (
            fluid.phases()[nbrField.internalField().group()]
        );

        const phaseModel& vapor(fluid.phases()[otherPhaseName_]);


        scalarField KDeltaLiqNbr;
        const fvPatchScalarField& alphal = liquid.boundaryField()[samplePatchi];
        KDeltaLiqNbr =
            alphal*(liquid.kappaEff(samplePatchi))*nbrPatch.deltaCoeffs();
        mpp.distribute(KDeltaLiqNbr);

        scalarField KDeltaVapNbr;
        const fvPatchScalarField& alphav = vapor.boundaryField()[samplePatchi];
        KDeltaVapNbr =
            alphav*(vapor.kappaEff(samplePatchi))*nbrPatch.deltaCoeffs();
        mpp.distribute(KDeltaVapNbr);

        scalarField TvNbr;
        const fvPatchScalarField& Tv =
            vapor.thermo().T().boundaryField()[samplePatchi];
        TvNbr = Tv.patchInternalField();
        mpp.distribute(TvNbr);

        // TcNbr: liquid Tp
        // TvNbr: vapour Tp
        scalarField c(TcNbr*KDeltaLiqNbr + TvNbr*KDeltaVapNbr);

        //valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
        //refValue() = c/KDeltaNbr;
        scalarField KDeltaLiqVapNbr(KDeltaLiqNbr + KDeltaVapNbr);
        valueFraction() = KDeltaLiqVapNbr/(KDeltaLiqVapNbr + KDelta);
        refValue() = c/KDeltaLiqVapNbr;
        refGrad() = (qr + qrNbr)/kappa(Tp);

        if (debug)
        {
            scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

                Info<< "T solid : " << nl << endl;

                Info
                    << " heat transfer rate from solid:" << Q
                    << " walltemperature "
                    << " min:" << gMin(Tp)
                    << " max:" << gMax(Tp)
                    << " avg:" << gAverage(Tp) << nl
                    << endl;
        }
    }
    else if (regionType_ == fluid)
    {
        const phaseSystem& fluid =
        (
            mesh.lookupObject<phaseSystem>("phaseProperties")
        );

        const phaseModel& liquid
        (
            fluid.phases()[internalField().group()]
        );

        const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

        const fvPatchScalarField& Tv =
            vapor.thermo().T().boundaryField()[patchi];

        const fvPatchScalarField& alphav = vapor.boundaryField()[patchi];

        const scalarField KdeltaVap
        (
            alphav*(vapor.kappaEff(patchi))*patch().deltaCoeffs()
        );

        const fvPatchScalarField& alphal = liquid.boundaryField()[patchi];

        const scalarField KdeltaLiq
        (
            alphal*(liquid.kappaEff(patchi))*patch().deltaCoeffs()
        );

        // TcNbr: solid Tp
        // Tv: vapour Tp
        const scalarField c(TcNbr*KDeltaNbr + Tv.patchInternalField()*KdeltaVap);

        const scalarField a(KdeltaVap + KDeltaNbr);

        valueFraction() = a/(a + KdeltaLiq);
        refValue() = c/a;
        refGrad() = (qr + qrNbr)/kappa(Tp);

        if (debug)
        {
            scalarField Tc(patchInternalField());
            scalarField qLiq((Tp - Tc)*KdeltaLiq);
            scalarField qVap((Tp - Tv.patchInternalField())*KdeltaVap);

            Info<< "T flow : " << nl << endl;

            Info<< "  qLiq: " << gMin(qLiq) << " - " << gMax(qLiq) << endl;
            Info<< "  qVap: " << gMin(qVap) << " - " << gMax(qVap) << endl;

            scalar QLiq = gSum(qLiq*patch().magSf());
            scalar QVap = gSum(qVap*patch().magSf());

            Info<<  " Heat transfer to Liq: " << QLiq << endl;
            Info<<  " Heat transfer to Vap: " << QVap << endl;

            Info<< " walltemperature "
                << " min:" << gMin(Tp)
                << " max:" << gMax(Tp)
                << " avg:" << gAverage(Tp)
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown phase type. Valid types are: "
            << regionTypeNames_ << nl << exit(FatalError);
    }

    UPstream::msgType(oldTag);  // Restore tag

    mixedFvPatchScalarField::updateCoeffs();
}


void turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntry("kappaMethod", KMethodTypeNames_[method_]);
    os.writeEntryIfDifferent<word>("kappa","none", kappaName_);

    os.writeEntry("Tnbr", TnbrName_);

    os.writeEntryIfDifferent<word>("qrNbr", "none", qrNbrName_);
    os.writeEntryIfDifferent<word>("qr", "none", qrName_);

    os.writeEntry("region", regionTypeNames_[regionType_]);
    os.writeEntry("otherPhase", otherPhaseName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureTwoPhaseRadCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
