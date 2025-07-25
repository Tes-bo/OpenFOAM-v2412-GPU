/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "humidityTemperatureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::humidityTemperatureCoupledMixedFvPatchScalarField::massTransferMode
>
Foam::humidityTemperatureCoupledMixedFvPatchScalarField::massModeTypeNames_
({
    { massTransferMode::mtConstantMass, "constantMass" },
    { massTransferMode::mtCondensation, "condensation" },
    { massTransferMode::mtEvaporation, "evaporation" },
    {
        massTransferMode::mtCondensationAndEvaporation,
        "condensationAndEvaporation"
    },
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::humidityTemperatureCoupledMixedFvPatchScalarField::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


Foam::scalar
Foam::humidityTemperatureCoupledMixedFvPatchScalarField::htcCondensation
(
    const scalar Tsat,
    const scalar Re
) const
{
    // (BLID:Eq. 10.52-10.53)
    if (Tsat > 295.15 && Tsat < 373.15)
    {
        return 51104 + 2044*(Tsat - 273.15);
    }
    else
    {
        return 255510;
    }
}


Foam::volScalarField&
Foam::humidityTemperatureCoupledMixedFvPatchScalarField::thicknessField
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    volScalarField* ptr = mesh.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        ptr->store();
    }

    return *ptr;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::humidityTemperatureCoupledMixedFvPatchScalarField::
humidityTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), temperatureCoupledBase::mtFluidThermo),
    mode_(mtConstantMass),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    muName_("thermo:mu"),
    TnbrName_("T"),
    qrNbrName_("none"),
    qrName_("none"),
    specieName_("none"),
    liquid_(nullptr),
    liquidDict_(nullptr),
    mass_(patch().size(), Zero),
    Tvap_(0.0),
    myKDelta_(patch().size(), Zero),
    dmHfg_(patch().size(), Zero),
    mpCpTp_(patch().size(), Zero),
    Mcomp_(0.0),
    L_(0.0),
    fluid_(false),
    cp_(patch().size(), Zero),
    thickness_(patch().size(), Zero),
    rho_(patch().size(), Zero),
    thicknessLayers_(0),
    kappaLayers_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::humidityTemperatureCoupledMixedFvPatchScalarField::
humidityTemperatureCoupledMixedFvPatchScalarField
(
    const humidityTemperatureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    mode_(psf.mode_),
    pName_(psf.pName_),
    UName_(psf.UName_),
    rhoName_(psf.rhoName_),
    muName_(psf.muName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    specieName_(psf.specieName_),
    liquid_(psf.liquid_.clone()),
    liquidDict_(psf.liquidDict_),
    mass_(psf.mass_, mapper),
    Tvap_(psf.Tvap_),
    myKDelta_(psf.myKDelta_, mapper),
    dmHfg_(psf.dmHfg_, mapper),
    mpCpTp_(psf.mpCpTp_, mapper),
    Mcomp_(psf.Mcomp_),
    L_(psf.L_),
    fluid_(psf.fluid_),
    cp_(psf.cp_, mapper),
    thickness_(psf.thickness_, mapper),
    rho_(psf.rho_, mapper),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_)
{}


Foam::humidityTemperatureCoupledMixedFvPatchScalarField::
humidityTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(mtCondensationAndEvaporation),
    pName_(dict.getOrDefault<word>("p", "p")),
    UName_(dict.getOrDefault<word>("U", "U")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    muName_(dict.getOrDefault<word>("mu", "thermo:mu")),
    TnbrName_(dict.getOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.getOrDefault<word>("qrNbr", "none")),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    specieName_(dict.getOrDefault<word>("specie", "none")),
    liquid_(nullptr),
    liquidDict_(),
    mass_(patch().size(), Zero),
    Tvap_(0.0),
    myKDelta_(patch().size(), Zero),
    dmHfg_(patch().size(), Zero),
    mpCpTp_(patch().size(), Zero),
    Mcomp_(0.0),
    L_(0.0),
    fluid_(false),
    cp_(patch().size(), Zero),
    thickness_(patch().size(), Zero),
    rho_(patch().size(), Zero),
    thicknessLayers_(0),
    kappaLayers_(0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalIOError);
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (massModeTypeNames_.readIfPresent("mode", dict, mode_))
    {
        fluid_ = true;
    }

    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }

    if (fluid_)
    {
        switch(mode_)
        {
            case mtConstantMass:
            {
                thickness_ = scalarField("thickness", dict, p.size());
                cp_ = scalarField("cp", dict, p.size());
                rho_ = scalarField("rho", dict, p.size());

                break;
            }
            case mtCondensation:
            case mtEvaporation:
            case mtCondensationAndEvaporation:
            {
                dict.readEntry("carrierMolWeight", Mcomp_);
                dict.readEntry("L", L_);
                dict.readEntry("Tvap", Tvap_);
                liquidDict_ = dict.subDict("liquid");
                liquid_ =
                    liquidProperties::New(liquidDict_.subDict(specieName_));

                if (dict.found("thickness"))
                {
                    scalarField& Tp = *this;
                    const scalarField& magSf = patch().magSf();

                    // Assume initially standard pressure for rho calculation
                    scalar pf = 1e5;
                    thickness_ = scalarField("thickness", dict, p.size());
                    forAll(thickness_, i)
                    {
                        mass_[i] =
                            thickness_[i]*liquid_->rho(pf, Tp[i])*magSf[i];
                    }
                }
                fluid_ = true;

                break;
            }
            default:
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find mode " << mode_
                    << " on  patch " << patch().name()
                    << nl
                    << "Please set 'mode' to one of "
                    << massModeTypeNames_.sortedToc()
                    << exit(FatalIOError);
            }
        }
    }


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


Foam::humidityTemperatureCoupledMixedFvPatchScalarField::
humidityTemperatureCoupledMixedFvPatchScalarField
(
    const humidityTemperatureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    mode_(psf.mode_),
    pName_(psf.pName_),
    UName_(psf.UName_),
    rhoName_(psf.rhoName_),
    muName_(psf.muName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    specieName_(psf.specieName_),
    liquid_(psf.liquid_.clone()),
    liquidDict_(psf.liquidDict_),
    mass_(psf.mass_),
    Tvap_(psf.Tvap_),
    myKDelta_(psf.myKDelta_),
    dmHfg_(psf.dmHfg_),
    mpCpTp_(psf.mpCpTp_),
    Mcomp_(psf.Mcomp_),
    L_(psf.L_),
    fluid_(psf.fluid_),
    cp_(psf.cp_),
    thickness_(psf.thickness_),
    rho_(psf.rho_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::humidityTemperatureCoupledMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    temperatureCoupledBase::autoMap(m);

    if (fluid_)
    {
        mass_.autoMap(m);
        myKDelta_.autoMap(m);
        dmHfg_.autoMap(m);
        mpCpTp_.autoMap(m);
        cp_.autoMap(m);
        thickness_.autoMap(m);
        rho_.autoMap(m);
    }
}


void Foam::humidityTemperatureCoupledMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const humidityTemperatureCoupledMixedFvPatchScalarField& tiptf =
        refCast<const humidityTemperatureCoupledMixedFvPatchScalarField>
        (
            ptf
        );

    temperatureCoupledBase::rmap(tiptf, addr);
    if (fluid_)
    {
        mass_.rmap(tiptf.mass_, addr);
        myKDelta_.rmap(tiptf.myKDelta_, addr);
        dmHfg_.rmap(tiptf.dmHfg_, addr);
        mpCpTp_.rmap(tiptf.mpCpTp_, addr);
        cp_.rmap(tiptf.cp_, addr);
        thickness_.rmap(tiptf.thickness_, addr);
        rho_.rmap(tiptf.rho_, addr);
    }
}


void Foam::humidityTemperatureCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const scalarField& magSf = patch().magSf();

    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    const auto& nbrField = refCast
        <
            const humidityTemperatureCoupledMixedFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);


    scalarField& Tp = *this;

    const volScalarField& T =
        static_cast<const volScalarField&>(internalField());

    const scalarField TpOld(T.oldTime().boundaryField()[patch().index()]);

    scalarField Tin(patchInternalField());

    const scalarField K(this->kappa(*this));

    // Neighbour kappa done separately because we need kappa solid for the
    // htc correlation
    scalarField nbrK(nbrField.kappa(*this));
    mpp.distribute(nbrK);

    scalarField nrbDeltaCoeffs(nbrPatch.deltaCoeffs());
    mpp.distribute(nrbDeltaCoeffs);

    scalarField KDeltaNbr(nbrField.kappa(*this)*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    myKDelta_ = K*patch().deltaCoeffs();

    if (thicknessLayers_.size() > 0)
    {
        myKDelta_ = 1.0/myKDelta_;
        forAll(thicknessLayers_, iLayer)
        {
            myKDelta_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
        }
        myKDelta_ = 1.0/myKDelta_;
    }

    scalarField dm(patch().size(), Zero);

    // Fluid Side
    if (fluid_)
    {
        scalarField Yvp(patch().size(), Zero);
        const scalar dt = mesh.time().deltaTValue();

        const scalarField myDelta(patch().deltaCoeffs());

        if (mode_ != mtConstantMass)
        {
            scalarField cp(patch().size(), Zero);
            scalarField hfg(patch().size(), Zero);
            scalarField htc(patch().size(), GREAT);
            scalarField liquidRho(patch().size(), Zero);
            scalarField Tdew(patch().size(), Zero);
            scalarField RH(patch().size(), Zero);

            fixedGradientFvPatchField<scalar>& Yp =
                const_cast<fixedGradientFvPatchField<scalar>&>
                (
                    refCast
                    <
                        const fixedGradientFvPatchField<scalar>
                    >
                    (
                        patch().lookupPatchField<volScalarField>(specieName_)
                    )
                );

            const auto& pp =
                patch().lookupPatchField<volScalarField>(pName_);

            const auto& Up =
                patch().lookupPatchField<volVectorField>(UName_);

            const auto& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            const auto& mup =
                patch().lookupPatchField<volScalarField>(muName_);

            const vectorField Ui(Up.patchInternalField());
            const scalarField Yi(Yp.patchInternalField());

            forAll(Tp, faceI)
            {
                const scalar Tf = Tp[faceI];
                const scalar Tint = Tin[faceI];
                const vector Uf = Ui[faceI];
                const scalar pf = pp[faceI];

                const scalar muf = mup[faceI];
                const scalar rhof = rhop[faceI];
                const scalar nuf = muf/rhof;
                const scalar pSat = liquid_->pv(pf, Tint);
                const scalar Mv = liquid_->W();
                const scalar TSat = liquid_->pvInvert(pSat);
                const scalar Re = mag(Uf)*L_/nuf;

                cp[faceI] = liquid_->Cp(pf, Tf);
                hfg[faceI] = liquid_->hl(pf, Tf);

                // Calculate relative humidity
                const scalar invMwmean =
                        Yi[faceI]/Mv + (1.0 - Yi[faceI])/Mcomp_;
                const scalar Xv = Yi[faceI]/invMwmean/Mv;
                RH[faceI] = min(Xv*pf/pSat, 1.0);

                scalar RHmin = 0.01;
                Tdew[faceI] = 0.0;

                if (RH[faceI] > RHmin)
                {
                    scalar b = 243.5;
                    scalar c = 17.65;
                    scalar TintDeg = Tint - 273;
                    Tdew[faceI] =
                        b*(log(RH[faceI]) + (c*TintDeg)/(b + TintDeg))
                       /(c - log(RH[faceI]) - ((c*TintDeg)/(b + TintDeg)))
                       + 273;
                }

                if
                (
                    Tf < Tdew[faceI]
                 && RH[faceI] > RHmin
                 && (
                        mode_ == mtCondensation
                     || mode_ == mtCondensationAndEvaporation
                    )
                )
                {
                    htc[faceI] = htcCondensation(TSat, Re)*nbrK[faceI]/L_;

                    scalar htcTotal =
                        1.0/((1.0/myKDelta_[faceI]) + (1.0/htc[faceI]));

                    // Heat flux [W] (>0 heat is converted into mass)
                    const scalar q = (Tint - Tf)*htcTotal*magSf[faceI];

                    // Mass flux rate [Kg/s/m2]
                    dm[faceI] = q/hfg[faceI]/magSf[faceI];

                    mass_[faceI] += q/hfg[faceI]*dt;

                    // -dYp/dn = q/Dab (fixedGradient)
                    const scalar Dab = liquid_->D(pf, Tf);
                    Yvp[faceI] =
                        -min(dm[faceI]/Dab/rhof, Yi[faceI]*myDelta[faceI]);
                }
                else if
                (
                    Tf > Tvap_
                 && mass_[faceI] > 0.0
                 && (
                        mode_ == mtEvaporation
                     || mode_ == mtCondensationAndEvaporation
                    )
                )
                {
                    const scalar Dab = liquid_->D(pf, Tf);

                    const scalar Sc = nuf/Dab;
                    const scalar Sh = this->Sh(Re, Sc);

                    const scalar Ys = Mv*pSat/(Mv*pSat + Mcomp_*(pf - pSat));

                    // Mass transfer coefficient [m/s]
                    const scalar hm = Dab*Sh/L_;

                    const scalar Yinf = max(Yi[faceI], 0.0);

                    // Mass flux rate [Kg/s/m2]
                    dm[faceI] = -rhof*hm*max((Ys - Yinf), 0.0)/(1.0 - Ys);

                    // Set fixedGradient for carrier species.
                    Yvp[faceI] = -dm[faceI]/Dab/rhof;

                    // Total mass accumulated [Kg]
                    mass_[faceI] += dm[faceI]*magSf[faceI]*dt;

                    htc[faceI] = htcCondensation(TSat, Re)*nbrK[faceI]/L_;
                }
                else if (Tf > Tdew[faceI] && Tf < Tvap_ && mass_[faceI] > 0.0)
                {
                    htc[faceI] = htcCondensation(TSat, Re)*nbrK[faceI]/L_;
                }
                else if (mass_[faceI] == 0.0)
                {
                    // Do nothing
                }

                liquidRho[faceI] = liquid_->rho(pf, Tf);
            }

            mass_ = max(mass_, scalar(0));

            Yp.gradient() = Yvp;

            // Output film delta (e.g. H2OThickness) [m]
            const word fieldName(specieName_ + "Thickness");

            scalarField& pDelta =
                thicknessField
                (
                    fieldName,
                    refCast<const fvMesh>(mesh)
                ).boundaryFieldRef()[patch().index()];


            pDelta = mass_/liquidRho/magSf;

            // Weight myKDelta and htc
            myKDelta_ = 1.0/((1.0/myKDelta_) + (1.0/htc));

            mpCpTp_ = mass_*cp/dt/magSf;

            // Heat flux due to change of phase [W/m2]
            dmHfg_ = dm*hfg;

            if (debug)
            {
                // Output RH and Tdew
                scalarField& bRH =
                    thicknessField
                    (
                        "RH",
                        refCast<const fvMesh>(mesh)
                    ).boundaryFieldRef()[patch().index()];
                bRH = RH;

                scalarField& bTdew =
                    thicknessField
                    (
                        "Tdew",
                        refCast<const fvMesh>(mesh)
                    ).boundaryFieldRef()[patch().index()];
                bTdew = Tdew;
            }
        }
        else
        {
            // Inertia term [W/K/m2]
            mpCpTp_ = thickness_*rho_*cp_/dt;
        }
    }

    scalarField mpCpTpNbr(patch().size(), Zero);
    scalarField dmHfgNbr(patch().size(), Zero);

    if (!fluid_)
    {
        mpCpTpNbr = nbrField.mpCpTp();
        mpp.distribute(mpCpTpNbr);

        dmHfgNbr = nbrField.dmHfg();
        mpp.distribute(dmHfgNbr);
    }

    // Obtain Rad heat (qr)
    scalarField qr(Tp.size(), Zero);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField>(qrName_);
    }

    scalarField qrNbr(Tp.size(), Zero);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField>(qrNbrName_);
        mpp.distribute(qrNbr);
    }

    const scalarField dmHfg(dmHfgNbr + dmHfg_);

    const scalarField mpCpdt(mpCpTpNbr + mpCpTp_);

    // qr > 0 (heat up the wall)
    scalarField alpha(KDeltaNbr + mpCpdt - (qr + qrNbr)/Tp);

    valueFraction() = alpha/(alpha + myKDelta_);

    refValue() = (KDeltaNbr*nbrIntFld + mpCpdt*TpOld + dmHfg)/alpha;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug && fluid_)
    {
        scalar Qdm = gSum(dm*magSf);
        scalar QMass = gSum(mass_);
        scalar Qt = gSum(myKDelta_*(Tp - Tin)*magSf);
        scalar QtSolid = gSum(KDeltaNbr*(Tp - nbrIntFld)*magSf);

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << internalField().name() << " :" << nl
            << "    Total mass flux   [Kg/s] : " << Qdm << nl
            << "    Total mass on the wall [Kg] : " << QMass << nl
            << "    Total heat (>0 leaving the wall to the fluid) [W] : "
            << Qt << nl
            << "     Total heat (>0 leaving the wall to the solid) [W] : "
            << QtSolid << nl
            << "     wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::humidityTemperatureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);
    os.writeEntryIfDifferent<word>("Tnbr", "T", TnbrName_);
    os.writeEntryIfDifferent<word>("qrNbr", "none", qrNbrName_);
    os.writeEntryIfDifferent<word>("qr", "none", qrName_);

    if (fluid_)
    {
        os.writeEntry("mode", massModeTypeNames_[mode_]);

        os.writeEntryIfDifferent<word>("specie", "none", specieName_);

        os.writeEntry("carrierMolWeight", Mcomp_);

        os.writeEntry("L", L_);
        os.writeEntry("Tvap", Tvap_);
        os.writeEntry("fluid", fluid_);
        mass_.writeEntry("mass", os);

        if (mode_ == mtConstantMass)
        {
            cp_.writeEntry("cp", os);
            rho_.writeEntry("rho", os);
        }

        thickness_.writeEntry("thickness", os);
        word liq = "liquid";
        os << token::TAB << token::TAB << liq;
        liquidDict_.write(os);
    }

    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }

    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        humidityTemperatureCoupledMixedFvPatchScalarField
    );
}


// ************************************************************************* //
