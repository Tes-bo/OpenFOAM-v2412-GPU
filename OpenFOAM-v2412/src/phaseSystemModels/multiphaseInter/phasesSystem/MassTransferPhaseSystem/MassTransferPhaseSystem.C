/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "MassTransferPhaseSystem.H"
#include "HashPtrTable.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "fundamentalConstants.H"

using namespace Foam::multiphaseInter;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::MassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels("massTransferModel", massTransferModels_);

    forAllConstIters(massTransferModels_, iterModel)
    {
        if (!dmdt_.found(iterModel()->pair()))
        {
            dmdt_.set
            (
                iterModel()->pair(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("dmdt",iterModel()->pair().name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        IOobject::REGISTER
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::calculateL
(
    const volScalarField& dmdtNetki,
    const phasePairKey& keyik,
    const phasePairKey& keyki,
    const volScalarField& T
) const
{
    auto tL = volScalarField::New
    (
        "tL",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero)
    );
    auto& L = tL.ref();

    if (massTransferModels_.found(keyik))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyik];

        word speciesName = interfacePtr->transferSpecie();

        const word species(speciesName.substr(0, speciesName.find('.')));

        L -= neg(dmdtNetki)*interfacePtr->L(species, T);
    }

    if (massTransferModels_.found(keyki))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyki];

        word speciesName = interfacePtr->transferSpecie();

        const word species(speciesName.substr(0, speciesName.find('.')));

        L -= pos(dmdtNetki)*interfacePtr->L(species, T);
    }

    return tL;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    auto tdmdt = volScalarField::New
    (
        "dmdt",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    auto& dmdt = tdmdt.ref();

    if (dmdt_.found(key))
    {
        dmdt = *dmdt_[key];
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    const volScalarField& T
)
{
    auto teqn = tmp<fvScalarMatrix>::New(T, dimEnergy/dimTime);
    auto& eqn = teqn.ref();

    forAllConstIters(this->phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != this->phaseModels_.end(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                // Net mass transfer from k to i phase
                auto tdmdtNetki = volScalarField::New
                (
                    "tdmdtYki",
                    IOobject::NO_REGISTER,
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                );
                auto& dmdtNetki = tdmdtNetki.ref();

                auto tSp = volScalarField::New
                (
                    "Sp",
                    IOobject::NO_REGISTER,
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime/dimTemperature, Zero)
                );
                auto& Sp = tSp.ref();

                auto tSu = volScalarField::New
                (
                    "Su",
                    IOobject::NO_REGISTER,
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                );
                auto& Su = tSu.ref();


                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    dmdtNetki -= *dmdt_[keyik];

                    tmp<volScalarField> KSp =
                        interfacePtr->KSp(interfaceCompositionModel::T, T);

                    if (KSp.valid())
                    {
                        Sp += KSp.ref();
                    }

                    tmp<volScalarField> KSu =
                        interfacePtr->KSu(interfaceCompositionModel::T, T);

                    if (KSu.valid())
                    {
                        Su += KSu.ref();
                    }

                    // If linearization is not provided used full explicit
                    if (!KSp.valid() && !KSu.valid())
                    {
                        Su += *dmdt_[keyik];
                    }
                }

                // Looking for mass transfer in the other direction (k to i)
                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    dmdtNetki += *dmdt_[keyki];


                    tmp<volScalarField> KSp =
                        interfacePtr->KSp(interfaceCompositionModel::T, T);

                    if (KSp.valid())
                    {
                        Sp += KSp.ref();
                    }

                    tmp<volScalarField> KSu =
                        interfacePtr->KSu(interfaceCompositionModel::T, T);

                    if (KSu.valid())
                    {
                        Su += KSu.ref();
                    }

                    // If linearization is not provided used full explicit
                    if (!KSp.valid() && !KSu.valid())
                    {
                        Su += *dmdt_[keyki];
                    }
                }

                tmp<volScalarField> L = calculateL(dmdtNetki, keyik, keyki, T);

                //eqn -= dmdtNetki*L;
                eqn -= fvm::Sp(Sp*L.ref(), T) + Su*L.ref();
            }
        }
    }
    return teqn;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::volTransfer
(
    const volScalarField& p
)
{
    auto teqn = tmp<fvScalarMatrix>::New(p, dimVolume/dimTime);
    auto& eqn = teqn.ref();

    auto tSp = volScalarField::New
    (
        "Sp",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimless/dimTime/dimPressure, Zero)
    );
    auto& Sp = tSp.ref();

    auto tSu = volScalarField::New
    (
        "Su",
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& Su = tSu.ref();

    forAllConstIters(this->totalPhasePairs(), iter)
    {
        const phasePair& pair = iter()();

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const phasePairKey key12
        (
            phase1.name(),
            phase2.name(),
            true
        );

        if (massTransferModels_.found(key12))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key12];

            tmp<volScalarField> KSp =
                interfacePtr->KSp(interfaceCompositionModel::P, p);

            if (KSp.valid())
            {
                Sp -=
                    KSp.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::P, p);

            if (KSu.valid())
            {
                Su -=
                    KSu.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            // If linearization is not provided used full explicit
            if (!KSp.valid() && !KSu.valid())
            {
                Su -=
                    *dmdt_[key12]
                    *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }
        }

        const phasePairKey key21
        (
            phase2.name(),
            phase1.name(),
            true
        );

        if (massTransferModels_.found(key21))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key21];

            tmp<volScalarField> KSp =
                interfacePtr->KSp(interfaceCompositionModel::P, p);

            if (KSp.valid())
            {
                Sp +=
                    KSp.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::P, p);

            if (KSu.valid())
            {
                Su +=
                    KSu.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            // If linearization is not provided used full explicit
            if (!KSp.valid() && !KSu.valid())
            {
                Su +=
                    *dmdt_[key21]
                    *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }
        }

    }

    eqn += fvm::Sp(Sp, p) + Su;
    return teqn;
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::correctMassSources
(
    const volScalarField& T
)
{
    forAllConstIters(this->phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != this->phaseModels_.end(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    tmp<volScalarField> Kexp = interfacePtr->Kexp(T);

                    *dmdt_[keyik] = Kexp.ref();

                }

                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    // Explicit temperature mass transfer rate
                    const tmp<volScalarField> Kexp = interfacePtr->Kexp(T);

                    *dmdt_[keyki] = Kexp.ref();
                }
            }
        }
    }
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::alphaTransfer
(
    SuSpTable& Su,
    SuSpTable& Sp
)
{
    // This term adds/subtracts alpha*div(U) as a source term
    // for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
    bool includeDivU(true);

    forAllConstIters(this->totalPhasePairs(), iter)
    {
        const phasePair& pair = iter()();

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& alpha1 = pair.phase1();
        const volScalarField& alpha2 = pair.phase2();

        tmp<volScalarField> tCoeffs1 = this->coeffs(phase1.name());
        const volScalarField&  coeffs1 = tCoeffs1();

        tmp<volScalarField> tCoeffs2 = this->coeffs(phase2.name());
        const volScalarField&  coeffs2 = tCoeffs2();

        // Phase 1 to phase 2
        const phasePairKey key12
        (
            phase1.name(),
            phase2.name(),
            true
        );

        tmp<volScalarField> tdmdt12(this->dmdt(key12));
        volScalarField& dmdt12 = tdmdt12.ref();

        if (massTransferModels_.found(key12))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key12];

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::alpha, phase1);

            if (KSu.valid())
            {
                dmdt12 = KSu.ref();
            }

            includeDivU = interfacePtr->includeDivU();
        }


        // Phase 2 to phase 1
        const phasePairKey key21
        (
            phase2.name(),
            phase1.name(),
            true
        );

        tmp<volScalarField> tdmdt21(this->dmdt(key21));
        volScalarField& dmdt21 = tdmdt21.ref();

        if (massTransferModels_.found(key21))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key21];

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::alpha, phase2);

            if (KSu.valid())
            {
                dmdt21 = KSu.ref();
            }

            includeDivU = interfacePtr->includeDivU();
        }

        volScalarField::Internal& SpPhase1 = Sp[phase1.name()];

        volScalarField::Internal& SuPhase1 = Su[phase1.name()];

        volScalarField::Internal& SpPhase2 = Sp[phase2.name()];

        volScalarField::Internal& SuPhase2 = Su[phase2.name()];

        const volScalarField dmdtNet(dmdt21 - dmdt12);

        const volScalarField coeffs12(coeffs1 - coeffs2);

        const surfaceScalarField& phi = this->phi();

        if (includeDivU)
        {
            SuPhase1 += fvc::div(phi)*clamp(alpha1, zero_one{});
            SuPhase2 += fvc::div(phi)*clamp(alpha2, zero_one{});
        }

        // NOTE: dmdtNet is distributed in terms =
        //  Source for phase 1 =
        //      dmdtNet/rho1
        //    - alpha1*dmdtNet(1/rho1 - 1/rho2)

        forAll(dmdtNet, celli)
        {
            scalar dmdt21 = dmdtNet[celli];
            scalar coeffs12Cell = coeffs12[celli];

            scalar alpha1Limited = clamp(alpha1[celli], zero_one{});

            // exp.
            SuPhase1[celli] += coeffs1[celli]*dmdt21;

            if (includeDivU)
            {
                if (dmdt21 > 0)
                {
                    if (coeffs12Cell > 0)
                    {
                        // imp
                        SpPhase1[celli] -= dmdt21*coeffs12Cell;
                    }
                    else if (coeffs12Cell < 0)
                    {
                        // exp
                        SuPhase1[celli] -=
                            dmdt21*coeffs12Cell*alpha1Limited;
                    }
                }
                else if (dmdt21 < 0)
                {
                    if (coeffs12Cell > 0)
                    {
                        // exp
                        SuPhase1[celli] -=
                            dmdt21*coeffs12Cell*alpha1Limited;
                    }
                    else if (coeffs12Cell < 0)
                    {
                        // imp
                        SpPhase1[celli] -= dmdt21*coeffs12Cell;
                    }
                }
            }
        }

        forAll(dmdtNet, celli)
        {
            scalar dmdt12 = -dmdtNet[celli];
            scalar coeffs21Cell = -coeffs12[celli];

            scalar alpha2Limited = clamp(alpha2[celli], zero_one{});

            // exp
            SuPhase2[celli] += coeffs2[celli]*dmdt12;

            if (includeDivU)
            {
                if (dmdt12 > 0)
                {
                    if (coeffs21Cell > 0)
                    {
                        // imp
                        SpPhase2[celli] -= dmdt12*coeffs21Cell;
                    }
                    else if (coeffs21Cell < 0)
                    {
                        // exp
                        SuPhase2[celli] -=
                            dmdt12*coeffs21Cell*alpha2Limited;
                    }
                }
                else if (dmdt12 < 0)
                {
                    if (coeffs21Cell > 0)
                    {
                        // exp
                        SuPhase2[celli] -=
                            coeffs21Cell*dmdt12*alpha2Limited;
                    }
                    else if (coeffs21Cell < 0)
                    {
                        // imp
                        SpPhase2[celli] -= dmdt12*coeffs21Cell;
                    }
                }
            }

        }

        // Update ddtAlphaMax
        this->ddtAlphaMax_ =
            max(gMax((dmdt21*coeffs1)()), gMax((dmdt12*coeffs2)()));
    }
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massSpeciesTransfer
(
    const phaseModel& phase,
    volScalarField::Internal& Su,
    volScalarField::Internal& Sp,
    const word speciesName
)
{
    // Fill the volumetric mass transfer for species
    forAllConstIters(massTransferModels_, iter)
    {
        if (iter()->transferSpecie() == speciesName)
        {
            // Explicit source
            Su +=
                  this->Su()[phase.name()]
                + this->Sp()[phase.name()]*phase.oldTime();
        }
    }
}


template<class BasePhaseSystem>
bool Foam::MassTransferPhaseSystem<BasePhaseSystem>::includeVolChange()
{
    bool includeVolChange(true);
    forAllIters(massTransferModels_, iter)
    {
        if (!iter()->includeVolChange())
        {
            includeVolChange = false;
        }
    }
    return includeVolChange;
}

// ************************************************************************* //
