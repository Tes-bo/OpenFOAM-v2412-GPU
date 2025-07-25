/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "aspectRatioModel.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"
#include "localEulerDdtScheme.H"

#include "dragModel.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::calcPhi
(
    const phaseModelList& phaseModels
) const
{
    auto tphi = surfaceScalarField::New
    (
        "phi",
        IOobject::NO_REGISTER,
        fvc::interpolate(phaseModels[0])*phaseModels[0].phi()
    );

    for (label phasei=1; phasei<phaseModels.size(); ++phasei)
    {
        tphi.ref() +=
            fvc::interpolate(phaseModels[phasei])*phaseModels[phasei].phi();
    }

    return tphi;
}


void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIters(modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),

    mesh_(mesh),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    phi_(calcPhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, Zero)
    ),

    MRF_(mesh_)
{
    // Groupings
    label movingPhasei = 0;
    label stationaryPhasei = 0;
    label anisothermalPhasei = 0;
    label multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        anisothermalPhasei += !phase.isothermal();
        multiComponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multiComponentPhaseModels_.resize(multiComponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    anisothermalPhasei = 0;
    multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            movingPhaseModels_.set(movingPhasei ++, &phase);
        }
        if (phase.stationary())
        {
            stationaryPhaseModels_.set(stationaryPhasei ++, &phase);
        }
        if (!phase.isothermal())
        {
            anisothermalPhaseModels_.set(anisothermalPhasei ++, &phase);
        }
        if (!phase.pure())
        {
            multiComponentPhaseModels_.set(multiComponentPhasei ++, &phase);
        }
    }

    // Write phi
    phi_.writeOpt(IOobject::AUTO_WRITE);

    // Blending methods
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().keyword(),
            blendingMethod::New
            (
                iter().keyword(),
                iter().dict(),
                phaseModels_.toc()
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);

    // Update motion fields
    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    auto phasei = movingPhaseModels_.cbegin();

    tmp<volScalarField> trho(phasei()*phasei().rho());

    for (++phasei; phasei != movingPhaseModels_.cend(); ++phasei)
    {
        trho.ref() += phasei()*phasei().rho();
    }

    if (stationaryPhaseModels_.empty())
    {
        return trho;
    }

    phasei = movingPhaseModels_.cbegin();

    volScalarField alpha(phasei());
    for (++phasei; phasei != movingPhaseModels_.cend(); ++phasei)
    {
        alpha += phasei();
    }

    return trho/alpha;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    auto phasei = movingPhaseModels_.cbegin();

    tmp<volVectorField> tU(phasei()*phasei().U());

    for (++phasei; phasei != movingPhaseModels_.cend(); ++phasei)
    {
        tU.ref() += phasei()*phasei().U();
    }

    if (stationaryPhaseModels_.empty())
    {
        return tU;
    }

    phasei = movingPhaseModels_.cbegin();

    volScalarField alpha(phasei());
    for (++phasei; phasei != movingPhaseModels_.cend(); ++phasei)
    {
        alpha += phasei();
    }


    return tU/alpha;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }

    return volScalarField::New
    (
        IOobject::scopedName(aspectRatioModel::typeName, "E"),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar("one", dimless, 1)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }

    return volScalarField::New
    (
        IOobject::scopedName
        (
            reactingMultiphaseEuler::surfaceTensionModel::typeName, "sigma"
        ),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar
        (
            reactingMultiphaseEuler::surfaceTensionModel::dimSigma,
            Zero
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::dmdt
(
    const phasePairKey& key
) const
{
    return volScalarField::New
    (
        IOobject::groupName("dmdt", phasePairs_[key]->name()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
}


Foam::PtrList<Foam::volScalarField> Foam::phaseSystem::dmdts() const
{
    PtrList<volScalarField> dmdts(this->phaseModels_.size());

    return dmdts;
}


void Foam::phaseSystem::solve()
{}


void Foam::phaseSystem::correct()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correct();
    }
}


void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    for (phaseModel& phase : phaseModels_)
    {
        phase.correctKinematics();

        updateDpdt = updateDpdt || phase.thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.cbegin()().thermo().p());
    }
}


void Foam::phaseSystem::correctThermo()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctThermo();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctTurbulence();
    }
}


void Foam::phaseSystem::correctEnergyTransport()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctEnergyTransport();
    }
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        for (phaseModel& phase : phaseModels_)
        {
            readOK &= phase.read();
        }

        // models ...

        return readOK;
    }

    return false;
}


Foam::tmp<Foam::volScalarField> Foam::byDt(const volScalarField& vf)
{
    if (fv::localEulerDdt::enabled(vf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaT(vf.mesh())*vf;
    }
    else
    {
        return vf/vf.mesh().time().deltaT();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::byDt(const surfaceScalarField& sf)
{
    if (fv::localEulerDdt::enabled(sf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaTf(sf.mesh())*sf;
    }
    else
    {
        return sf/sf.mesh().time().deltaT();
    }
}


// ************************************************************************* //
