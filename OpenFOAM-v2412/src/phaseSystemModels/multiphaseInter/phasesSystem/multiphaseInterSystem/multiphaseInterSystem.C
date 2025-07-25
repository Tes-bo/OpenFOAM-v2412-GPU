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

#include "multiphaseInterSystem.H"
#include "surfaceTensionModel.H"
#include "porousModel.H"

#include "HashPtrTable.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseInterSystem, 0);
}

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

const Foam::word
Foam::multiphaseInterSystem::phasePropertiesName("phaseProperties");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multiphaseInterSystem::calcMu()
{
    mu_ = mu()();
}


Foam::multiphaseInterSystem::phaseModelTable
Foam::multiphaseInterSystem::generatePhaseModels
(
    const wordList& phaseNames
) const
{
    phaseModelTable phaseModels;

    for (const word& phaseName : phaseNames)
    {
        phaseModels.insert
        (
            phaseName,
            multiphaseInter::phaseModel::New
            (
                *this,
                phaseName
            )
        );
    }

    return phaseModels;
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseInterSystem::generatePhi
(
    const phaseModelTable& phaseModels
) const
{
    auto iter = phaseModels.cbegin();

    auto tmpPhi = surfaceScalarField::New
    (
        "phi",
        IOobject::NO_REGISTER,
        fvc::interpolate(iter()()) * iter()->phi()
    );

    for (++iter; iter != phaseModels.cend(); ++iter)
    {
        tmpPhi.ref() += fvc::interpolate(iter()()) * iter()->phi();
    }

    return tmpPhi;
}


void Foam::multiphaseInterSystem::generatePairs(const dictTable& modelDicts)
{
    forAllConstIters(modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {
            // do nothing ...
        }
        else if (key.ordered())
        {
            // New ordered pair
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
        else
        {
            // New unordered pair
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


void Foam::multiphaseInterSystem::generatePairsTable()
{
    forAllConstIters(phaseModels_, phaseIter1)
    {
        forAllConstIters(phaseModels_, phaseIter2)
        {
            if (phaseIter1()->name() != phaseIter2()->name())
            {
                phasePairKey key
                (
                    phaseIter1()->name(),
                    phaseIter2()->name(),
                    true
                );

                phasePairKey keyInverse
                (
                    phaseIter2()->name(),
                    phaseIter1()->name(),
                    true
                );

                if
                (
                    !totalPhasePairs_.found(key)
                 && !totalPhasePairs_.found(keyInverse)
                )
                {
                    totalPhasePairs_.set
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
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInterSystem::multiphaseInterSystem
(
    const fvMesh& mesh
)
:
    basicThermo(mesh, word::null, phasePropertiesName),
    mesh_(mesh),
    mu_
    (
        IOobject
        (
            "mu",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimViscosity*dimDensity, Zero)
    ),
    phaseNames_(get<wordList>("phases")),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVolume/dimTime, Zero)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),
    phaseModels_(generatePhaseModels(phaseNames_)),
    phasePairs_(),
    totalPhasePairs_(),
    Prt_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Prt", *this, 1.0
        )
    )
{
    rhoPhi_.setOriented();
    phi_.setOriented();

    // sub models
    if (found("surfaceTension"))
    {
        generatePairsAndSubModels
        (
            "surfaceTension",
            surfaceTensionModels_
        );
    }
    if (found("interfacePorous"))
    {
        generatePairsAndSubModels
        (
            "interfacePorous",
            mesh_,
            interfacePorousModelTable_
        );
    }

    // Total phase pair
    generatePairsTable();

    // Update mu_
    calcMu();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseInterSystem::~multiphaseInterSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::hc() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tAlphaHc
    (
        iter()() * iter()->hc()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tAlphaHc.ref() += iter()() * iter()->hc();
    }

    return tAlphaHc;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::rho() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpRho
    (
        iter()() * iter()->rho()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpRho.ref() += iter()() * iter()->rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::rho
(
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpRho
    (
        iter()().boundaryField()[patchI]
      * iter()->rho()().boundaryField()[patchI]
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpRho.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->rho()().boundaryField()[patchI]
        );
    }

    return tmpRho;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::Cp() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCp
    (
        iter()() * iter()->Cp()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp();
    }

    return tmpCp;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCp
    (
        iter()() * iter()->Cp(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp(p, T, patchI);
    }

    return tmpCp;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::Cv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCv
    (
        iter()() * iter()->Cv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCv.ref() += iter()() * iter()->Cv();
    }

    return tmpCv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCv
    (
        iter()() * iter()->Cv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCv.ref() += iter()() * iter()->Cv(p, T, patchI);
    }

    return tmpCv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::rhoEoS
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::gamma() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCp
    (
        iter()() * iter()->Cp()
    );

    tmp<volScalarField> tmpCv
    (
        iter()() * iter()->Cv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCp.ref() += iter()() * iter()->Cp();
        tmpCv.ref() += iter()() * iter()->Cv();
    }

    return (tmpCp/tmpCv);
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return
    (
        gamma()().boundaryField()[patchI]
    );
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::Cpv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCpv
    (
        iter()() * iter()->Cpv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() += iter()() * iter()->Cpv();
    }

    return tmpCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCpv
    (
        iter()() * iter()->Cpv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() += iter()() * iter()->Cpv(p, T, patchI);
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::CpByCpv() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpCpByCpv
    (
        iter()() * iter()->CpByCpv()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpByCpv.ref() += iter()() * iter()->CpByCpv();
    }

    return tmpCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpCpv
    (
        iter()().boundaryField()[patchI]
      * iter()->CpByCpv(p, T, patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpCpv.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->CpByCpv(p, T, patchI)
        );
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::W() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::kappa() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpkappa
    (
        iter()() * iter()->kappa()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpkappa.ref() += iter()() * iter()->kappa();
    }

    return tmpkappa;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::kappa
(
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpKappa
    (
        iter()().boundaryField()[patchI]
      * iter()->kappa(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpKappa.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->kappa(patchI)
        );
    }

    return tmpKappa;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::alphahe() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> talphaEff
    (
        phaseModelIter()()*phaseModelIter()->alphahe()
    );

    for (; phaseModelIter != phaseModels_.end(); ++phaseModelIter)
    {
        talphaEff.ref() += phaseModelIter()()*phaseModelIter()->alphahe();
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::alphahe
(
    const label patchi
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> talphaEff
    (
        phaseModelIter()().boundaryField()[patchi]
       *phaseModelIter()->alphahe(patchi)
    );

    for (; phaseModelIter != phaseModels_.end(); ++phaseModelIter)
    {
        talphaEff.ref() +=
            phaseModelIter()().boundaryField()[patchi]
           *phaseModelIter()->alphahe(patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField>Foam::multiphaseInterSystem::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::kappaEff
(
    const scalarField& kappat,
    const label patchI
) const
{
    return kappa(patchI) + kappat;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpAlpha
    (
        iter()() * iter()->alpha()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpAlpha.ref() += iter()() * iter()->alpha();
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::alphaEff
(
    const scalarField& alphat,
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpAlpha
    (
        iter()().boundaryField()[patchI]
      * iter()->alpha(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpAlpha.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->alpha(patchI)
        );
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


const Foam::dimensionedScalar& Foam::multiphaseInterSystem::Prt() const
{
    return Prt_;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::mu() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpMu
    (
        iter()() * iter()->mu()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpMu.ref() += iter()() * iter()->mu();
    }

    return tmpMu;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::mu
(
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpMu
    (
        iter()().boundaryField()[patchI]
      * iter()->mu(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpMu.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->mu(patchI)
        );
    }

    return tmpMu;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::nu() const
{
    auto iter = phaseModels_.cbegin();

    tmp<volScalarField> tmpNu
    (
        iter()() * iter()->nu()
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpNu.ref() += iter()() * iter()->nu();
    }

    return tmpNu;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseInterSystem::nu
(
    const label patchI
) const
{
    auto iter = phaseModels_.cbegin();

    tmp<scalarField> tmpNu
    (
        iter()().boundaryField()[patchI]
      * iter()->nu(patchI)
    );

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        tmpNu.ref() +=
        (
            iter()().boundaryField()[patchI]
          * iter()->nu(patchI)
        );
    }

    return tmpNu;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::mut() const
{
    return turb_->mut();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::muEff() const
{
    return turb_->muEff();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::nut() const
{
    return turb_->nut();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::nuEff() const
{
    return turb_->nuEff();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::kappaEff() const
{
    return
    (
        this->kappa() + this->Cp()*turb_->mut()/Prt_
    );
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInterSystem::kappaEff(const label patchi) const
{
    const scalarField Cp(this->Cp()().boundaryField()[patchi]);
    const scalarField kappaEffp
    (
         this->kappa(patchi) + Cp*turb_->mut(patchi)/Prt_.value()
    );

    return tmp<scalarField>::New(kappaEffp);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::alphaEff() const
{
    return this->alpha() + turb_->mut()/Prt_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseInterSystem::alphaEff(const label patchi) const
{
    return (this->alpha(patchi) + turb_->mut(patchi))/Prt_.value();
}


const Foam::surfaceScalarField& Foam::multiphaseInterSystem::phi() const
{
    return phi_;
}


Foam::surfaceScalarField& Foam::multiphaseInterSystem::phi()
{
    return phi_;
}


const Foam::surfaceScalarField& Foam::multiphaseInterSystem::rhoPhi() const
{
    return rhoPhi_;
}


Foam::surfaceScalarField& Foam::multiphaseInterSystem::rhoPhi()
{
    return rhoPhi_;
}


void Foam::multiphaseInterSystem::correct()
{
    forAllIters(phaseModels_, iter)
    {
        iter()->correct();
    }

    calcMu();
}


void Foam::multiphaseInterSystem::correctTurbulence()
{
    forAllIters(phaseModels_, iter)
    {
        iter()->correctTurbulence();
    }
}


const Foam::multiphaseInterSystem::phaseModelTable&
Foam::multiphaseInterSystem::phases() const
{
    return phaseModels_;
}


Foam::multiphaseInterSystem::phaseModelTable&
Foam::multiphaseInterSystem::phases()
{
    return phaseModels_;
}


const Foam::multiphaseInterSystem::phasePairTable&
Foam::multiphaseInterSystem::totalPhasePairs() const
{
    return totalPhasePairs_;
}


Foam::multiphaseInterSystem::phasePairTable&
Foam::multiphaseInterSystem::totalPhasePairs()
{
    return totalPhasePairs_;
}


bool Foam::multiphaseInterSystem::incompressible() const
{
    forAllConstIters(phaseModels_, iter)
    {
        if (!iter()->thermo().incompressible())
        {
            return false;
        }
    }

    return true;
}


bool Foam::multiphaseInterSystem::incompressible(const word phaseName) const
{
    return phaseModels_[phaseName]->thermo().incompressible();
}


bool Foam::multiphaseInterSystem::isochoric() const
{
    forAllConstIters(phaseModels_, iter)
    {
        if (!iter()->thermo().isochoric())
        {
            return false;
        }
    }

    return true;
}


const Foam::fvMesh& Foam::multiphaseInterSystem::mesh() const
{
    return mesh_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseInterSystem::surfaceTensionForce() const
{
    auto tstf = surfaceScalarField::New
    (
        "surfaceTensionForce",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar({1, -2, -2, 0, 0, 0}, Zero)
    );

    auto& stf = tstf.ref();
    stf.setOriented();

    if (surfaceTensionModels_.size())
    {
        forAllConstIters(phaseModels_, iter1)
        {
            const volScalarField& alpha1 = iter1()();

            auto iter2 = iter1;

            for (++iter2; iter2 != phaseModels_.cend(); ++iter2)
            {
                const volScalarField& alpha2 = iter2()();

                stf +=
                    fvc::interpolate
                    (
                        surfaceTensionCoeff
                        (
                            phasePairKey(iter1()->name(), iter2()->name())
                        )
                    )
                * fvc::interpolate(K(alpha1, alpha2))*
                    (
                        fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                      - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                    );
            }
        }
    }

    return tstf;
}


Foam::tmp<Foam::volVectorField> Foam::multiphaseInterSystem::U() const
{
    auto tstf = volVectorField::New
    (
        "U",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(dimVelocity, Zero)
    );

    auto& stf = tstf.ref();

    forAllConstIters(phaseModels_, iter)
    {
        stf += iter()() * iter()->U();
    }

    return tstf;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseInterSystem::surfaceTensionCoeff(const phasePairKey& key) const
{
    return surfaceTensionModels_[key]->sigma();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::coeffs
(
    const word& key
) const
{
    return 1.0/(phaseModels_[key]->thermo().rho());
}


void Foam::multiphaseInterSystem::addInterfacePorosity(fvVectorMatrix& UEqn)
{
    const scalarField& Vc = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    forAllConstIters(phaseModels_, iteri)
    {
        const multiphaseInter::phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != phaseModels_.cend(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const multiphaseInter::phaseModel& phasek = iterk()();

                // Phase i and k
                const phasePairKey keyik
                (
                    phasei.name(),
                    phasek.name(),
                    false
                );

                if (interfacePorousModelTable_.found(keyik))
                {
                    autoPtr<porousModel>& interfacePtr =
                        interfacePorousModelTable_[keyik];

                    Udiag += Vc*interfacePtr->S();
                }
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::K
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    // Simple expression for curvature
    return -fvc::div(tnHatfv.ref() & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseInterSystem::nearInterface
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return
    (
        pos(alpha1 - 0.1)*pos(0.9 - alpha1)
       *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseInterSystem::nearInterface() const
{
    auto tnearInt = volScalarField::New
    (
        "nearInterface",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    auto& nearInt = tnearInt.ref();

    forAllConstIters(phaseModels_, iter1)
    {
        const volScalarField& alpha1 = iter1()();

        auto iter2 = iter1;

        for (++iter2; iter2 != phaseModels_.cend(); ++iter2)
        {
            const volScalarField& alpha2 = iter2()();

            nearInt +=
            (
                pos(alpha1 - 0.1)*pos(0.9 - alpha1)
               *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
            );
        }
    }

    return tnearInt;
}


Foam::tmp<Foam::volVectorField> Foam::multiphaseInterSystem::nVolHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    const volScalarField alpha1m
    (
        clamp(alpha1, zero_one{})
    );

    const volScalarField alpha2m
    (
        clamp(alpha2, zero_one{})
    );

    const volVectorField gradAlphaf
    (
        alpha2m*(fvc::grad(alpha1m))
      - alpha1m*(fvc::grad(alpha2m))
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseInterSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{

    const volScalarField alpha1b
    (
        clamp(alpha1, zero_one{})
    );

    const volScalarField alpha2b
    (
        clamp(alpha2, zero_one{})
    );

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2b)*fvc::interpolate(fvc::grad(alpha1b))
      - fvc::interpolate(alpha1b)*fvc::interpolate(fvc::grad(alpha2b))
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseInterSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


bool Foam::multiphaseInterSystem::read()
{
    if (regIOobject::read())
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
