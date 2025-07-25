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

#include "multiphaseSystem.H"

#include "fixedValueFvsPatchFields.H"
#include "Time.H"
#include "subCycle.H"
#include "fvcMeshPhi.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvcAverage.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "CMULES.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseInter
{
    defineTypeNameAndDebug(multiphaseSystem, 0);
    defineRunTimeSelectionTable(multiphaseSystem, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInter::multiphaseSystem::multiphaseSystem
(
     const fvMesh& mesh
)
:
    multiphaseInterSystem(mesh),
    cAlphas_(),
    ddtAlphaMax_(0.0),
    limitedPhiAlphas_(phaseModels_.size()),
    Su_(phaseModels_.size()),
    Sp_(phaseModels_.size())
{
    label phasei = 0;
    phases_.setSize(phaseModels_.size());
    forAllIters(phaseModels_, iter)
    {
        multiphaseInter::phaseModel& pm = iter()();
        phases_.set(phasei++, &pm);
    }

    mesh.solverDict("alpha").readEntry("cAlphas", cAlphas_);

    // Initiate Su and Sp
    forAllConstIters(phaseModels_, iter)
    {
        const multiphaseInter::phaseModel& pm = iter()();

        Su_.emplace
        (
            pm.name(),
            // volScalarField::Internal
            mesh_.newIOobject("Su" + pm.name()),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        Sp_.emplace
        (
            pm.name(),
            // volScalarField::Internal
            mesh_.newIOobject("Sp" + pm.name()),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseInter::multiphaseSystem::calculateSuSp()
{
    this->alphaTransfer(Su_, Sp_);
}


void Foam::multiphaseInter::multiphaseSystem::solve()
{
    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

    volScalarField& alpha = phases_.first();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            mesh_.newIOobject("rhoPhiSum"),
            mesh_,
            dimensionedScalar(rhoPhi_.dimensions(), Zero)
        );

        dimensionedScalar totalDeltaT = mesh_.time().deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();
            rhoPhiSum += (mesh_.time().deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas();
    }

}

void Foam::multiphaseInter::multiphaseSystem::solveAlphas()
{

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    alphaControls.readEntry("cAlphas", cAlphas_);
    label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

    PtrList<surfaceScalarField> phiAlphaCorrs(phases_.size());

    const surfaceScalarField& phi = this->phi();

    surfaceScalarField phic(mag((phi)/mesh_.magSf()));

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        label phasei = 0;
        for (multiphaseInter::phaseModel& phase1 : phases_)
        {
            const volScalarField& alpha1 = phase1;

            phiAlphaCorrs.set
            (
                phasei,
                new surfaceScalarField
                (
                    "phi" + alpha1.name() + "Corr",
                    fvc::flux
                    (
                        phi,
                        alpha1,
                        "div(phi," + alpha1.name() + ')'
                    )
                )
            );

            surfaceScalarField& phiAlphaCorr = phiAlphaCorrs[phasei];

            for (multiphaseInter::phaseModel& phase2 : phases_)
            {
                const volScalarField& alpha2 = phase2;

                if (&phase2 == &phase1) continue;

                const phasePairKey key12(phase1.name(), phase2.name());

                if (!cAlphas_.found(key12))
                {
                    FatalErrorInFunction
                        << "Phase compression factor (cAlpha) not found for : "
                        << key12
                        << exit(FatalError);
                }
                scalar cAlpha = cAlphas_.find(key12)();

                phic = min(cAlpha*phic, max(phic));

                surfaceScalarField phir(phic*nHatf(alpha1, alpha2));

                word phirScheme
                (
                    "div(phir," + alpha2.name() + ',' + alpha1.name() + ')'
                );

                phiAlphaCorr += fvc::flux
                (
                   -fvc::flux(-phir, alpha2, phirScheme),
                    alpha1,
                    phirScheme
                );
            }

            // Ensure that the flux at inflow BCs is preserved
            forAll(phiAlphaCorr.boundaryField(), patchi)
            {
                fvsPatchScalarField& phiAlphaCorrp =
                    phiAlphaCorr.boundaryFieldRef()[patchi];

                if (!phiAlphaCorrp.coupled())
                {
                    const scalarField& phi1p = phi.boundaryField()[patchi];
                    const scalarField& alpha1p =
                        alpha1.boundaryField()[patchi];

                    forAll(phiAlphaCorrp, facei)
                    {
                        if (phi1p[facei] < 0)
                        {
                            phiAlphaCorrp[facei] = alpha1p[facei]*phi1p[facei];
                        }
                    }
                }
            }

            ++phasei;
        }

        // Set Su and Sp to zero
        for (const multiphaseInter::phaseModel& phase : phases_)
        {
            Su_[phase.name()] = dimensionedScalar("Su", dimless/dimTime, Zero);
            Sp_[phase.name()] = dimensionedScalar("Sp", dimless/dimTime, Zero);

            // Add alpha*div(U)
            //const volScalarField& alpha = phase;
            //Su_[phase.name()] += fvc::div(phi)*clamp(alpha, zero_one{});
        }

        // Fill Su and Sp
        calculateSuSp();

        // Limit phiAlphaCorr on each phase
        phasei = 0;
        for (multiphaseInter::phaseModel& phase : phases_)
        {
            volScalarField& alpha1 = phase;

            surfaceScalarField& phiAlphaCorr = phiAlphaCorrs[phasei];

            volScalarField::Internal& Su = Su_[phase.name()];
            volScalarField::Internal& Sp = Sp_[phase.name()];

            MULES::limit
            (
                1.0/mesh_.time().deltaTValue(),
                geometricOneField(),
                alpha1,
                phi,
                phiAlphaCorr,
                Sp,
                Su,
                oneField(),
                zeroField(),
                true
            );
            ++phasei;
        }

        MULES::limitSum(phiAlphaCorrs);

        volScalarField sumAlpha
        (
            mesh_.newIOobject("sumAlpha"),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );

        phasei = 0;
        for (multiphaseInter::phaseModel& phase : phases_)
        {
            volScalarField& alpha1 = phase;

            const volScalarField::Internal& Su = Su_[phase.name()];

            const volScalarField::Internal& Sp = Sp_[phase.name()];

            surfaceScalarField& phiAlpha = phiAlphaCorrs[phasei];

            // Add a bounded upwind U-mean flux
            //phiAlpha += upwind<scalar>(mesh_, phi).flux(alpha1);
            fvScalarMatrix alpha1Eqn
            (
                fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(alpha1)
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh_,
                    phi,
                    upwind<scalar>(mesh_, phi)
                ).fvmDiv(phi, alpha1)
              ==
                 Su + fvm::Sp(Sp, alpha1)
            );

            alpha1Eqn.boundaryManipulate(alpha1.boundaryFieldRef());

            alpha1Eqn.solve();

            phiAlpha += alpha1Eqn.flux();

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                phiAlpha,
                Sp,
                Su,
                oneField(),
                zeroField()
            );

            phase.alphaPhi() = phiAlpha;

            ++phasei;
        }

        if (acorr == nAlphaCorr - 1)
        {
            volScalarField sumAlpha
            (
                mesh_.newIOobject("sumAlpha"),
                mesh_,
                dimensionedScalar(dimless, Zero)
            );

            // Reset rhoPhi
            rhoPhi_ = dimensionedScalar("rhoPhi", dimMass/dimTime, Zero);

            for (multiphaseInter::phaseModel& phase : phases_)
            {
                volScalarField& alpha1 = phase;
                sumAlpha += alpha1;

                // Update rhoPhi
                rhoPhi_ += fvc::interpolate(phase.rho()) * phase.alphaPhi();
            }

            Info<< "Phase-sum volume fraction, min, max = "
                << sumAlpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(sumAlpha).value()
                << ' ' << max(sumAlpha).value()
                << endl;

            volScalarField sumCorr(1.0 - sumAlpha);

            for (multiphaseInter::phaseModel& phase : phases_)
            {
                volScalarField& alpha = phase;
                alpha += alpha*sumCorr;

                Info<< alpha.name() << " volume fraction = "
                    << alpha.weightedAverage(mesh_.V()).value()
                    << "  Min(alpha) = " << min(alpha).value()
                    << "  Max(alpha) = " << max(alpha).value()
                    << endl;
            }
        }
    }
}


const Foam::UPtrList<Foam::multiphaseInter::phaseModel>&
Foam::multiphaseInter::multiphaseSystem::phases() const
{
    return phases_;
}


Foam::UPtrList<Foam::multiphaseInter::phaseModel>&
Foam::multiphaseInter::multiphaseSystem::phases()
{
    return phases_;
}


const Foam::multiphaseInter::phaseModel&
Foam::multiphaseInter::multiphaseSystem::phase(const label i) const
{
    return phases_[i];
}


Foam::multiphaseInter::phaseModel&
Foam::multiphaseInter::multiphaseSystem::phase(const label i)
{
    return phases_[i];
}


Foam::dimensionedScalar
Foam::multiphaseInter::multiphaseSystem::ddtAlphaMax() const
{
    return ddtAlphaMax_;
}


Foam::scalar
Foam::multiphaseInter::multiphaseSystem::maxDiffNo() const
{
    auto iter = phaseModels_.cbegin();

    scalar maxVal = max(iter()->diffNo()).value();

    for (++iter; iter != phaseModels_.cend(); ++iter)
    {
        maxVal = max(maxVal, max(iter()->diffNo()).value());
    }

    return maxVal * mesh_.time().deltaTValue();
}


const Foam::multiphaseInter::multiphaseSystem::compressionFluxTable&
Foam::multiphaseInter::multiphaseSystem::limitedPhiAlphas() const
{
    return limitedPhiAlphas_;
}


Foam::multiphaseInter::multiphaseSystem::SuSpTable&
Foam::multiphaseInter::multiphaseSystem::Su()
{
    return Su_;
}


Foam::multiphaseInter::multiphaseSystem::SuSpTable&
Foam::multiphaseInter::multiphaseSystem::Sp()
{
    return Sp_;
}


bool Foam::multiphaseInter::multiphaseSystem::read()
{
    return true;
}


// ************************************************************************* //
