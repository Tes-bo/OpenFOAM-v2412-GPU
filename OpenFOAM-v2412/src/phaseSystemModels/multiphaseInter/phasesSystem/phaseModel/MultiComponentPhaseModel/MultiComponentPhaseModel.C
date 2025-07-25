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

#include "MultiComponentPhaseModel.H"

#include "multiphaseInterSystem.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcDDt.H"
#include "fvMatrix.H"
#include "fvcFlux.H"
#include "CMULES.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
MultiComponentPhaseModel
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    species_(),
    inertIndex_(-1),
    addDiffusion_(false)
{
    thermoPtr_.reset
    (
        phaseThermo::New
        (
            fluid.mesh(),
            phaseName,
            basicThermo::phasePropertyName(basicThermo::dictName, phaseName)
        )
    );

    if (thermoPtr_->composition().species().empty())
    {
        FatalErrorInFunction
            << " The selected thermo is pure. Use a multicomponent thermo."
            << exit(FatalError);
    }

    species_ = thermoPtr_->composition().species();

    inertIndex_ = species_.find(thermoPtr_().template get<word>("inertSpecie"));

    addDiffusion_ =
        thermoPtr_().template getOrDefault<bool>("addDiffusion", false);

    Sct_ = thermoPtr_().template getOrDefault<scalar>("Sct", 1.0);

    X_.setSize(thermoPtr_->composition().species().size());

    // Initiate X's using Y's to set BC's
    forAll(species_, i)
    {
        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("X" + species_[i], phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Y()[i]
            )
        );
    }

    // Init vol fractions from mass fractions
    calculateVolumeFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>
::calculateVolumeFractions()
{
    volScalarField Xtotal(0.0*X_[0]);
    const volScalarField W(thermo().W());

    forAll(X_, i)
    {
        const dimensionedScalar Wi
        (
            "Wi",
            dimMass/dimMoles,
            thermo().composition().W(i)
        );

        if (i != inertIndex_)
        {
            X_[i] = W*Y()[i]/Wi;
            Xtotal += X_[i];
            X_[i].correctBoundaryConditions();

            const volScalarField::Boundary& YBf = Y()[i].boundaryField();

            forAll(YBf, patchi)
            {
                const fvPatchScalarField& YPf = YBf[patchi];
                if (YPf.fixesValue())
                {
                    scalarField& xbf = X_[i].boundaryFieldRef()[patchi];
                    const scalarField& ybf = Y()[i].boundaryField()[patchi];
                    const scalarField& Wbf = W.boundaryField()[patchi];
                    forAll(xbf, facei)
                    {
                        xbf[facei] = Wbf[facei]*ybf[facei]/Wi.value();
                    }
                }
            }
        }
    }
    X_[inertIndex_] = 1.0 - Xtotal;
    X_[inertIndex_].correctBoundaryConditions();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
calculateMassFractions()
{
    volScalarField W(X_[0]*thermo().composition().W(0));
    for(label i=1; i< species_.size(); i++)
    {
        W += X_[i]*thermo().composition().W(i);
    }

    forAll(Y(), i)
    {
        Y()[i] = X_[i]*thermo().composition().W(i)/W;

        Info<< Y()[i].name() << " mass fraction = "
            << "  Min(Y) = " << min(Y()[i]).value()
            << "  Max(Y) = " << max(Y()[i]).value()
            << endl;
    }
}


template<class BasePhaseModel, class phaseThermo>
const phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::correct()
{
    return thermo().correct();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::solveYi
(
    PtrList<volScalarField::Internal>& Su,
    PtrList<volScalarField::Internal>& Sp
)
{
    const volScalarField& alpha1 = *this;

    const fvMesh& mesh = alpha1.mesh();

    const dictionary& MULEScontrols = mesh.solverDict(alpha1.name());

    scalar cAlpha(MULEScontrols.get<scalar>("cYi"));

    PtrList<surfaceScalarField> phiYiCorrs(species_.size());
    const surfaceScalarField& phi = this->fluid().phi();

    surfaceScalarField phic(mag((phi)/mesh.magSf()));

    phic = min(cAlpha*phic, max(phic));

    surfaceScalarField phir(0.0*phi);

    forAllConstIter
    (
        multiphaseInterSystem::phaseModelTable,this->fluid().phases(),iter2
    )
    {
        const volScalarField& alpha2 = iter2()();
        if (&alpha2 == &alpha1)
        {
            continue;
        }

        phir += phic*this->fluid().nHatf(alpha1, alpha2);
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phirBf = phir.boundaryFieldRef();
    forAll(phir.boundaryField(), patchi)
    {
        fvsPatchScalarField& phirp = phirBf[patchi];

        if (!phirp.coupled())
        {
            phirp == 0;
        }
    }

    word phirScheme("div(Yiphir," + alpha1.name() + ')');

    forAll(X_, i)
    {
        if (inertIndex_ != i)
        {
            volScalarField& Yi = X_[i];

            phiYiCorrs.set
            (
                i,
                new surfaceScalarField
                (
                    "phi" + Yi.name() + "Corr",
                    fvc::flux
                    (
                        phi,
                        Yi,
                        "div(phi," + Yi.name() + ')'
                    )
                )
            );

            surfaceScalarField& phiYiCorr = phiYiCorrs[i];

            forAllConstIter
            (
                multiphaseInterSystem::phaseModelTable,
                this->fluid().phases(),
                iter2
            )
            {
                //const volScalarField& alpha2 = iter2()().oldTime();
                const volScalarField& alpha2 = iter2()();

                if (&alpha2 == &alpha1)
                {
                    continue;
                }

                phiYiCorr += fvc::flux
                (
                   -fvc::flux(-phir, alpha2, phirScheme),
                    Yi,
                    phirScheme
                );
            }

            // Ensure that the flux at inflow BCs is preserved
            forAll(phiYiCorr.boundaryField(), patchi)
            {
                fvsPatchScalarField& phiYiCorrp =
                    phiYiCorr.boundaryFieldRef()[patchi];

                if (!phiYiCorrp.coupled())
                {
                    const scalarField& phi1p = phi.boundaryField()[patchi];
                    const scalarField& Yip = Yi.boundaryField()[patchi];

                    forAll(phiYiCorrp, facei)
                    {
                        if (phi1p[facei] < 0)
                        {
                            phiYiCorrp[facei] = Yip[facei]*phi1p[facei];
                        }
                    }
                }
            }

            MULES::limit
            (
                1.0/mesh.time().deltaTValue(),
                geometricOneField(),
                Yi,
                phi,
                phiYiCorr,
                Sp[i],
                Su[i],
                oneField(),
                zeroField(),
                true
            );
        }
    }

    volScalarField Yt(0.0*X_[0]);

    scalar nYiSubCycles
    (
        MULEScontrols.getOrDefault<scalar>("nYiSubCycles", 1)
    );

    forAll(X_, i)
    {
        if (inertIndex_ != i)
        {
            volScalarField& Yi = X_[i];

            fvScalarMatrix YiEqn
            (
                fv::EulerDdtScheme<scalar>(mesh).fvmDdt(Yi)
              + fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phi,
                    upwind<scalar>(mesh, phi)
                ).fvmDiv(phi, Yi)
                 ==
                Su[i] + fvm::Sp(Sp[i], Yi)
            );

            YiEqn.solve();

            surfaceScalarField& phiYiCorr = phiYiCorrs[i];

            // Add a bounded upwind U-mean flux
            phiYiCorr += YiEqn.flux();

            if (nYiSubCycles > 1)
            {
                for
                (
                    subCycle<volScalarField> YiSubCycle(Yi, nYiSubCycles);
                    !(++YiSubCycle).end();
                )
                {
                    MULES::explicitSolve
                    (
                        geometricOneField(),
                        Yi,
                        phi,
                        phiYiCorr,
                        Sp[i],
                        Su[i],
                        oneField(),
                        zeroField()
                    );
                }
            }
            else
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    Yi,
                    phi,
                    phiYiCorr,
                    Sp[i],
                    Su[i],
                    oneField(),
                    zeroField()
                );
            }

            if (addDiffusion_)
            {
                const volScalarField& alpha = *this;
                fvScalarMatrix YiDiffEqn
                (
                    fvm::ddt(Yi) - fvc::ddt(Yi)
                  - fvm::laplacian
                    (
                        alpha*this->fluid().turbulence()->nut()/Sct_,
                        Yi
                    )
                );

                YiDiffEqn.solve("diffusion" + Yi.name());
            }

            Yt += Yi;
        }
    }

    X_[inertIndex_] = scalar(1) - Yt;
    X_[inertIndex_].clamp_min(0);

    calculateMassFractions();
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::label
Foam::MultiComponentPhaseModel<BasePhaseModel, phaseThermo>::inertIndex() const
{
    return inertIndex_;
}


// ************************************************************************* //
