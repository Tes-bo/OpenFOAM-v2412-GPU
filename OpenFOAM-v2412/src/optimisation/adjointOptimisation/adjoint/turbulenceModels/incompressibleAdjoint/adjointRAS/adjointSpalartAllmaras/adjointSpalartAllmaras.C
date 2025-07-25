/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "adjointSpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "nutUSpaldingWallFunctionFvPatchScalarField.H"
#include "boundaryAdjointContribution.H"
#include "coupledFvPatch.H"
#include "ATCModel.H"
#include "fvOptions.H"
#include "sensitivityTopO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{
namespace adjointRASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointSpalartAllmaras, 0);
addToRunTimeSelectionTable
(
    adjointRASModel,
    adjointSpalartAllmaras,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * Primal Spalart - Allmaras * * * * * * * * * * * * //

tmp<volScalarField> adjointSpalartAllmaras::chi() const
{
    return nuTilda()/nu();
}


tmp<volScalarField> adjointSpalartAllmaras::fv1(const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> adjointSpalartAllmaras::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


tmp<volScalarField> adjointSpalartAllmaras::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField Omega(::sqrt(2.0)*mag(skew(gradU_)));

    return
    (
        max
        (
            Omega
          + fv2(chi, fv1)*nuTilda()/sqr(kappa_*y_),
            Cs_*Omega
        )
    );
}


tmp<volScalarField> adjointSpalartAllmaras::r
(
    const volScalarField& Stilda
) const
{
    auto tr =
        tmp<volScalarField>::New
        (
            min
            (
                nuTilda()/(max(Stilda, minStilda_)*sqr(kappa_*y_)),
                scalar(10)
            )
        );
    tr.ref().boundaryFieldRef() == Zero;

    return tr;
}


tmp<volScalarField> adjointSpalartAllmaras::fw
(
    const volScalarField& Stilda
) const
{
    const volScalarField g(r_ + Cw2_*(pow6(r_) - r_));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


tmp<volScalarField> adjointSpalartAllmaras::DnuTildaEff() const
{
    return
        tmp<volScalarField>::New
            ("DnuTildaEff", (nuTilda() + this->nu())/sigmaNut_);
}


const volScalarField& adjointSpalartAllmaras::nuTilda() const
{
    return primalVars_.RASModelVariables()().TMVar1();
}


const volScalarField& adjointSpalartAllmaras::nutRef() const
{
    return primalVars_.RASModelVariables()().nutRef();
}


// * * * * * * * * * * *  Adjoint Spalart - Allmaras * * * * * * * * * * * * //

tmp<volScalarField> adjointSpalartAllmaras::dFv1_dChi
(
    const volScalarField& chi
) const
{
    volScalarField chi3(pow3(chi));

    return 3.0*pow3(Cv1_)*sqr(chi/(chi3 + pow3(Cv1_)));
}


tmp<volScalarField> adjointSpalartAllmaras::dFv2_dChi
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volScalarField& dFv1dChi
) const
{
    return (chi*chi*dFv1dChi - 1.)/sqr(1. + chi*fv1);
}


tmp<volScalarField> adjointSpalartAllmaras::dStilda_dOmega
(
    const volScalarField& Omega,
    const volScalarField& fv2
) const
{
    volScalarField fieldSwitch
    (
        Omega + fv2*nuTilda()/sqr(kappa_*y_) - Cs_*Omega
    );

    return pos(fieldSwitch) + neg(fieldSwitch)*Cs_;
}


tmp<volScalarField> adjointSpalartAllmaras::dStilda_dNuTilda
(
    const volScalarField& Omega,
    const volScalarField& fv2,
    const volScalarField& dFv2dChi
) const
{
    volScalarField invDenom(1./sqr(kappa_*y_));
    volScalarField fieldSwitch(Omega + fv2*nuTilda()*invDenom - Cs_*Omega);

    return pos(fieldSwitch)*(dFv2dChi*nuTilda()*invDenom/nu() + fv2*invDenom);
}


tmp<volScalarField> adjointSpalartAllmaras::dStilda_dDelta
(
    const volScalarField& Omega,
    const volScalarField& fv2
) const
{
    volScalarField aux(fv2*nuTilda()/sqr(kappa_*y_));
    volScalarField fieldSwitch(Omega + aux - Cs_*Omega);

    return - 2.*pos(fieldSwitch)*aux/y_;
}


tmp<volScalarField> adjointSpalartAllmaras::dr_dNuTilda
(
    const volScalarField& Stilda
) const
{
    tmp<volScalarField> tdrdNutilda
    (
        1./(max(Stilda, minStilda_)*sqr(kappa_*y_))
        *(scalar(10) - r_)/(scalar(10) - r_ + SMALL)
    );
    tdrdNutilda.ref().boundaryFieldRef() == Zero;

    return tdrdNutilda;
}


tmp<volScalarField> adjointSpalartAllmaras::dr_dStilda
(
    const volScalarField& Stilda
) const
{
    tmp<volScalarField> tdrdStilda
    (
      - nuTilda()/sqr(max(Stilda, minStilda_)*kappa_*y_)
       *(scalar(10) - r_)/(scalar(10) - r_ + SMALL)
    );
    tdrdStilda.ref().boundaryFieldRef() == Zero;

    return tdrdStilda;
}


tmp<volScalarField> adjointSpalartAllmaras::dr_dDelta
(
    const volScalarField& Stilda
) const
{
    tmp<volScalarField> tdrdDelta
    (
      - 2.*nuTilda()/(max(Stilda, minStilda_)*sqr(kappa_*y_)*y_)
       *(scalar(10) - r_)/(scalar(10) - r_ + SMALL)
    );
    tdrdDelta.ref().boundaryFieldRef() == Zero;

    return tdrdDelta;
}


tmp<volScalarField> adjointSpalartAllmaras::dfw_dr
(
    const volScalarField& Stilda
) const
{
    volScalarField g(r_ + Cw2_*(pow6(r_) - r_));

    dimensionedScalar pow6Cw3 = pow6(Cw3_);
    volScalarField pow6g(pow6(g));

    return
        pow6Cw3/(pow6g + pow6Cw3)
       *pow((1.0 + pow6Cw3)/(pow6g + pow6Cw3), 1.0/6.0)
       *(1.0 + Cw2_*(6.0*pow5(r_) - 1.0));
}


tmp<volScalarField> adjointSpalartAllmaras::dfw_dNuTilda
(
    const volScalarField& Stilda,
    const volScalarField& dfwdr,
    const volScalarField& dStildadNuTilda
) const
{
    volScalarField invDenom(1./sqr(kappa_*y_));

    return
        dfwdr*(dr_dNuTilda(Stilda) + dr_dStilda(Stilda)*dStildadNuTilda);
}


tmp<volScalarField> adjointSpalartAllmaras::dfw_dOmega
(
    const volScalarField& Stilda,
    const volScalarField& dfwdr,
    const volScalarField& dStildadOmega
) const
{
    return dfwdr*dr_dStilda(Stilda)*dStildadOmega;
}


tmp<volScalarField> adjointSpalartAllmaras::dfw_dDelta
(
    const volScalarField& Stilda,
    const volScalarField& dfwdr,
    const volScalarField& dStildadDelta
) const
{
    return dfwdr*(dr_dDelta(Stilda) + dr_dStilda(Stilda)*dStildadDelta);
}


tmp<volScalarField> adjointSpalartAllmaras::dD_dNuTilda
(
    const volScalarField& fw,
    const volScalarField& dfwdNuTilda
) const
{
    return Cw1_*(nuTilda()*dfwdNuTilda + fw)/sqr(y_);
}


tmp<volScalarField> adjointSpalartAllmaras::dP_dNuTilda
(
    const volScalarField& dStildadNuTilda
) const
{
    return - Cb1_*dStildadNuTilda;
}


tmp<volScalarField> adjointSpalartAllmaras::dnut_dNuTilda
(
    const volScalarField& fv1,
    const volScalarField& dFv1dChi
) const
{
    return dFv1dChi*nuTilda()/nu() + fv1;
}


tmp<volVectorField> adjointSpalartAllmaras::conservativeMomentumSource()
{
    // Store boundary field of the conservative part,
    // for use in adjoint outlet boundary conditions
    forAll(mesh_.boundary(), pI)
    {
        const fvPatch& patch = mesh_.boundary()[pI];
        if(!isA<coupledFvPatch>(patch))
        {
            tmp<vectorField> tnf = patch.nf();
            adjMomentumBCSourcePtr_()[pI] =
                (tnf & momentumSourceMult_.boundaryField()[pI])
               *nuaTilda().boundaryField()[pI];
        }
    }

    return fvc::div(momentumSourceMult_*nuaTilda());
}


void adjointSpalartAllmaras::updatePrimalRelatedFields()
{
    if (changedPrimalSolution_)
    {
        Info<< "Updating primal-based fields of the adjoint turbulence "
            << "model ..." << endl;

        // Grab references
        const volVectorField& U = primalVars_.U();

        // Update gradient fields
        gradU_ = mask_*fvc::grad(U, "gradUStilda");
        gradNuTilda_ = fvc::grad(nuTilda());

        const volScalarField Omega(::sqrt(2.0)*mag(skew(gradU_)));

        // Primal SA fields
        volScalarField chi(this->chi());
        volScalarField fv1(this->fv1(chi));
        volScalarField fv2(this->fv2(chi, fv1));
        Stilda_ = Stilda(chi, fv1);
        r_ = r(Stilda_);
        fw_ = this->fw(Stilda_);

        // Derivatives of primal fields wrt to nuTilda
        volScalarField dFv1_dChi(this->dFv1_dChi(chi));
        volScalarField dFv2_dChi(this->dFv2_dChi(chi, fv1, dFv1_dChi));
        volScalarField dStilda_dNuTilda
            (this->dStilda_dNuTilda(Omega, fv2, dFv2_dChi));
        volScalarField dfw_dr(this->dfw_dr(Stilda_));
        volScalarField dfw_dNuTilda
            (this->dfw_dNuTilda(Stilda_, dfw_dr, dStilda_dNuTilda));

        // Fields to be used in the nuaTilda equation
        symmAdjointProductionU_ =
            symm(mask_*fvc::grad(U, "adjointProductionU"));

        productionDestructionSource_ =
            nuTilda()
           *(
                dD_dNuTilda(fw_, dfw_dNuTilda)
              + dP_dNuTilda(dStilda_dNuTilda)
            );

        Cdnut_ = dnut_dNuTilda(fv1, dFv1_dChi);

        // Constant multiplier in the adjoint momentum source term
        volScalarField dStilda_dOmega(this->dStilda_dOmega(Omega, fv2));
        volScalarField dfw_dOmega
            (this->dfw_dOmega(Stilda_, dfw_dr, dStilda_dOmega));

        momentumSourceMult_ =
            2.*skew(gradU_)
           /(Omega + dimensionedScalar("SMALL", Omega.dimensions(), SMALL))
           *(
              - Cb1_*nuTilda()*dStilda_dOmega
              + Cw1_*sqr(nuTilda()/y_)*dfw_dOmega
            );

        // Set changedPrimalSolution_ to false to avoid recomputing these
        // fields unless the primal has changed
        changedPrimalSolution_ = false;
    }
}


tmp<volScalarField> adjointSpalartAllmaras::allocateMask()
{
    if (limitAdjointProduction_)
    {
        return ATCModel::createLimiter(mesh_, coeffDict_);
    }

    return volScalarField::New
    (
        "unitMask",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar("unit", dimless, scalar(1))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointSpalartAllmaras::adjointSpalartAllmaras
(
    incompressibleVars& primalVars,
    incompressibleAdjointMeanFlowVars& adjointVars,
    objectiveManager& objManager,
    const word& adjointTurbulenceModelName,
    const word& modelName
)
:
    adjointRASModel
    (
        modelName,
        primalVars,
        adjointVars,
        objManager,
        adjointTurbulenceModelName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),

    limitAdjointProduction_
    (
        coeffDict_.getOrDefault("limitAdjointProduction", true)
    ),

    y_(primalVars_.RASModelVariables()().d()),

    mask_(allocateMask()),

    symmAdjointProductionU_
    (
        IOobject
        (
            "symmAdjointProductionU",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor(dimless/dimTime, Zero)
    ),

    productionDestructionSource_
    (
        IOobject
        (
            "productionDestructionSource",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),

    Stilda_
    (
        IOobject
        (
            "Stilda",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),

    r_
    (
        IOobject
        (
            "r",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    fw_
    (
        IOobject
        (
            "fw",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    Cdnut_
    (
        IOobject
        (
            "Cdnut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    momentumSourceMult_
    (
        IOobject
        (
            "momentumSourceMult",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimLength)/dimTime, Zero)
    ),

    gradU_(fvc::grad(primalVars.U())),
    gradNuTilda_(fvc::grad(nuTilda())),
    minStilda_("SMALL", Stilda_.dimensions(), SMALL)
{
    adjointTMVariablesBaseNames_.setSize(1);
    adjointTMVariablesBaseNames_[0] = "nuaTilda";
    // Read nuaTilda field and reset pointer to the first
    // adjoint turbulence model variable
    variablesSet::setField
    (
        adjointTMVariable1Ptr_,
        mesh_,
        "nuaTilda",
        adjointVars.solverName(),
        adjointVars.useSolverNameForFields()
    );

    setMeanFields();

    // Set the includeDistance to true, to allow for the automatic solution
    // of the adjoint eikonal equation when computing sensitivities
    includeDistance_ = true;

    // Update the primal related fields here so that functions computing
    // sensitivities have the updated fields in case of continuation
    updatePrimalRelatedFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> adjointSpalartAllmaras::devReff() const
{
    const volVectorField& Ua = adjointVars_.UaInst();
    return devReff(Ua);
}


tmp<volSymmTensorField> adjointSpalartAllmaras::devReff
(
    const volVectorField& U
) const
{
    return volSymmTensorField::New
    (
        "devRhoReff",
        IOobject::NO_REGISTER,
        -nuEff()*devTwoSymm(fvc::grad(U))
    );
}


tmp<fvVectorMatrix> adjointSpalartAllmaras::divDevReff(volVectorField& Ua) const
{
    tmp<volScalarField> tnuEff(nuEff());
    const volScalarField& nuEff = tnuEff();

    return
    (
      - fvm::laplacian(nuEff, Ua)
      - fvc::div(nuEff*dev(fvc::grad(Ua)().T()))
    );
}


tmp<volVectorField> adjointSpalartAllmaras::adjointMeanFlowSource()
{
    addProfiling
    (
        adjointSpalartAllmaras,
        "adjointSpalartAllmaras::addMomentumSource"
    );
    // cm formulation
    //return (- nuTilda()*fvc::grad(nuaTilda() - conservativeMomentumSource());

    // ncm formulation
    return (nuaTilda()*gradNuTilda_ - conservativeMomentumSource());
}


tmp<volScalarField> adjointSpalartAllmaras::nutJacobianTMVar1() const
{
    volScalarField chi(this->chi());
    volScalarField fv1(this->fv1(chi));
    volScalarField dFv1_dChi(this->dFv1_dChi(chi));

    return dnut_dNuTilda(fv1, dFv1_dChi);
}


tmp<scalarField> adjointSpalartAllmaras::diffusionCoeffVar1(label patchI) const
{
    auto tdiffCoeff =
        tmp<scalarField>::New(mesh_.boundary()[patchI].size(), Zero);

    scalarField& diffCoeff = tdiffCoeff.ref();

    diffCoeff =
        (nuTilda().boundaryField()[patchI] + nu()().boundaryField()[patchI])
       /sigmaNut_.value();

    return tdiffCoeff;
}


const boundaryVectorField&
adjointSpalartAllmaras::adjointMomentumBCSource() const
{
    // Computed in conservativeMomentumSource
    return adjMomentumBCSourcePtr_();
}


const boundaryVectorField& adjointSpalartAllmaras::wallShapeSensitivities()
{
    boundaryVectorField& wallShapeSens = wallShapeSensitivitiesPtr_();

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        tmp<vectorField> tnf = patch.nf();
        if (isA<wallFvPatch>(patch) && patch.size() != 0)
        {
            wallShapeSens[patchI] =
              - nuaTilda().boundaryField()[patchI].snGrad()
               *diffusionCoeffVar1(patchI)
               *nuTilda().boundaryField()[patchI].snGrad()*tnf;
        }
    }

    return wallShapeSens;
}


const boundaryVectorField& adjointSpalartAllmaras::wallFloCoSensitivities()
{
    boundaryVectorField& wallFloCoSens = wallFloCoSensitivitiesPtr_();

    forAll(mesh_.boundary(), patchI)
    {
        tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();

        wallFloCoSens[patchI] =
            nuaTilda().boundaryField()[patchI]
           *nuTilda().boundaryField()[patchI]*tnf;
    }

    return wallFloCoSens;
}


tmp<volScalarField> adjointSpalartAllmaras::distanceSensitivities()
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.Ua();

    // Primal SA fields
    volScalarField chi(this->chi());
    volScalarField fv1(this->fv1(chi));
    volScalarField fv2(this->fv2(chi, fv1));
    volScalarField Omega(::sqrt(2.0)*mag(gradU_));

    // Derivatives of primal fields wrt to nuTilda
    volScalarField dFv1_dChi(this->dFv1_dChi(chi));
    volScalarField dFv2_dChi(this->dFv2_dChi(chi, fv1, dFv1_dChi));
    volScalarField dStilda_dDelta(this->dStilda_dDelta(Omega, fv2));
    volScalarField dfw_dr(this->dfw_dr(Stilda_));
    volScalarField dfw_dDelta
        (this->dfw_dDelta(Stilda_, dfw_dr, dStilda_dDelta));

    auto tadjointEikonalSource =
        tmp<volScalarField>::New
        (
            "adjointEikonalSource" + type(),
            (
              - Cb1_*nuTilda()*dStilda_dDelta
              + Cw1_*sqr(nuTilda()/y_)*(dfw_dDelta - 2.*fw_/y_)
            )*nuaTilda()
        );
    volScalarField& adjointEikonalSource = tadjointEikonalSource.ref();

    // if wall functions are used, add appropriate source terms
    typedef nutUSpaldingWallFunctionFvPatchScalarField
        SAwallFunctionPatchField;

    const volScalarField::Boundary& nutBoundary = nutRef().boundaryField();
    const scalarField& V = mesh_.V().field();

    tmp<volScalarField> tnuEff = nuEff();
    const volScalarField& nuEff = tnuEff();

    forAll(nutBoundary, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if
        (
            isA<SAwallFunctionPatchField>(nutBoundary[patchi])
         && patch.size() != 0
        )
        {
            const scalar kappa_(0.41);
            const scalar E_(9.8);
            tmp<vectorField> tnf = patch.nf();
            const vectorField& nf = tnf.ref();
            const scalarField& magSf = patch.magSf();

            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvPatchVectorField& Uap = Ua.boundaryField()[patchi];
            const vectorField Uc(Up.patchInternalField());
            const vectorField Uc_t(Uc - (Uc & nf)*nf);

            // By convention, tf has the direction of the tangent
            // PRIMAL velocity at the first cell off the wall
            const vectorField tf(Uc_t/mag(Uc_t));

            const scalarField nuw(nuEff.boundaryField()[patchi]);
            const scalarField nu(this->nu()().boundaryField()[patchi]);
            const fvPatchScalarField& yC = y()[patchi];

            const scalarField magGradU(mag(Up.snGrad()));

            // Note: What happens in separation?? sign change needed
            const scalarField vtau(sqrt(nuw*magGradU));

            // Note: mag for positive uPlus
            const scalarField uPlus(mag(Uc)/vtau);

            const scalarField yPlus(yC*vtau/nu);
            const scalarField kUu(min(kappa_*uPlus, scalar(50)));
            const scalarField auxA
                ((kappa_/E_)*(exp(kUu) - 1 - kUu - 0.5*kUu*kUu));
            const scalarField Cwf_d(sqr(vtau)/nu/(yPlus+uPlus*(1 + auxA)));

            // Tangential components are according to tf
            autoPtr<boundaryAdjointContribution> boundaryContrPtr
            (
                boundaryAdjointContribution::New
                (
                    "objectiveManager" + objectiveManager_.adjointSolverName(),
                    objectiveManager_.adjointSolverName(),
                    "incompressible",
                    patch
                )
            );
            tmp<vectorField> tsource(boundaryContrPtr->normalVelocitySource());

            const scalarField rt(tsource() & tf);
            const scalarField Uap_t(Uap & tf);

            const labelList& faceCells = patch.faceCells();
            forAll(faceCells, faceI)
            {
                label cellI = faceCells[faceI];
                adjointEikonalSource[cellI] -=
                    2.*( rt[faceI] + Uap_t[faceI] )
                   *vtau[faceI]*Cwf_d[faceI]*magSf[faceI]
                   /V[cellI]; // Divide with cell volume since the term
                              // will be used as a source term in the
                              // adjoint eikonal equation
            }
        }
    }

    return tadjointEikonalSource;
}


tmp<volTensorField> adjointSpalartAllmaras::FISensitivityTerm()
{
    const volVectorField& U = primalVars_.U();

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();
    tmp<volVectorField> tgradNuTilda = fvc::grad(nuTilda());
    volVectorField& gradNuTilda = tgradNuTilda.ref();
    tmp<volVectorField> tgradNuaTilda = fvc::grad(nuaTilda());
    volVectorField::Boundary& gradNuaTildab =
        tgradNuaTilda.ref().boundaryFieldRef();

    // Explicitly correct the boundary gradient to get rid of
    // the tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = patch.nf();
            const vectorField& nf = tnf();
            // gradU:: can cause problems in zeroGradient patches for U
            // and zero fixedValue for nuTilda.
            // S becomes 0 and is used as a denominator in G
            //gradU.boundaryField()[patchI] =
            //  nf * U_.boundaryField()[patchI].snGrad();
            gradNuTilda.boundaryFieldRef()[patchI]  =
                nf*nuTilda().boundaryField()[patchI].snGrad();
            gradNuaTildab[patchI] =
                nf*nuaTilda().boundaryField()[patchI].snGrad();
        }
    }

    // delta vorticity
    volScalarField Omega(::sqrt(2.0)*mag(skew(gradU)));
    volTensorField deltaOmega
    (
        (
           (gradU & gradU)().T() //jk
         - (gradU & gradU.T())   //symmetric
        )
       /(Omega + dimensionedScalar("SMALL", Omega.dimensions(), SMALL))
    );
    tgradU.clear();

    volScalarField chi(this->chi());
    volScalarField fv1(this->fv1(chi));
    volScalarField fv2(this->fv2(chi, fv1));

    volScalarField dfw_dr(this->dfw_dr(Stilda_));
    volScalarField dStilda_dOmega(this->dStilda_dOmega(Omega, fv2));
    volScalarField dfw_dOmega
        (this->dfw_dOmega(Stilda_, dfw_dr, dStilda_dOmega));

    // Assemply of the return field
    auto tFISens
    (
        tmp<volTensorField>::New
        (
            IOobject
            (
                type() + "flowTerm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    volTensorField& FISens = tFISens.ref();
    FISens =
            // jk, cm formulation for the TM model convection
          - (nuaTilda()*(U*gradNuTilda))
            // jk, symmetric in theory
          + nuaTilda()*fvc::grad(DnuTildaEff()*gradNuTilda)().T()
            // jk
          - DnuTildaEff()*(tgradNuaTilda*gradNuTilda)
            // symmetric
          + 2.*nuaTilda()*Cb2_/sigmaNut_*(gradNuTilda*gradNuTilda)
          + (
              - Cb1_*nuTilda()*dStilda_dOmega
              + Cw1_*sqr(nuTilda()/y_)*dfw_dOmega
            )
           *nuaTilda()*deltaOmega; // jk
    FISens.correctBoundaryConditions();

    return tFISens;
}


tmp<scalarField> adjointSpalartAllmaras::topologySensitivities
(
    const word& designVarsName
) const
{
    auto tres(tmp<scalarField>::New(nuTilda().primitiveField().size(), Zero));
    scalarField nuTildaSens
        (nuTilda().primitiveField()*nuaTilda().primitiveField());

    fv::options& fvOptions(fv::options::New(mesh_));
    sensitivityTopO::postProcessSens
    (
        tres.ref(), nuTildaSens, fvOptions, nuTilda().name(), designVarsName
    );

    return tres;
}


void adjointSpalartAllmaras::nullify()
{
    variablesSet::nullifyField(nuaTilda());
}


void adjointSpalartAllmaras::correct()
{
    addProfiling
    (
        adjointSpalartAllmaras,
        "adjointSpalartAllmaras::correct"
    );

    if (!adjointTurbulence_)
    {
        return;
    }

    adjointTurbulenceModel::correct();

    updatePrimalRelatedFields();

    const surfaceScalarField& phi = primalVars_.phi();
    const volVectorField& Ua = adjointVars_.UaInst();

    volScalarField gradNua(gradNuTilda_ & fvc::grad(nuaTilda()));
    volScalarField gradUaR
    (
        2.0*fvc::grad(Ua) && symmAdjointProductionU_
    );

    dimensionedScalar oneOverSigmaNut = 1./sigmaNut_;

    nuaTilda().storePrevIter();

    fv::options& fvOptions(fv::options::New(mesh_));

    tmp<fvScalarMatrix> nuaTildaEqn
    (
        fvm::ddt(nuaTilda())
      + fvm::div(-phi, nuaTilda())
      - fvm::laplacian(DnuTildaEff(), nuaTilda())
        // Note: Susp
      + fvm::SuSp(productionDestructionSource_, nuaTilda())
      + fvc::laplacian(2.0*Cb2_*oneOverSigmaNut*nuaTilda(), nuTilda())
      + gradNua*oneOverSigmaNut
     ==
        // Always a negative contribution to the lhs. No Sp used!
        Cb1_*Stilda_*nuaTilda()
        // Always a positive contribution to the lhs. no need for SuSp
      - fvm::Sp(Cw1_*fw_*nuTilda()/sqr(y_), nuaTilda())
      - Cdnut_*gradUaR
      + fvOptions(nuaTilda())
    );

    // Add sources from the objective functions
    objectiveManager_.addSource(nuaTildaEqn.ref());

    nuaTildaEqn.ref().relax();
    fvOptions.constrain(nuaTildaEqn.ref());
    solve(nuaTildaEqn);
    nuaTilda().correctBoundaryConditions();
    nuaTilda().relax();

    if (adjointVars_.getSolverControl().printMaxMags())
    {
        dimensionedScalar maxDeltaNuaTilda =
            max(mag(nuaTilda() - nuaTilda().prevIter())());
        dimensionedScalar maxNuaTilda = max(mag(nuaTilda()));
        Info<< "Max mag (" << nuaTilda().name() << ") = "
            << maxNuaTilda.value() << endl;
        Info<< "Max mag (Delta " << nuaTilda().name() << ") = "
            << maxDeltaNuaTilda.value()
            << endl;
    }
}


bool adjointSpalartAllmaras::read()
{
    if (adjointRASModel::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adjointRASModels
} // End namespace incompressibleAdjoint
} // End namespace Foam

// ************************************************************************* //
