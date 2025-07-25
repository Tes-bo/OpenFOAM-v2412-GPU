/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "thermoSingleLayer.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || isA<mixedFvPatchScalarField>(T_.boundaryField()[patchi])
         || isA<mappedFieldFvPatchField<scalar>>(T_.boundaryField()[patchi])
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), Zero);
}


void thermoSingleLayer::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();
}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            hsBf[patchi] == hs(Tp, patchi);
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchi];
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();
}


void thermoSingleLayer::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}


void thermoSingleLayer::updateSubmodels()
{
    DebugInFunction << endl;

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_)/(2*rhoPrimary_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Update source fields
    rhoSp_ += rMagSfDt*(cloudMassTrans_ + primaryMassTrans_);
    hsSp_ += rMagSfDt*(cloudMassTrans_*hs_ + primaryEnergyTrans_);

    turbulence_->correct();
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{
    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp(htcs_->h()/Cp_, hs)
      + htcs_->h()*(hs/Cp_ + alpha_*(TPrimary_ - T_))

        // Heat-transfer to the wall
      - fvm::Sp(htcw_->h()/Cp_, hs)
      + htcw_->h()*(hs/Cp_ + alpha_*(Tw_- T_))
    );
}


void thermoSingleLayer::solveEnergy()
{
    DebugInFunction << endl;

    dimensionedScalar residualDeltaRho
    (
        "residualDeltaRho",
        deltaRho_.dimensions(),
        1e-10
    );

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
      - hsSp_
      + q(hs_)
      + radiation_->Shs()
    );

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);

    // Update film wall and surface temperatures
    updateSurfaceTemperatures();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    Cp_
    (
        IOobject
        (
            "Cp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    T_
    (
        IOobject
        (
            "Tf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        fvPatchFieldBase::zeroGradientType()
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        fvPatchFieldBase::zeroGradientType()
    ),
    hs_
    (
        IOobject
        (
            "hf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero),
        hsBoundaryTypes()
    ),

    primaryEnergyTrans_
    (
        IOobject
        (
            "primaryEnergyTrans",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    deltaWet_(coeffs_.get<scalar>("deltaWet")),
    hydrophilic_(coeffs_.get<bool>("hydrophilic")),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, Zero),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // Must have same name as hSp_ to enable mapping
            primaryMesh().time().timeName(),
            primaryMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero)
    ),

    TPrimary_
    (
        IOobject
        (
            "T", // Same name as T on primary region to enable mapping
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    YPrimary_(),

    viscosity_(filmViscosityModel::New(*this, coeffs(), mu_)),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    radiation_(filmRadiationModel::New(*this, coeffs())),
    withTbounds_(limitType::CLAMP_NONE),
    Tbounds_(0, 5000)
{
    unsigned userLimits(limitType::CLAMP_NONE);

    if (coeffs().readIfPresent("Tmin", Tbounds_.min()))
    {
        userLimits |= limitType::CLAMP_MIN;
        Info<< "    limiting minimum temperature to " << Tbounds_.min() << nl;
    }

    if (coeffs().readIfPresent("Tmax", Tbounds_.max()))
    {
        userLimits |= limitType::CLAMP_MAX;
        Info<< "    limiting maximum temperature to " << Tbounds_.max() << nl;
    }
    withTbounds_ = limitType(userLimits);

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        regionMesh().time().timeName(),
                        regionMesh().thisDb(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, Zero),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.readEntry("hydrophilicDryScale", hydrophilicDryScale_);
        coeffs_.readEntry("hydrophilicWetScale", hydrophilicWetScale_);
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_);

        deltaRho_ == delta_*rho_;

        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                regionMesh().time().timeName(),
                regionMesh().thisDb(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    DebugInfo
        << "    energy   = " << energySource << nl << nl;

    hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    DebugInFunction << endl;

    kinematicSingleLayer::preEvolveRegion();
    primaryEnergyTrans_ == dimensionedScalar(dimEnergy, Zero);
}


void thermoSingleLayer::evolveRegion()
{
    DebugInFunction << endl;

    // Solve continuity for deltaRho_
    solveContinuity();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> tUEqn = solveMomentum(tpu(), tpp());
        fvVectorMatrix& UEqn = tUEqn.ref();

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn);
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
    T_ == T(hs_);
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = T_;

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho() const
{
    auto tSrho = volScalarField::Internal::New
    (
        IOobject::scopedName(typeName, "Srho"),
        IOobject::NO_REGISTER,
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
    scalarField& Srho = tSrho.ref();

    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];
        const label primaryPatchi = primaryPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMass);

        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    auto tSrho = volScalarField::Internal::New
    (
        IOobject::scopedName(typeName, "Srho(" + Foam::name(i) + ")"),
        IOobject::NO_REGISTER,
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
    scalarField& Srho = tSrho.ref();

    if (vapId == i)
    {
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];
            const label primaryPatchi = primaryPatchIDs()[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Sh() const
{
    auto tSh = volScalarField::Internal::New
    (
        IOobject::scopedName(typeName, "Sh"),
        IOobject::NO_REGISTER,
        primaryMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
    scalarField& Sh = tSh.ref();

    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];
        const label primaryPatchi = primaryPatchIDs()[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace surfaceFilmModels

// ************************************************************************* //
