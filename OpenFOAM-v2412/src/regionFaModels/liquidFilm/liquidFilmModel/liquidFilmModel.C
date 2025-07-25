/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "liquidFilmModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmModel, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void liquidFilmModel::correctThermoFields()
{
    scalarField X(thermo_.size(), 1);

    forAll(rho_, faceI)
    {
        rho_[faceI] = thermo_.rho(pRef_, Tf_[faceI], X);
        mu_[faceI] = thermo_.mu(pRef_, Tf_[faceI], X);
        sigma_[faceI] = thermo_.sigma(pRef_, Tf_[faceI], X);
        Cp_[faceI] = thermo_.Cp(pRef_, Tf_[faceI], X);
    }

    forAll(regionMesh().boundary(), patchI)
    {
        const scalarField& patchTf = Tf_.boundaryFieldRef()[patchI];

        scalarField& patchRho = rho_.boundaryFieldRef()[patchI];
        scalarField& patchmu = mu_.boundaryFieldRef()[patchI];
        scalarField& patchsigma = sigma_.boundaryFieldRef()[patchI];
        scalarField& patchCp = Cp_.boundaryFieldRef()[patchI];

        forAll(patchRho, edgeI)
        {
            patchRho[edgeI] = thermo_.rho(pRef_, patchTf[edgeI], X);
            patchmu[edgeI] = thermo_.mu(pRef_, patchTf[edgeI], X);
            patchsigma[edgeI] = thermo_.sigma(pRef_, patchTf[edgeI], X);
            patchCp[edgeI] = thermo_.Cp(pRef_, patchTf[edgeI], X);
        }
    }

    //Initialize pf_
    pf_ = rho_*gn_*h_ - sigma_*fac::laplacian(h_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmModel::liquidFilmModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    liquidFilmBase(modelType, mesh, dict),
    thermo_(dict.subDict("thermo")),
    rho_
    (
        IOobject
        (
            "rhof",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, Zero)
    ),
    mu_
    (
        IOobject
        (
            "muf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimViscosity, Zero)
    ),
    Tf_
    (
        IOobject
        (
            "Tf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero)
    ),
    Cp_
    (
        IOobject
        (
            "Cp_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTemperature, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero)
    ),
    hRho_
    (
        IOobject
        (
            h_.name() + "*" + rho_.name(),
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(h_.dimensions()*rho_.dimensions(), Zero)
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimVelocity, Zero)
    ),
    USp_
    (
        IOobject
        (
            "USp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(sqr(dimVelocity), Zero)
    ),
    pnSp_
    (
        IOobject
        (
            "pnSp",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    cloudMassTrans_
    (
        IOobject
        (
            "cloudMassTrans",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, Zero)
    ),
    cloudDiameterTrans_
    (
        IOobject
        (
            "cloudDiameterTrans",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    turbulence_(filmTurbulenceModel::New(*this, dict)),
    availableMass_(regionMesh().faces().size(), Zero),
    injection_(*this, dict),
    forces_(*this, dict)
{
    if (dict.readIfPresent("T0", Tref_))
    {
        Info<< "Initialised film temperature to T0" << endl;

        Tf_ = dimensionedScalar("T0", dimTemperature, dict);
    }
    correctThermoFields();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& liquidFilmModel::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const volScalarField& liquidFilmModel::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}


void liquidFilmModel::preEvolveRegion()
{
    liquidFilmBase::preEvolveRegion();


    cloudMassTrans_ == dimensionedScalar(dimMass, Zero);
    cloudDiameterTrans_ == dimensionedScalar(dimLength, Zero);

    const scalar deltaT = primaryMesh().time().deltaTValue();
    const scalarField rAreaDeltaT(scalar(1)/deltaT/regionMesh().S().field());

    // Map the total mass, mom [kg.m/s] and pnSource from particles

    vsm().mapToSurface(massSource_, rhoSp_.primitiveFieldRef());

    vsm().mapToSurface(momentumSource_, USp_.primitiveFieldRef());

    vsm().mapToSurface(pnSource_, pnSp_.primitiveFieldRef());


    // Calculate rate per area
    rhoSp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    USp_.primitiveFieldRef() *= rAreaDeltaT/rho_;
    pnSp_.primitiveFieldRef() *= rAreaDeltaT/rho_;

    rhoSp_.relax();
    pnSp_.relax();
    USp_.relax();
}


void liquidFilmModel::postEvolveRegion()
{
    availableMass_ = (h() - h0_)*rho()*regionMesh().S();
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);
    liquidFilmBase::postEvolveRegion();
}


void liquidFilmModel::info()
{
    Info<< "\nSurface film: " << type() << " on patch: ";

    for (const label patchi : this->primaryPatchIDs())
    {
        Info<< ' ' << patchi;
    }
    Info<< endl;

    const DimensionedField<scalar, areaMesh>& sf = regionMesh().S();

    Info<< indent << "min/max(mag(Uf))    = "
        << gMinMaxMag(Uf_.field()) << nl
        << indent << "min/max(delta)     = "
        << gMinMax(h_.field()) << nl
        << indent << "coverage           = "
        << gSum(alpha()().field()*mag(sf.field()))/gSumMag(sf.field()) <<  nl
        << indent << "total mass         = "
        << gSum(availableMass_) << nl;

    injection_.info(Info);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
