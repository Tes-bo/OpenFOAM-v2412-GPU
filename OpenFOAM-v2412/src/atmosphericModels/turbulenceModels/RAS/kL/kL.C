/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "kL.H"
#include "fvOptions.H"
#include "bound.H"
#include "gravityMeshObject.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::Cmu() const
{
    // (A:Eq. 31)
    return (0.556 + 0.108*Rt_)/(1.0 + 0.308*Rt_ + 0.00837*sqr(Rt_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::CmuPrime() const
{
    // (A:Eq. 32)
    return 0.556/(1.0 + 0.277*Rt_);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::nutPrime() const
{
    // (A:Eq. 12)
    return CmuPrime()*sqrt(k_)*L_;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::epsilonCanopy() const
{
    const auto* CdPtr = this->mesh_.template findObject<volScalarField>("Cd");
    const auto* LADPtr = this->mesh_.template findObject<volScalarField>("LAD");
    const volVectorField& U = this->U_;

    if (CdPtr && LADPtr)
    {
        const auto& Cd = *CdPtr;
        const auto& LAD = *LADPtr;

        // (W:Eq. 13)
        return Cd*LAD*mag(U)*k_;
    }

    return volScalarField::New
    (
        IOobject::groupName("epsilonCanopy", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::epsilon() const
{
    // (W:Eq. 13)
    tmp<volScalarField> tepsilonCanopy = epsilonCanopy();

    // (A:Eq. 19)
    tmp<volScalarField> tepsilonPlain = pow3(Cmu0_)*pow(k_, 1.5)/L_;

    // (W:Eq. 13)
    tmp<volScalarField> tepsilon = max(tepsilonPlain, tepsilonCanopy);
    volScalarField& epsilon = tepsilon.ref();
    bound(epsilon, this->epsilonMin_);

    return tepsilon;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::canopyHeight() const
{
    tmp<volScalarField> tcanopy;

    tcanopy.cref
    (
        this->mesh_.template cfindObject<volScalarField>("canopyHeight")
    );

    if (tcanopy)
    {
        return tcanopy;
    }

    return volScalarField::New
    (
        IOobject::groupName("canopyHeight", this->alphaRhoPhi_.group()),
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimLength, Zero)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kL<BasicTurbulenceModel>::L() const
{
    // (A:Eq. 22)
    const volScalarField Lplain(kappa_*y_);

    // Increase roughness for canopy (forest, vegetation etc)
    tmp<volScalarField> tLcanopy = kappa_*canopyHeight();
    const volScalarField& Lcanopy = tLcanopy;

    // (W:Eq. 16)
    return max(Lcanopy, Lplain);
}


template<class BasicTurbulenceModel>
void kL<BasicTurbulenceModel>::stratification(const volScalarField& fVB)
{
    tmp<volScalarField> tLg = L();
    const volScalarField& Lg = tLg.cref();

    tmp<volScalarField> tcanopyHeight = canopyHeight();
    const volScalarField& canopyHeight = tcanopyHeight;

    tmp<volScalarField> tLcanopy = kappa_*canopyHeight;
    const volScalarField& Lcanopy = tLcanopy;

    const scalar Cmu0 = Cmu0_.value();
    const scalar CbStable = CbStable_.value();
    const scalar CbUnstable = CbUnstable_.value();

    forAll(L_, i)
    {
        if (y_[i] > canopyHeight[i])
        {
            if (fVB[i] > 0)
            {
                // (A:Eq. 23)
                const scalar Lb = CbStable*sqrt(k_[i])/sqrt(fVB[i]);

                // (A:Eq. 26)
                L_[i] = sqrt(sqr(Lg[i]*Lb)/(sqr(Lg[i]) + sqr(Lb)));
            }
            else
            {
                // For unstable/neutral boundary layer (A:p. 80)
                // Smoothing function for turbulent Richardson
                // number to ensure gentle transition into
                // the regime of strong convection
                Rt_[i] =
                    min
                    (
                        max(Rt_[i], -1.0),
                        Rt_[i] - sqr(Rt_[i] + 1.0)/(Rt_[i] - 1.0)
                    );

                // (A:Eq. 28)
                L_[i] =
                    Lg[i]
                   *sqrt(1.0 - pow(Cmu0, 6.0)*pow(CbUnstable, -2.0)*Rt_[i]);
            }
        }
        else
        {
            L_[i] = Lcanopy[i];
        }
    }

    // Limit characteristic length scale
    L_ = min(L_, Lmax_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kL<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu()*sqrt(k_)*L_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kL<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kL<BasicTurbulenceModel>::kL
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
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
    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            dimless/dimTemperature,
            3.3e-03
        )
    ),
    Cmu0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu0",
            this->coeffDict_,
            0.556
        )
    ),
    Lmax_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Lmax",
            this->coeffDict_,
            dimLength,
            GREAT
        )
    ),
    CbStable_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CbStable",
            this->coeffDict_,
            0.25
        )
    ),
    CbUnstable_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CbUnstable",
            this->coeffDict_,
            0.35
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    L_
    (
        IOobject
        (
            IOobject::groupName("L", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimLength, scalar(1))
    ),
    Rt_
    (
        IOobject
        (
            IOobject::groupName("Rt", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    g_(meshObjects::gravity::New(this->mesh_.time())),
    y_(wallDist::New(this->mesh_).y())
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kL<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        kappa_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        Cmu0_.readIfPresent(this->coeffDict());
        Lmax_.readIfPresent(this->coeffDict());
        CbStable_.readIfPresent(this->coeffDict());
        CbUnstable_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kL<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    // Turbulent kinetic energy production rate
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal G
    (
        this->GName(),
        nut.v()*2*magSqr(devSymm(tgradU.cref().v()))
    );
    tgradU.clear();

    // Square of Brunt-Vaisala (buoyancy) frequency
    const auto& T = U.mesh().lookupObject<volScalarField>("T");
    tmp<volScalarField> tfBV = -beta_*(fvc::grad(T) & g_);
    const volScalarField& fBV = tfBV.cref();

    // Sink or source of TKE depending on stratification type (A:Eq. 15)
    tmp<volScalarField> tPb = -fBV*nutPrime();
    const volScalarField& Pb = tPb.cref();

    // Turbulent kinetic energy dissipation rate due to plains and canopy
    tmp<volScalarField> tepsilon = epsilon();
    const volScalarField& epsilon = tepsilon.cref();

    // Divergence of velocity
    tmp<volScalarField> tdivU = fvc::div(fvc::absolute(this->phi(), U));
    const volScalarField::Internal& divU = tdivU.cref().v();

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      + fvm::SuSp((Pb - epsilon)/k_, k_)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    tdivU.clear();
    tPb.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // Turbulent Richardson number (A:Eq. 29)
    Rt_ = fBV*sqr(k_/tepsilon);

    stratification(fBV);
    tfBV.clear();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
