/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020,2023 OpenCFD Ltd.
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

#include "kEpsilonPhitF.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonPhitF<BasicTurbulenceModel>::correctNut()
{
    // (LUU:p. 173)
    this->nut_ = Cmu_*phit_*k_*T_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonPhitF<BasicTurbulenceModel>::Ts() const
{
    // (LUU:Eq. 7)
    return
        max
        (
            k_/epsilon_,
            CT_*sqrt
            (
                max
                (
                    this->nu(),
                    dimensionedScalar(this->nu()().dimensions(), Zero)
                )/epsilon_
            )
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonPhitF<BasicTurbulenceModel>::Ls() const
{
    // (LUU:Eq. 7)
    return
        CL_*max
        (
            pow(k_, 1.5)/epsilon_,
            Ceta_*pow025
            (
                pow3
                (
                    max
                    (
                        this->nu(),
                        dimensionedScalar(this->nu()().dimensions(), Zero)
                    )
                )/epsilon_
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonPhitF<BasicTurbulenceModel>::kEpsilonPhitF
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

    includeNu_
    (
        Switch::getOrAddToDict
        (
            "includeNu",
            this->coeffDict_,
            true
        )
    ),
    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    Ceps1a_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1a",
            this->coeffDict_,
            1.4
        )
    ),
    Ceps1b_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1b",
            this->coeffDict_,
            1.0
        )
    ),
    Ceps1c_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1c",
            this->coeffDict_,
            0.05
        )
    ),
    Ceps2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    Cf1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cf1",
            this->coeffDict_,
            1.4
        )
    ),
    Cf2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cf2",
            this->coeffDict_,
            0.3
        )
    ),
    CL_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.25
        )
    ),
    Ceta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            110.0
        )
    ),
    CT_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CT",
            this->coeffDict_,
            6.0
        )
    ),
    sigmaK_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    sigmaPhit_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaPhit",
            this->coeffDict_,
            1.0
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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    phit_
    (
        IOobject
        (
            IOobject::groupName("phit", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    f_
    (
        IOobject
        (
            IOobject::groupName("f", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    T_
    (
        IOobject
        (
            "T",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        this->mesh_,
        dimensionedScalar(dimTime, Zero)
    ),

    phitMin_(dimensionedScalar("phitMin", phit_.dimensions(), SMALL)),
    fMin_(dimensionedScalar("fMin", f_.dimensions(), SMALL)),
    TMin_(dimensionedScalar("TMin", dimTime, SMALL)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), SMALL))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(phit_, phitMin_);
    bound(f_, fMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    if
    (
        mag(sigmaK_.value()) < VSMALL
     || mag(sigmaEps_.value()) < VSMALL
     || mag(sigmaPhit_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmaK = " << sigmaK_ << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << "sigmaPhit = " << sigmaPhit_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonPhitF<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        includeNu_.readIfPresent("includeNu", this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        Ceps1a_.readIfPresent(this->coeffDict());
        Ceps1b_.readIfPresent(this->coeffDict());
        Ceps1c_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Cf1_.readIfPresent(this->coeffDict());
        Cf2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        CT_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaPhit_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kEpsilonPhitF<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volSymmTensorField> tS(symm(fvc::grad(U)));
    volScalarField G(this->GName(), nut*(2.0*(dev(tS()) && tS())));
    tS.clear();

    T_ = Ts();
    bound(T_, TMin_);

    const volScalarField L2(type() + "L2", sqr(Ls()) + L2Min_);

    const volScalarField::Internal Ceps1
    (
        "Ceps1",
        Ceps1a_*(Ceps1b_ + Ceps1c_*sqrt(1.0/phit_()))
    );


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // Turbulent kinetic energy dissipation rate equation (LUU:Eq. 4)
    // k/T ~ epsilon
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
        alpha()*rho()*Ceps1*G()/T_()
      - fvm::SuSp
        (
            (2.0/3.0*Ceps1)*(alpha()*rho()*divU),
            epsilon_
        )
      - fvm::Sp(alpha()*rho()*Ceps2_/T_(), epsilon_)
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation (LUU:Eq. 3)
    // epsilon/k ~ 1/Ts
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
      ==
        alpha()*rho()*G()
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*(1.0/T_()), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // Elliptic relaxation function equation (LUU:Eq. 18)
    // All source terms are non-negative functions (LUU:p. 176)
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
      ==
      - fvm::Sp(1.0/L2(), f_)
      - (
            (Cf1_ - 1.0)*(phit_() - 2.0/3.0)/T_()
           -(Cf2_*G())/k_()
           +(Cf2_*(2.0/3.0)*divU)
           -(2.0*this->nu()*(fvc::grad(phit_) & fvc::grad(k_)))()/k_()
           -(this->nu()*fvc::laplacian(phit_))()
        )/L2()
    );

    fEqn.ref().relax();
    solve(fEqn);
    bound(f_, fMin_);


    // Normalised wall-normal fluctuating velocity scale equation (LUU:Eq. 17)
    // All source terms are non-negative functions (LUU:p. 176)
    tmp<fvScalarMatrix> phitEqn
    (
        fvm::ddt(alpha, rho, phit_)
      + fvm::div(alphaRhoPhi, phit_)
      - fvm::laplacian(alpha*rho*DphitEff(), phit_)
      ==
        alpha()*rho()*f_()
      - fvm::SuSp
        (
            alpha()*rho()*
            (
                G()/k_()
              - (2.0/3.0)*divU
              - (2.0*nut*(fvc::grad(phit_) & fvc::grad(k_)))()
                /(k_()*sigmaPhit_*phit_())
            )
          , phit_
        )
      + fvOptions(alpha, rho, phit_)
    );

    phitEqn.ref().relax();
    fvOptions.constrain(phitEqn.ref());
    solve(phitEqn);
    fvOptions.correct(phit_);
    bound(phit_, phitMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
