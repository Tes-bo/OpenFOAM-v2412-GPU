/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "EBRSM.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> EBRSM<BasicTurbulenceModel>::calcL() const
{
    // (M:Eq. C.13)
    return
        Cl_*max
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


template<class BasicTurbulenceModel>
tmp<volVectorField> EBRSM<BasicTurbulenceModel>::calcN() const
{
    const volVectorField gradf(fvc::grad(f_));

    // (M:Eq. C.9)
    return gradf/max(mag(gradf), dimensionedScalar(dimless/dimLength, SMALL));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> EBRSM<BasicTurbulenceModel>::calcTau() const
{
    // (M:Eq. C.12)
    return
        max
        (
            k_/epsilon_,
            Ct_*sqrt
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
tmp<volSymmTensorField> EBRSM<BasicTurbulenceModel>::D
(
    const volScalarField& tau,
    const dimensionedScalar& sigma
) const
{
    // (M:Eq. C.10, C.14)
    return (Cmu_/sigma*tau)*this->R_ + this->nu()*I;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> EBRSM<BasicTurbulenceModel>::D
(
    const dimensionedScalar& sigma
) const
{
    // (LM:p. 2)
    return this->nut_/sigma + this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> EBRSM<BasicTurbulenceModel>::Ceps1Prime
(
    const volScalarField& G
) const
{
    // (M:Eq. C.15)
    return Ceps1_*(scalar(1) + A1_*(scalar(1) - pow3(f_))*G/this->epsilon_);
}


template<class BasicTurbulenceModel>
void EBRSM<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*k_*calcTau();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EBRSM<BasicTurbulenceModel>::EBRSM
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
    ReynoldsStress<RASModel<BasicTurbulenceModel>>
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

    simpleGradientDiffusion_
    (
        Switch::getOrAddToDict
        (
            "simpleGradientDiffusion",
            this->coeffDict_,
            false
        )
    ),
    g1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g1",
            this->coeffDict_,
            3.4
        )
    ),
    g1star_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g1star",
            this->coeffDict_,
            1.8
        )
    ),
    g3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g3",
            this->coeffDict_,
            0.8
        )
    ),
    g3star_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g3star",
            this->coeffDict_,
            1.3
        )
    ),
    g4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g4",
            this->coeffDict_,
            1.25
        )
    ),
    g5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "g5",
            this->coeffDict_,
            0.2
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.21
        )
    ),
    Ceps1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.83
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.15
        )
    ),
    A1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1",
            this->coeffDict_,
            0.065
        )
    ),
    Ct_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct",
            this->coeffDict_,
            6.0
        )
    ),
    Cl_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cl",
            this->coeffDict_,
            0.133
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            80.0
        )
    ),
    Cstability_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cstability",
            this->coeffDict_,
            sqr(dimLength)/pow3(dimTime),
            10.0
        )
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
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
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
    )
{
    this->boundNormalStress(this->R_);
    bound(epsilon_, this->epsilonMin_);
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool EBRSM<BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicTurbulenceModel>>::read())
    {
        simpleGradientDiffusion_.readIfPresent
        (
            "simpleGradientDiffusion",
            this->coeffDict()
        );
        g1_.readIfPresent(this->coeffDict());
        g1star_.readIfPresent(this->coeffDict());
        g3_.readIfPresent(this->coeffDict());
        g3star_.readIfPresent(this->coeffDict());
        g4_.readIfPresent(this->coeffDict());
        g5_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        A1_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());
        Cl_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Cstability_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void EBRSM<BasicTurbulenceModel>::correct()
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
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    ReynoldsStress<RASModel<BasicTurbulenceModel>>::correct();


    // Calculate the velocity gradient tensor in Hessian form (delta_i u_j)
    // Tranpose of the classical Jacobian form (delta_j u_i)
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    // Calculate the production tensor
    tmp<volSymmTensorField> tP = -twoSymm(R & gradU);
    const volSymmTensorField& P = tP.cref();

    // Calculate turbulent kinetic energy production rate
    const volScalarField G(this->GName(), 0.5*mag(tr(P)));


    // Calculate elliptic blending function
    // between near-wall and weakly-inhomogeneous regions
    {
        // (M:Eq. C.13)
        tmp<volScalarField> tinvLsqr = scalar(1)/sqr(calcL());
        const volScalarField& invLsqr = tinvLsqr.cref();

        // (M:Eq. C.8)
        tmp<fvScalarMatrix> fEqn
        (
            fvm::Sp(invLsqr, f_)
          - fvm::laplacian(f_)
          ==
            invLsqr
        );

        tinvLsqr.clear();

        fEqn.ref().relax();
        solve(fEqn);
    }


    // Calculate approximate wall-normal vector field (M:Eq. C.9)
    tmp<volVectorField> tn = calcN();
    const volVectorField& n = tn.cref();

    // Calculate turbulent time scale (M:Eq. C.12)
    tmp<volScalarField> ttau = calcTau();
    const volScalarField& tau = ttau.cref();


    // Calculate turbulent dissipation rate field
    {
        // Dissipation-production stimulator in the buffer layer (M:Eq. C.15)
        tmp<volScalarField> tCeps1Prime = Ceps1Prime(G);
        const volScalarField& Ceps1Prime = tCeps1Prime.cref();

        // Update epsilon and G at the wall
        epsilon_.boundaryFieldRef().updateCoeffs();
        // Push any changed cell values to coupled neighbours
        epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

        // (M:Eq. C.14)
        tmp<fvScalarMatrix> epsEqn
        (
            fvm::ddt(alpha, rho, epsilon_)
          + fvm::div(alphaRhoPhi, epsilon_)
          - (
                simpleGradientDiffusion_
              ? fvm::laplacian(alpha*rho*D(sigmaEps_), epsilon_)
              : fvm::laplacian(alpha*rho*D(tau, sigmaEps_), epsilon_)
            )
          + fvm::Sp(Cstability_/k_, epsilon_)
          ==
            alpha()*rho()*Ceps1Prime()*G()/tau()
          - fvm::Sp(alpha*rho*Ceps2_/tau, epsilon_)
          + Cstability_*epsilon_()/k_()
          + fvOptions(alpha, rho, epsilon_)
        );

        tCeps1Prime.clear();

        epsEqn.ref().relax();
        fvOptions.constrain(epsEqn.ref());
        epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
        solve(epsEqn);
        fvOptions.correct(epsilon_);

        bound(epsilon_, this->epsilonMin_);
    }


    // Calculate Reynolds-stress field
    {
        // Homogeneous form of the redistribution term (M:Eq. C.3)
        tmp<volSymmTensorField> tPhiH;
        {
            // Reynolds stress anisotropy tensor (M:Eq. C.4)
            const volSymmTensorField B(R/(scalar(2)*k_) - oneThirdI);

            // Rate-of-strain tensor (M:Eq. C.5)
            const volSymmTensorField S(symm(gradU));

            // Rate-of-rotation tensor (M:Eq. C.6)
            // Note the Hessian form of gradient
            const volTensorField W(gradU.T() - gradU);

            tPhiH =
                k_
               *(
                    (g3_ - g3star_*mag(B))*S
                  + g4_*devTwoSymm(B & S)
                  + g5_*twoSymm(B & W.T())
                );
        }


        // Near-wall form of the redistribution model (M:Eq. C.7)
        tmp<volSymmTensorField> tPhiW;
        {
            tmp<volSymmTensorField> tnn = symm(n*n);
            const volSymmTensorField& nn = tnn.cref();

            tn.clear();

            tPhiW =
              - scalar(5)*epsilon_/k_*
                (
                    twoSymm(R & nn)
                  - 0.5*(R && nn)*(nn + I)
                );
        }


        tmp<fvSymmTensorMatrix> REqn;
        {
            const volScalarField fCube(pow3(f_));

            // Velocity-pressure gradient correlation (M:Eq. C.2)
            const volSymmTensorField Phi
            (
                (scalar(1) - fCube)*tPhiW + fCube*tPhiH
            );

            // Near-wall part of the dissipation tensor (M:Eq. C.11)
            const volScalarField epsilonW
            (
                (scalar(1) - fCube)*epsilon_/k_
            );

            // Homogeneous part of the dissipation tensor (M:Eq. C.11)
            const volSphericalTensorField epsilonH
            (
                fCube*epsilon_*twoThirdsI
            );

            // Implicit part of the g1-term (M:Eq. C.3)
            const volScalarField Phi1Implicit
            (
                fCube*(g1_*epsilon_ + g1star_*G)/(scalar(2)*k_)
            );

            // Explicit part of the g1-term (M:Eq. C.3)
            const volSphericalTensorField Phi1Explicit
            (
                fCube*(g1_*epsilon_ + g1star_*G)*oneThirdI
            );


            // (M:Eq. C.1)
            REqn =
            (
                fvm::ddt(alpha, rho, R)
              + fvm::div(alphaRhoPhi, R)
              - (
                    simpleGradientDiffusion_
                  ? fvm::laplacian(alpha*rho*D(sigmaK_), R)
                  : fvm::laplacian(alpha*rho*D(tau, sigmaK_), R)
                )
              + fvm::Sp(alpha*rho*Phi1Implicit, R)
              + fvm::Sp(alpha*rho*epsilonW, R)
              ==
                alpha()*rho()*
                (
                    P()
                  + Phi()
                  + Phi1Explicit()
                  - epsilonH()
                )
              + fvOptions(alpha, rho, R)
            );

            ttau.clear();
            tP.clear();
        }


        REqn.ref().relax();
        fvOptions.constrain(REqn.ref());
        solve(REqn);
        fvOptions.correct(R);

        this->boundNormalStress(R);

        #ifdef FULLDEBUG
        this->checkRealizabilityConditions(R);
        #endif

        k_ == 0.5*tr(R);
        bound(k_, this->kMin_);
    }

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
