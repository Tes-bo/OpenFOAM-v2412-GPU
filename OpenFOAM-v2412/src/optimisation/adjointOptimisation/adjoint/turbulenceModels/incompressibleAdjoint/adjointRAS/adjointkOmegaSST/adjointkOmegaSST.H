/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022 PCOpt/NTUA
    Copyright (C) 2014-2022 FOSS GP
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


Class
    Foam::incompressibleAdjoint::adjointRASModels::adjointkOmegaSST

Description
    Continuous adjoint to the kOmegaSST turbulence model for incompressible
    flows.

    Reference:
    \verbatim
        The code is based on the following reference, with a number of
        changes in the numerical implementation

            Kavvadias, I., Papoutsis-Kiachagias, E., Dimitrakopoulos, G., &
            Giannakoglou, K. (2014).
            The continuous adjoint approach to the k–ω SST turbulence model
            with applications in shape optimization
            Engineering Optimization, 47(11), 1523-1542.
            https://doi.org/10.1080/0305215X.2014.979816
    \endverbatim

SourceFiles
    adjointkOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adjointkOmegaSST_H
#define Foam_adjointkOmegaSST_H

#include "adjointRASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{
namespace adjointRASModels
{

/*---------------------------------------------------------------------------*\
                       Class adjointkOmegaSST Declaration
\*---------------------------------------------------------------------------*/

class adjointkOmegaSST
:
    public adjointRASModel
{
    // Private Member Functions

        //- No copy construct
        adjointkOmegaSST(const adjointkOmegaSST&) = delete;

        //- No copy assignment
        void operator=(const adjointkOmegaSST&) = delete;


protected:

    // Protected data

        // Primal Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            //- Flag to include the F3 term
            Switch F3_;


        // Fields

            // Primal Fields

                //- Wall distance
                //  Note: reference to the distance known by the primal model
                const volScalarField& y_;

                //- Cached primal gradient fields
                volTensorField gradU_;
                volVectorField gradOmega_;
                volVectorField gradK_;

                //- Primal cached fields involved in the solution of the
                // adjoint equations
                // Cached to reduce the computational cost
                volScalarField S2_;
                volScalarField S_;
                volScalarField GbyNu0_;
                volScalarField CDkOmega_;
                volScalarField CDkOmegaPlus_;
                volScalarField F1_;
                volScalarField F2_;

                // Model Field coefficients
                volScalarField alphaK_;
                volScalarField alphaOmega_;
                volScalarField beta_;
                volScalarField gamma_;


         // Switch fields

            // Switch fields for the differentiation of F1
            volScalarField case_1_F1_;
            volScalarField case_2_F1_;
            volScalarField case_3_F1_;
            volScalarField case_4_F1_;

            //- Switch fields for the production in the k Eqn
            volScalarField case_1_Pk_;
            volScalarField case_2_Pk_;
            volScalarField case_3_Pk_;

            // Switch fields for the differentiation of nut
            // Holds also for the differentiation of the second branch of
            // GbyNu
            volScalarField case_1_nut_;
            volScalarField case_2_nut_;
            volScalarField case_3_nut_;

            // Switch fields for GPrime
            volScalarField case_1_GPrime_;
            volScalarField case_2_GPrime_;

            // Zero first cell field and IDs
            // Since the omega equation is a two-zonal one, some
            // of the terms in the adjoint equations need to ba canceled
            // at the cells next to omegaWallFunction BCs
            labelList firstCellIDs_;
            volScalarField zeroFirstCell_;


        // Turbulence model multipliers

            //- Nut Jacobian w.r.t. omega
            volScalarField dnut_domega_;

            //- Nut Jacobian w.r.t. k
            volScalarField dnut_dk_;

            //- Diffusivity of the omega equation
            volScalarField DOmegaEff_;

            //- Diffusivity of the k equation
            volScalarField DkEff_;


    // Protected Member Functions

        // Primal functions

            virtual tmp<volScalarField> F1() const;
            virtual tmp<volScalarField> F2() const;
            virtual tmp<volScalarField> GbyNu
            (
                const volScalarField& GbyNu0,
                const volScalarField& F2,
                const volScalarField& S2
            ) const;

            //- Return G/nu
            virtual tmp<volScalarField::Internal> GbyNu
            (
                const volScalarField::Internal& GbyNu0,
                const volScalarField::Internal& F2,
                const volScalarField::Internal& S2
            ) const;

            tmp<volScalarField> blend
            (
                const volScalarField& F1,
                const dimensionedScalar& psi1,
                const dimensionedScalar& psi2
            ) const
            {
                return F1*(psi1 - psi2) + psi2;
            }

            tmp<volScalarField::Internal> blend
            (
                const volScalarField::Internal& F1,
                const dimensionedScalar& psi1,
                const dimensionedScalar& psi2
            ) const
            {
                return F1*(psi1 - psi2) + psi2;
            }

            tmp<volScalarField> alphaK(const volScalarField& F1) const
            {
                return blend(F1, alphaK1_, alphaK2_);
            }

            tmp<volScalarField> alphaOmega(const volScalarField& F1) const
            {
                return blend(F1, alphaOmega1_, alphaOmega2_);
            }

            tmp<volScalarField::Internal> beta
            (
                const volScalarField::Internal& F1
            ) const
            {
                return tmp<volScalarField::Internal>::New
                (
                    IOobject::scopedName(this->type(), "beta"),
                    blend(F1, beta1_, beta2_)
                );
            }

            tmp<volScalarField> beta
            (
                const volScalarField& F1
            ) const
            {
                return tmp<volScalarField>::New
                (
                    IOobject::scopedName(this->type(), "beta"),
                    blend(F1, beta1_, beta2_)
                );
            }

            tmp<volScalarField::Internal> gamma
            (
                const volScalarField::Internal& F1
            ) const
            {
                return tmp<volScalarField::Internal>::New
                (
                    IOobject::scopedName(this->type(), "gamma"),
                    blend(F1, gamma1_, gamma2_)
                );
            }

            tmp<volScalarField> gamma
            (
                const volScalarField& F1
            ) const
            {
                return tmp<volScalarField>::New
                (
                    IOobject::scopedName(this->type(), "gamma"),
                    blend(F1, gamma1_, gamma2_)
                );
            }

            tmp<volScalarField> zeroFirstCell();


            //  References to the primal turbulence model variables

                inline const volScalarField& k() const
                {
                    return primalVars_.RASModelVariables()().TMVar1();
                }

                inline volScalarField& k()
                {
                    return primalVars_.RASModelVariables()().TMVar1();
                }

                inline const volScalarField& omega() const
                {
                    return primalVars_.RASModelVariables()().TMVar2();
                }

                inline volScalarField& omega()
                {
                    return primalVars_.RASModelVariables()().TMVar2();
                }

                inline const volScalarField& nutRef() const
                {
                    return primalVars_.RASModelVariables()().nutRef();
                }

                inline volScalarField& nutRef()
                {
                    return primalVars_.RASModelVariables()().nutRef();
                }


        // Adjoint related protected member functions

            //- Derivative of the primal equations wrt nut
            tmp<volScalarField> dR_dnut();

            //- Nut Jacobian wrt omega
            tmp<volScalarField> dnut_domega
            (
                const volScalarField& F2,
                const volScalarField& S,
                const volScalarField& case_1_nut,
                const volScalarField& case_2_nut,
                const volScalarField& case_3_nut
            ) const;

            //- Nut Jacobian wrt k
            tmp<volScalarField> dnut_dk
            (
                const volScalarField& F2,
                const volScalarField& S,
                const volScalarField& case_2_nut
            ) const;

            //- F2 Jacobian wrt omega
            tmp<volScalarField> dF2_domega
            (
                const volScalarField& F2,
                const volScalarField& case_2_nut,
                const volScalarField& case_3_nut
            ) const;

            //- F2 Jacobian wrt k
            tmp<volScalarField> dF2_dk
            (
                const volScalarField& F2,
                const volScalarField& case_2_nut
            ) const;

            //- GbyNu Jacobian wrt omega
            tmp<volScalarField> dGPrime_domega() const;

            //- GbyNu Jacobian wrt k
            tmp<volScalarField> dGPrime_dk() const;

            //- Derivative of the primal equations wrt F1
            tmp<volScalarField> dR_dF1() const;

            //- F1 Jacobian wrt omega (no contributions from grad(omega))
            tmp<volScalarField> dF1_domega(const volScalarField& arg1) const;

            //- F1 Jacobian wrt grad(omega)
            tmp<volVectorField> dF1_dGradOmega
            (
                const volScalarField& arg1
            ) const;

            //- Source to waEqn from the differentiation of F1
            tmp<volScalarField> waEqnSourceFromF1() const;

            //- Source to waEqn from the differentiation of CDkOmega
            tmp<fvScalarMatrix> waEqnSourceFromCDkOmega() const;

            //- F1 Jacobian wrt k (no contributions from grad(k))
            tmp<volScalarField> dF1_dk(const volScalarField& arg1) const;

            //- F1 Jacobian wrt grad(k)
            tmp<volVectorField> dF1_dGradK(const volScalarField& arg1) const;

            //- Source to kaEqn from the differentiation of F1
            tmp<volScalarField> kaEqnSourceFromF1() const;

            //- Source to kaEqn from the differentiation of CDkOmega
            tmp<volScalarField> kaEqnSourceFromCDkOmega() const;

            //- Differentiation of the turbulence model diffusion coefficients
            tmp<volScalarField> coeffsDifferentiation
            (
                const volScalarField& primalField,
                const volScalarField& adjointField,
                const word& schemeName
            ) const;

            //- Term multiplying dnut/db, coming from the turbulence model
            tmp<volScalarField> dNutdbMult
            (
                const volScalarField& primalField,
                const volScalarField& adjointField,
                const volScalarField& coeffField,
                const volScalarField& bcField,
                const word& schemeName
            ) const;

            //- Term multiplying dnut/db, coming from the momentum equations
            tmp<volScalarField> dNutdbMult
            (
                const volVectorField& primalField,
                const volVectorField& adjointField,
                const volScalarField& bcField,
                const word& schemeName
            ) const;

            // Functions computing the adjoint mean flow source

                //- Contributions from the turbulence model convection terms
                tmp<volVectorField> convectionMeanFlowSource
                (
                    const volScalarField& primalField,
                    const volScalarField& adjointField
                ) const;

                //- Contributions from the G
                tmp<volVectorField> GMeanFlowSource
                (
                    tmp<volSymmTensorField>& GbyNuMult
                ) const;

                //- Contributions from the divU
                tmp<volVectorField> divUMeanFlowSource
                (
                    tmp<volScalarField>& divUMult
                ) const;

                //- Contributions from nut(U), in the diffusion coefficients
                //- of the turbulence model
                tmp<volScalarField> diffusionNutMeanFlowMult
                (
                    const volScalarField& primalField,
                    const volScalarField& adjointField,
                    const volScalarField& coeffField
                ) const;

                //- Contributions from nut(U)
                tmp<volVectorField> nutMeanFlowSource
                (
                    tmp<volScalarField>& mult,
                    const volScalarField& F2,
                    const volScalarField& S,
                    const volScalarField& case_1_nut,
                    const volTensorField& gradU
                ) const;


            //- Contributions from the differentiation of k existing in
            //- nutkWallFunction.
            //  This could also be implemented in kaqRWallFunction but all
            //  the fields required for the computation already exist here,
            //  hence the code complexity is reduced
            void addWallFunctionTerms
            (
                fvScalarMatrix& kaEqn,
                const volScalarField& dR_dnut
            );

            // References to the adjoint turbulence model fields

                inline volScalarField& ka()
                {
                    return adjointTMVariable1Ptr_();
                };

                inline const volScalarField& ka() const
                {
                    return adjointTMVariable1Ptr_();
                };

                inline volScalarField& wa()
                {
                    return adjointTMVariable2Ptr_();
                };

                inline const volScalarField& wa() const
                {
                    return adjointTMVariable2Ptr_();
                };


            //- Update of the primal cached fields
            void updatePrimalRelatedFields();

            //- Return the requested interpolation scheme if it exists,
            //- otherwise return a reverseLinear scheme
            template<class Type>
            tmp<surfaceInterpolationScheme<Type>> interpolationScheme
            (
                const word& schemeName
            ) const;

            //- Return the interpolation scheme used by the primal convection
            //- term of the equation corresponding to the argument
            tmp<surfaceInterpolationScheme<scalar>> convectionScheme
            (
                const word& varName
            ) const;


public:

    //- Runtime type information
    TypeName("adjointkOmegaSST");


    // Constructors

        //- Construct from components
        adjointkOmegaSST
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName
                = adjointTurbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~adjointkOmegaSST() = default;


    // Member Functions

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK(F1)*this->nut() + this->nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega(F1)*this->nut() + this->nu()
                )
            );
        }

        virtual tmp<volSymmTensorField> devReff() const;

        virtual tmp<volSymmTensorField> devReff(const volVectorField& U) const;

        //- Return the transpose part of the adjoint momentum stresses
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Non-conservative part of the terms added to the mean flow equations
        virtual tmp<volVectorField> nonConservativeMomentumSource() const;

        //- Source term added to the adjoint mean flow due to the
        //- differentiation of the turbulence model
        virtual tmp<volVectorField> adjointMeanFlowSource();

        //- Jacobian of nut wrt to k
        //  Needs to be implemented for objectives related to nut, defined  in
        //  the internal field
        virtual tmp<volScalarField> nutJacobianTMVar1() const;

        //- Jacobian of nut wrt to omega
        //  Needs to be implemented for objectives related to nut, defined  in
        //  the internal field
        virtual tmp<volScalarField> nutJacobianTMVar2() const;

        //- Jacobian of nut wrt the flow velocity
        //  Assumes we want to get contributions of mult*d(nut)/dU
        //  Since the dependency of nut to U is usually through a differential
        //  operator, the term multiplying d(nut)/dU is passed as an argument
        //  to this function; the latter should then compute the
        //  contribution of the afforementioned term to the adjoint mean flow
        //  equations
        virtual tmp<volVectorField> nutJacobianU
        (
            tmp<volScalarField>& dNutdUMult
        ) const;

        //- Diffusion coeff at the boundary for k
        virtual tmp<scalarField> diffusionCoeffVar1(label patchI) const;

        //- Diffusion coeff at the boundary for omega
        virtual tmp<scalarField> diffusionCoeffVar2(label patchI) const;

        virtual const boundaryVectorField& adjointMomentumBCSource() const;

        //- Sensitivity derivative contributions when using the (E)SI approach
        virtual const boundaryVectorField& wallShapeSensitivities();

        virtual const boundaryVectorField& wallFloCoSensitivities();

        //- Contributions to the adjoint eikonal equation (zero for now)
        virtual tmp<volScalarField> distanceSensitivities();

        //- Sensitivity derivative contributions when using the FI approach
        virtual tmp<volTensorField> FISensitivityTerm();

        virtual tmp<scalarField> topologySensitivities
        (
            const word& designVarsName
        ) const;

        //- Nullify all adjoint turbulence model fields and their old times
        virtual void nullify();

        //- Solve the adjoint turbulence equations
        virtual void correct();

        //- Read adjointRASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adjointRASModels
} // End namespace incompressibleAdjoint
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adjointkOmegaSSTTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
