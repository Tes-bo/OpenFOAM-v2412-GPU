/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::incompressibleAdjoint::adjointRASModels::adjointSpalartAllmaras

Description
    Continuous adjoint to the Spalart-Allmaras one-eqn mixing-length model for
    incompressible flows.

    Reference:
    \verbatim
        For the adjoint to the Spalart-Allmaras PDE

            Zymaris, A., Papadimitriou, D., Giannakoglou, K., &
            Othmer, C. (2009).
            Continuous adjoint approach to the Spalart-Allmaras turbulence
            model for incompressible flows.
            Computers & Fluids, 38(8), 1528-538.
            http://doi.org/10.1016/j.compfluid.2008.12.006

        For the FI sensitivity terms

            Papoutsis-Kiachagias, E. M., Asouti, V. G., Giannakoglou, K. C.,
            Gkagkas, K., Shimokawa, S., & Itakura, E. (2018).
            Multi-point aerodynamic shape optimization of cars based on
            continuous adjoint.
            Structural and Multidisciplinary Optimization, 59(2), 675-694.
            http://doi.org/10.1007/s00158-018-2091-3

    \endverbatim

    Both of the above-mentioned papers develop the adjoint to the
    Spalart-Allmaras PDE that includes the fv3 term.  The current
    implementation corresponds to the Spalart-Allmaras PDE as programmed within
    OpenFOAM and is, thus, slightly different than the one developed in the
    cited papers

SourceFiles
    adjointSpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef adjointSpalartAllmaras_H
#define adjointSpalartAllmaras_H

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
                    Class adjointSpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class adjointSpalartAllmaras
:
    public adjointRASModel
{
    // Private Member Functions

        //- No copy construct
        adjointSpalartAllmaras(const adjointSpalartAllmaras&) = delete;

        //- No copy assignment
        void operator=(const adjointSpalartAllmaras&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;

            //- Whether to limit the adjoint production term to enhance
            //- stability
            bool limitAdjointProduction_;


        // Fields

            //- Wall distance
            //  Note: reference to the distance known by the primal model
            const volScalarField& y_;

            //- Field for masking (convergence enhancement)
            volScalarField mask_;

            // Fields that depend only on primal fields and are very expensive
            // to compute in each iteration.
            // Store and update them when the primal solution has been updated
            volSymmTensorField symmAdjointProductionU_;
            volScalarField productionDestructionSource_;
            volScalarField Stilda_;
            volScalarField r_;
            volScalarField fw_;
            volScalarField Cdnut_;
            volTensorField momentumSourceMult_;
            volTensorField gradU_;
            volVectorField gradNuTilda_;

            // Useful quantities for bounding
            dimensionedScalar minStilda_;


    // Protected Member Functions

        // Primal Spalart - Allmaras

            tmp<volScalarField> chi() const;

            tmp<volScalarField> fv1(const volScalarField& chi) const;

            tmp<volScalarField> fv2
            (
                const volScalarField& chi,
                const volScalarField& fv1
            ) const;

            tmp<volScalarField> Stilda
            (
                const volScalarField& chi,
                const volScalarField& fv1
            ) const;

            tmp<volScalarField> r(const volScalarField& Stilda) const;

            tmp<volScalarField> fw(const volScalarField& Stilda) const;

            tmp<volScalarField> DnuTildaEff() const;

            //- References to the primal turbulence model variables
            const volScalarField& nuTilda() const;

            const volScalarField& nutRef() const;


        // Adjoint Spalart - Allmaras

            // Differentiation of the primal Spalart - Allmaras terms

                tmp<volScalarField> dFv1_dChi
                (
                    const volScalarField& chi
                ) const;

                tmp<volScalarField> dFv2_dChi
                (
                    const volScalarField& chi,
                    const volScalarField& fv1,
                    const volScalarField& dFv1dChi
                ) const;

                tmp<volScalarField> dStilda_dOmega
                (
                    const volScalarField& Omega,
                    const volScalarField& fv2
                ) const;

                tmp<volScalarField> dStilda_dNuTilda
                (
                    const volScalarField& Omega,
                    const volScalarField& fv2,
                    const volScalarField& dFv2dChi
                ) const;

                tmp<volScalarField> dStilda_dDelta
                (
                    const volScalarField& Omega,
                    const volScalarField& fv2
                ) const;

                tmp<volScalarField> dr_dNuTilda
                (
                    const volScalarField& Stilda
                ) const;

                tmp<volScalarField> dr_dStilda
                (
                    const volScalarField& Stilda
                ) const;

                tmp<volScalarField> dr_dDelta
                (
                    const volScalarField& Stilda
                ) const;

                tmp<volScalarField> dfw_dr
                (
                    const volScalarField& Stilda
                ) const;

                tmp<volScalarField> dfw_dNuTilda
                (
                    const volScalarField& Stilda,
                    const volScalarField& dfwdr,
                    const volScalarField& dStildadNuTilda
                ) const;

                tmp<volScalarField> dfw_dOmega
                (
                    const volScalarField& Stilda,
                    const volScalarField& dfwdr,
                    const volScalarField& dStildadOmega
                ) const;

                tmp<volScalarField> dfw_dDelta
                (
                    const volScalarField& Stilda,
                    const volScalarField& dfwdr,
                    const volScalarField& dStildadDelta
                ) const;

                tmp<volScalarField> dD_dNuTilda
                (
                    const volScalarField& fw,
                    const volScalarField& dfwdNuTilda
                ) const;

                tmp<volScalarField> dP_dNuTilda
                (
                    const volScalarField& dStildadNuTilda
                ) const;

                tmp<volScalarField> dnut_dNuTilda
                (
                    const volScalarField& fv1,
                    const volScalarField& dFv1dChi
                ) const;


            //- Conservative source term for the adjoint momentum equations
            //  Sets also the adjointMomentumBCSource
            tmp<volVectorField> conservativeMomentumSource();

            //- Access to the adjoint Spalart - Allmaras field
            inline volScalarField& nuaTilda()
            {
                return adjointTMVariable1Ptr_();
            };

            //- Constant access to the adjoint Spalart - Allmaras field
            inline const volScalarField& nuaTilda() const
            {
                return adjointTMVariable1Ptr_();
            };

            //- Update the constant primal-related fields
            void updatePrimalRelatedFields();

            //- Allocate the mask field
            tmp<volScalarField> allocateMask();


public:

    //- Runtime type information
    TypeName("adjointSpalartAllmaras");


    // Constructors

        //- Construct from components
        adjointSpalartAllmaras
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName
                = adjointTurbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~adjointSpalartAllmaras() = default;


    // Member Functions

        virtual tmp<volSymmTensorField> devReff() const;

        virtual tmp<volSymmTensorField> devReff(const volVectorField& U) const;

        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        virtual tmp<volVectorField> adjointMeanFlowSource();

        virtual tmp<volScalarField> nutJacobianTMVar1() const;

        virtual tmp<scalarField> diffusionCoeffVar1(label patchI) const;

        virtual const boundaryVectorField& adjointMomentumBCSource() const;

        virtual const boundaryVectorField& wallShapeSensitivities();

        virtual const boundaryVectorField& wallFloCoSensitivities();

        virtual tmp<volScalarField> distanceSensitivities();

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

#endif

// ************************************************************************* //
