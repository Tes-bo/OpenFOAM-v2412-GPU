/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::RASModels::EBRSM

Group
    grpRASTurbulence

Description
    Manceau (2015)'s elliptic-blending Reynolds-stress turbulence model
    for incompressible and compressible flows.

    References:
    \verbatim
        Standard model (tag:M):
            Manceau, R. (2015).
            Recent progress in the development of the elliptic
            blending Reynolds-stress model.
            International Journal of Heat and Fluid Flow, 51, 195-220.
            DOI:10.1016/j.ijheatfluidflow.2014.09.002

        Simple gradient diffusion hypothesis (tag:LM):
            Lardeau, S., & Manceau, R. (2014).
            Computations of complex flow configurations using
            a modified elliptic-blending Reynolds-stress model.
            10th International ERCOFTAC Symposium on Engineering
            Turbulence Modelling and Measurements. Marbella, Spain.
            https://hal.archives-ouvertes.fr/hal-01051799
    \endverbatim

    The default model coefficients are (M:p. 219):
    \verbatim
        EBRSMCoeffs
        {
            g1          3.4;
            g1star      1.8;
            g3          0.8;
            g3star      1.3;
            g4          1.25;
            g5          0.2;
            Cmu         0.21;
            Ceps1       1.44;
            Ceps2       1.83;
            sigmaK      1.0;
            sigmaEps    1.15;
            A1          0.065;
            Ct          6.0;
            Cl          0.133;
            Ceta        80.0;
            Cstability  10.0;

            simpleGradientDiffusion  false;
        }
    \endverbatim

Note
  - \c g5 coefficient has been changed from its original value of 0.4 to 0.2
    after obtaining better results in smooth-wall plane channel flow cases.

SourceFiles
    EBRSM.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_RASModels_EBRSM_H
#define Foam_RASModels_EBRSM_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class EBRSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class EBRSM
:
    public ReynoldsStress<RASModel<BasicTurbulenceModel>>
{
    // Private Data

            Switch simpleGradientDiffusion_;

        // Coefficients

            dimensionedScalar g1_;
            dimensionedScalar g1star_;
            dimensionedScalar g3_;
            dimensionedScalar g3star_;
            dimensionedScalar g4_;
            dimensionedScalar g5_;
            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar A1_;
            dimensionedScalar Ct_;
            dimensionedScalar Cl_;
            dimensionedScalar Ceta_;
            dimensionedScalar Cstability_;

        // Fields

            //- Elliptic blending factor [-]
            volScalarField f_;

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;


    // Private Member Functions

        //- Return the Durbin-limited length-scale field
        //  Thickness of region of influence of near-wall model
        tmp<volScalarField> calcL() const;

        //- Return the approximate wall-normal vector field
        tmp<volVectorField> calcN() const;

        //- Return the turbulent time scale bounded by Kolmogorov time scale
        tmp<volScalarField> calcTau() const;

        //- Return the effective diffusivity for epsilon or R
        //- by using the generalised gradient diffusion hypothesis
        tmp<volSymmTensorField> D
        (
            const volScalarField& tau,
            const dimensionedScalar& sigma
        ) const;

        //- Return the effective diffusivity for epsilon or R
        //- by using the simple gradient diffusion hypothesis
        tmp<volScalarField> D
        (
            const dimensionedScalar& sigma
        ) const;

        //- Return dissipation-production stimulator in the buffer layer
        tmp<volScalarField> Ceps1Prime(const volScalarField& G) const;

        //- Update the eddy-viscosity
        virtual void correctNut();


        //- No copy construct
        EBRSM(const EBRSM&) = delete;

        //- No copy assignment
        void operator=(const EBRSM&) = delete;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("EBRSM");


    // Constructors

        //- Construct from components
        EBRSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~EBRSM() = default;


    // Member Functions

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EBRSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
