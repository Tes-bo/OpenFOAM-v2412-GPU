/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::RASModels::kEpsilonPhitF

Group
    grpRASTurbulence

Description
    The k-epsilon-phit-f turbulence closure model for incompressible and
    compressible flows.

    The model is a three-transport-equation linear-eddy-viscosity turbulence
    closure model alongside an elliptic relaxation equation.

    \heading Input fields
    \plaintable
        k         | Turbulent kinetic energy [m2/s2]
        epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
        phit      | Normalised wall-normal fluctuating velocity scale [-]
        f         | Elliptic relaxation factor [1/s]
    \endplaintable

    Reference:
    \verbatim
        Standard model (Tag:LUU):
            Laurence, D. R., Uribe, J. C., & Utyuzhnikov, S. V. (2005).
            A robust formulation of the v2-f model.
            Flow, Turbulence and Combustion, 73(3-4), 169–185.
            DOI:10.1007/s10494-005-1974-8
    \endverbatim

    The default model coefficients are (LUU:Eqs. 19-20):
    \verbatim
        kEpsilonPhitFCoeffs
        {
            includeNu   true;    // include nu in (LUU: Eq. 17), see Notes
            Cmu         0.22;    // Turbulent viscosity constant
            Ceps1a      1.4;     // Model constant for epsilon
            Ceps1b      1.0;     // Model constant for epsilon
            Ceps1c      0.05;    // Model constant for epsilon
            Ceps2       1.9;     // Model constant for epsilon
            Cf1         1.4;     // Model constant for f
            Cf2         0.3;     // Model constant for f
            CL          0.25;    // Model constant for L
            Ceta        110.0;   // Model constant for L
            CT          6.0;     // Model constant for T
            sigmaK      1.0;     // Turbulent Prandtl number for k
            sigmaEps    1.3;     // Turbulent Prandtl number for epsilon
            sigmaPhit   1.0;     // Turbulent Prandtl number for phit = sigmaK
        }
    \endverbatim

Note
    The name of the original variable replacing 'v2' is 'phi' (LUU:Eq. 14).
    However, the name 'phi' preexisted in OpenFOAM; therefore, this name was
    replaced by 'phit' herein.

    Including \c nu in \c DphitEff even though it is not present in (LUU:Eq. 17)
    provided higher level of resemblance to benchmarks for the tests considered,
    particularly for the peak skin friction (yet, pressure-related predictions
    were unaffected). Users can switch off \c nu in \c DphitEff by using
    \c includeNu entry in \c kEpsilonPhitFCoeffs as shown above in order to
    follow the reference paper thereat. \c includeNu is left \c true by default.
    See GitLab issue #1560.

SourceFiles
    kEpsilonPhitF.C

SeeAlso
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonPhitF_H
#define kEpsilonPhitF_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kEpsilonPhitF Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonPhitF
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kEpsilonPhitF(const kEpsilonPhitF&) = delete;

        //- No copy assignment
        void operator=(const kEpsilonPhitF&) = delete;


protected:

    // Protected Data

            Switch includeNu_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1a_;
            dimensionedScalar Ceps1b_;
            dimensionedScalar Ceps1c_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cf1_;
            dimensionedScalar Cf2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar CT_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaPhit_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //- Normalised wall-normal fluctuating velocity scale [-]
            volScalarField phit_;

            //- Elliptic relaxation factor [1/s]
            volScalarField f_;

            //- Turbulent time scale [s]
            volScalarField T_;


        // Bounding values

            dimensionedScalar phitMin_;
            dimensionedScalar fMin_;
            dimensionedScalar TMin_;
            dimensionedScalar L2Min_;


    // Protected Member Functions

        //- Update nut with the latest available k, phit, and T
        virtual void correctNut();

        //- Return the turbulent time scale, T
        tmp<volScalarField> Ts() const;

        //- Return the turbulent length scale, L
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonPhitF");


    // Constructors

        //- Construct from components
        kEpsilonPhitF
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
    virtual ~kEpsilonPhitF() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k (LUU:Eq. 3)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for epsilon (LUU:Eq. 4)
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for phit (LUU:Eq. 17)
        tmp<volScalarField> DphitEff() const
        {
            auto tfld =
                tmp<volScalarField>::New("DphitEff", this->nut_/sigmaPhit_);

            if (includeNu_)
            {
                tfld.ref() += this->nu();
            }

            return tfld;
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> phit() const
        {
            return phit_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonPhitF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
