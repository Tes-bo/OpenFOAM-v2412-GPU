/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::atmBuoyancyTurbSource

Group
    grpFvOptionsSources

Description
    Applies sources on \c k and either \c epsilon or \c omega to incorporate
    effects of buoyancy for atmospheric boundary layer modelling.

    Corrections applied to:
    \verbatim
      k         | Turbulent kinetic energy                  [m2/s2]
    \endverbatim

    Corrections applied to either of the below, if exist:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    Required fields:
    \verbatim
      k             | Turbulent kinetic energy                   [m2/s2]
      epsilon/omega | Dissipation rate OR Spec. dissipation rate [m2/s3]/[1/s]
      T             | Temperature                                [K]
      alphat        | Kinematic turbulent thermal conductivity   [m2/s]
    \endverbatim

    References:
    \verbatim
        Buoyancy effects (tags:SKL, ARAL):
            Sogachev, A., Kelly, M., & Leclerc, M. Y. (2012).
            Consistent two-equation closure modelling for atmospheric
            research: buoyancy and vegetation implementations.
            Boundary-layer meteorology, 145(2), 307-327.
            DOI:10.1007/s10546-012-9726-5

            Alletto, M., Radi, A., Adib, J., Langner, J.,
            Peralta, C., Altmikus, A., & Letzel, M. (2018).
            E-Wind: Steady state CFD approach for stratified flows
            used for site assessment at Enercon.
            In Journal of Physics: Conference Series, 1037(7).
            DOI:10.1088/1742-6596/1037/7/072020

        Mixing-length scale limiter for omega (tag:L):
            Langner, J. (2016).
            Implementierung und validierung von RANS-modellen der
            thermisch geschichteten, atmosphärischen grenzschicht.
            Masterarbeit zum thema, Technische Universität Berlin.

        Mixing-length scale estimation (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmBuoyancyTurbSource1
    {
        // Mandatory entries (unmodifiable)
        type             atmBuoyancyTurbSource;

        atmBuoyancyTurbSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Optional (unmodifiable)
            rho          rho;
            Lmax         41.575;
            beta         3.3e-03;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Req'd  | Dflt
      type      | Type name: atmBuoyancyTurbSource  | word   | yes    | -
      rho       | Name of density field             | word   | no     | rho
      Lmax      | Maximum mixing-length scale       | scalar | no     | 41.575
      beta      | Thermal expansion coefficient     | scalar | no     | 3.3e-03
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmBuoyancyTurbSource.C
    atmBuoyancyTurbSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmBuoyancyTurbSource_H
#define fv_atmBuoyancyTurbSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class atmBuoyancyTurbSource Declaration
\*---------------------------------------------------------------------------*/

class atmBuoyancyTurbSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Internal flag to determine the working field is epsilon or omega
        Switch isEpsilon_;

        //- Name of density field
        const word rhoName_;

        //- Maximum mixing-length scale [m]
        const dimensionedScalar Lmax_;

        //- Thermal expansion coefficient [1/K]
        const dimensionedScalar beta_;

        //- Required turbulence model coefficients (copied from turb model)
        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;

        //- Gravitational acceleration vector [m/s2]
        const dimensionedVector g_;


    // Fields

        //- Buoyancy production term [m2/s3]
        volScalarField::Internal B_;


     // Private Member Functions

        //- Compute the buoyancy production term B_
        void calcB();

        //- Return the buoyancy-effect modifier for epsilon-based models
        tmp<volScalarField::Internal> calcC3
        (
            const volScalarField::Internal& k,
            const volScalarField::Internal& epsilon,
            const volScalarField::Internal& G
        ) const;

        //- Return the buoyancy-effect modifier for omega-based models
        tmp<volScalarField::Internal> calcC3
        (
            const volScalarField::Internal& k,
            const volScalarField::Internal& omega,
            const volScalarField::Internal& G,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;

        //- Apply atmBuoyancyTurbSource to epsilon field
        template<class AlphaFieldType, class RhoFieldType>
        void atmBuoyancyTurbSourceEpsilon
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmBuoyancyTurbSource to omega field
        template<class AlphaFieldType, class RhoFieldType>
        void atmBuoyancyTurbSourceOmega
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmBuoyancyTurbSource to k field
        template<class AlphaFieldType, class RhoFieldType>
        void atmBuoyancyTurbSourceK
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;


public:

    //- Runtime type information
    TypeName("atmBuoyancyTurbSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmBuoyancyTurbSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmBuoyancyTurbSource(const atmBuoyancyTurbSource&) = delete;

        //- No copy assignment
        void operator=(const atmBuoyancyTurbSource&) = delete;


    // Member Functions

        //- Add explicit contribution to epsilon or omega equation
        //- for incompressible flow computations
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read source dictionary (effectively no-op)
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "atmBuoyancyTurbSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
