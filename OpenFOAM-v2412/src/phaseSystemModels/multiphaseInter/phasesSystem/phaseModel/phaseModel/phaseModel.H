/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::phaseModel

Description

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_multiphaseInter_phaseModel_H
#define Foam_multiphaseInter_phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class multiphaseInterSystem;

namespace multiphaseInter
{
/*---------------------------------------------------------------------------*\
                         Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private Data

        //- Reference to the multiphaseInterSystem to which this phase belongs
        const multiphaseInterSystem& fluid_;

        //- Name of phase
        word name_;


public:

    //- Runtime type information
    ClassName("phaseModel");

    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseModel,
        multiphaseInterSystem,
        (
            const multiphaseInterSystem& fluid,
            const word& phaseName
        ),
        (fluid, phaseName)
    );


    // Constructors

        //- Construct from multiphaseInterSystem and phaseName
        phaseModel(const multiphaseInterSystem& fluid, const word& phaseName);


    //- Destructor
    virtual ~phaseModel() = default;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const multiphaseInterSystem& fluid,
            const word& phaseName
        );


    // Member Functions

        //- The name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the system to which this phase belongs
        const multiphaseInterSystem& fluid() const;

        //- Correct phase thermo
        virtual void correct();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Solve species fraction equation
        virtual void solveYi
        (
            PtrList<volScalarField::Internal>& Su,
            PtrList<volScalarField::Internal>& Sp
        ) = 0;

        //- Read phase properties dictionary
        virtual bool read();


    // Thermo

        //- Access const to phase thermo
        virtual const rhoThermo& thermo() const = 0;

        //- Access to phase thermo
        virtual rhoThermo& thermo() = 0;

        //- Return the phase density
        tmp<volScalarField> rho() const;

        //- Return phase density on a patch
        tmp<scalarField> rho(const label patchi) const;

        //- Chemical enthalpy for phase [J/kg]
        tmp<volScalarField> hc() const;

        //- Return phase Cp
        tmp<volScalarField> Cp() const;

        //- Heat capacity of the phase at constant pressure for patch
        // [J/kg/K]
        tmp<scalarField> Cp
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Return Cv of the phase
        tmp<volScalarField> Cv() const;

        //- Heat capacity at constant volume for phase for a patch [J/kg/K]
        tmp<scalarField> Cv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchI
        ) const;

        //- Gamma = Cp/Cv of phase[]
        tmp<volScalarField> gamma() const;

        //- Gamma = Cp/Cv for phase on patch []
        tmp<scalarField> gamma
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant pressure/volume for phase [J/kg/K]
        tmp<volScalarField> Cpv() const;

        //- Heat capacity at constant pressure/volume for phase at patch
        // [J/kg/K]
        tmp<scalarField> Cpv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity ratio for phase []
        tmp<volScalarField> CpByCpv() const;

        //- Heat capacity ratio for phase at patch []
        tmp<scalarField> CpByCpv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Query thermo for dpdt
        bool dpdt() const
        {
             return thermo().dpdt();
        }


    // Transport

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        const volScalarField& alpha() const;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        const scalarField& alpha(const label patchi) const;

        //- Thermal diffusivity for temperature of phase [J/m/s/K]
        tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature of phase for patch [J/m/s/K]
        tmp<scalarField> kappa(const label patchi) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        tmp<scalarField> alphahe(const label patchi) const;

        //- Effective thermal diffusivity for temperature of phase [J/m/s/K]
        tmp<volScalarField> kappaEff(const volScalarField&) const;

        //- Effective thermal diffusivity for temperature
        //  of phase for patch [J/m/s/K]
        tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Effective thermal diffusivity of phase [kg/m/s]
        tmp<volScalarField> alphaEff(const volScalarField& alphat) const;

        //- Effective thermal diffusivity of phase for patch [kg/m/s]
        tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Return the mixture kinematic viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the mixture kinematic viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return the mixture dymanic viscosity
        virtual tmp<volScalarField> mu() const;

        //- Return the mixture dymanic viscosity on patchi
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Diffusion number
        virtual tmp<surfaceScalarField> diffNo() const = 0;


    // Species

        //- Constant access the species mass fractions
        virtual const PtrList<volScalarField>& Y() const = 0;

        //- Access the species mass fractions
        virtual PtrList<volScalarField>& Y() = 0;


    // Momentum

        //- Constant access the volumetric flux
        virtual tmp<surfaceScalarField> phi() const = 0;

        //- Access the volumetric flux
        virtual const surfaceScalarField& phi() = 0;

        //- Constant access the volumetric flux of the phase
        virtual tmp<surfaceScalarField> alphaPhi() const = 0;

        //- Access the volumetric flux of the phase
        virtual surfaceScalarField& alphaPhi() = 0;

        //- Access const reference to U
        virtual tmp<volVectorField> U() const = 0;


    // Turbulence (WIP: possible to add turbulence on each phase)

        /*
        //- Return the turbulent dynamic viscosity
        virtual tmp<volScalarField> mut() const = 0;

        //- Return the turbulent dynamic viscosity on a patch
        virtual tmp<scalarField> mut(const label patchI) const = 0;

        //- Return the turbulent kinematic viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the turbulent kinematic viscosity on a patch
        virtual tmp<scalarField> nut(const label patchI) const = 0;

        //- Return the kinetic pressure derivative w.r.t. volume fraction
        virtual tmp<volScalarField> pPrime() const = 0;

        //- Return the turbulent kinetic energy
        virtual tmp<volScalarField> k() const = 0;
        */
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multiphaseInter
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
