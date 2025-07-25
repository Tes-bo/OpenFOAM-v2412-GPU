/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::twoPhaseMixtureEThermo

Description

SourceFiles
    twoPhaseMixtureEThermo.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixtureEThermo_H
#define twoPhaseMixtureEThermo_H

#include "volFields.H"

#include "basicThermo.H"
#include "thermoIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class twoPhaseMixtureEThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixtureEThermo
:
    public basicThermo,
    public thermoIncompressibleTwoPhaseMixture
{

protected:

    // Protected Data

        //- Saturation Temperature
        dimensionedScalar TSat_;

public:

    TypeName("twoPhaseMixtureEThermo");

    // Constructor
    twoPhaseMixtureEThermo
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );


    //- Destructor
    virtual ~twoPhaseMixtureEThermo() = default;


    // Member Functions


        //- Return access to the internal energy field [J/Kg]
        virtual volScalarField& he()
        {
            NotImplemented;
            return const_cast<volScalarField&>(volScalarField::null());
        }

        //- Return access to the internal energy field [J/Kg]
        virtual const volScalarField& he() const
        {
            NotImplemented;
            return volScalarField::null();
        }

        //- Enthalpy/Internal energy
        //- for given pressure and temperature [J/kg]
        virtual tmp<volScalarField> he
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Enthalpy/Internal energy for cell-set [J/kg]
        virtual tmp<scalarField> he
        (
            const scalarField& p,
            const scalarField& T,
            const labelList& cells
        ) const;

        //- Enthalpy/Internal energy for patch [J/kg]
        virtual tmp<scalarField> he
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Chemical enthalpy [J/kg]
        virtual tmp<volScalarField> hc() const;

        //- Temperature from enthalpy/internal energy for cell-set
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& p,
            const scalarField& T0,      // starting temperature
            const labelList& cells
        ) const;

        //- Temperature from enthalpy/internal energy for patch
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& p,
            const scalarField& T0,      // starting temperature
            const label patchi
        ) const;


        //- Return true if the equation of state is incompressible
        //- i.e. rho != f(p)
        bool incompressible() const
        {
            return true;
        }

        //- Return true if the equation of state is isochoric
        //- i.e. rho = const
        bool isochoric() const
        {
            return false;
        }

        //- Return rho of the mixture
        virtual tmp<volScalarField> rho() const;

        //- Return rho for patch
        virtual tmp<scalarField> rho(const label patchi) const;

        //- Return Cp of the mixture
        virtual tmp<volScalarField> Cp() const;

        //- Heat capacity at constant pressure for patch [J/kg/K]
        virtual tmp<scalarField> Cp
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity using pressure and temperature
        virtual tmp<scalarField> Cp
        (
            const scalarField& p,
            const scalarField& T,
            const labelList& cells
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Return Cv of the mixture
        virtual tmp<volScalarField> Cv() const;

        //- Heat capacity at constant volume for patch [J/kg/K]
        virtual tmp<scalarField> Cv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchI
        ) const;

        //- Density from pressure and temperature
        virtual tmp<scalarField> rhoEoS
        (
            const scalarField& p,
            const scalarField& T,
            const labelList& cells
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Gamma = Cp/Cv []
        virtual tmp<volScalarField> gamma() const;

        //- Gamma = Cp/Cv for patch []
        virtual tmp<scalarField> gamma
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

         //- Heat capacity at constant pressure/volume [J/kg/K]
        virtual tmp<volScalarField> Cpv() const;

        //- Heat capacity at constant pressure/volume for patch [J/kg/K]
        virtual tmp<scalarField> Cpv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity ratio []
        virtual tmp<volScalarField> CpByCpv() const;

        //- Heat capacity ratio for patch []
        virtual tmp<scalarField> CpByCpv
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Molecular weight [kg/kmol]
        virtual tmp<volScalarField> W() const;

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature
        //- of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappa
        (
            const label patchi
        ) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Effective thermal diffusivity for temperature
        //- of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField&
        ) const;

        //- Effective thermal diffusivity for temperature
        //- of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Effective thermal diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const;

        //- Effective thermal diffusivity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;


        //- Correct the thermo fields
        virtual void correct();

        //- Return the name of the thermo physics
        virtual word thermoName() const;

        //- Read properties
        virtual bool read();


        // Access to thermodynamic state variables

            //- Return const-access to the saturation temperature
            const dimensionedScalar& TSat() const
            {
                return TSat_;
            }

            //- Return transport properties dictionary
            const incompressibleTwoPhaseMixture& transportPropertiesDict()
            {
                return *this;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
