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
    Foam::multiphaseInterSystem

Description

SourceFiles
    multiphaseInterSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseInterSystem_H
#define multiphaseInterSystem_H

#include "basicThermo.H"

#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "compressibleTransportModel.H"
#include "localMin.H"

#include "compressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
namespace multiphaseInter
{
class surfaceTensionModel;
}
class porousModel;

/*---------------------------------------------------------------------------*\
                    Class multiphaseInterSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseInterSystem
:
    public basicThermo,
    public compressibleTransportModel
{
public:

    // Public typedefs

        typedef
            HashTable
            <
                autoPtr<phasePair>, phasePairKey, phasePairKey::hash
            >
            phasePairTable;


        typedef HashTable<autoPtr<multiphaseInter::phaseModel>> phaseModelTable;


        typedef HashTable<volScalarField::Internal> SuSpTable;


protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash> dictTable;


        typedef
            HashTable
            <
                autoPtr<multiphaseInter::surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;


        typedef
            HashTable
            <
                autoPtr<porousModel>,
                phasePairKey,
                phasePairKey::hash
            >
            interfacePorousModelTable;



    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dynamic viscocity
        volScalarField mu_;

        //- Phase names
        wordList phaseNames_;

        //- Mixture total volumetric flux
        surfaceScalarField phi_;

        //- Mixture total mass flux
        surfaceScalarField rhoPhi_;

        //- Phase models
        phaseModelTable phaseModels_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Total ordered phase pairs in the system
        phasePairTable totalPhasePairs_;

        //- Turbulent Prandt number
        dimensionedScalar Prt_;

        //- Turbulence model
        compressibleTurbulenceModel* turb_;


        // Sub Models

            //- Surface tension models
            surfaceTensionModelTable surfaceTensionModels_;

            //- Interface porous models
            interfacePorousModelTable interfacePorousModelTable_;


    // Protected member functions


        //- Calculate and return the laminar viscosity
        void calcMu();

        //- Generate the phases
        HashTable<autoPtr<multiphaseInter::phaseModel>> generatePhaseModels
        (
            const wordList& names
        ) const;

        //- Generate the mixture flux
        tmp<surfaceScalarField> generatePhi
        (
            const HashTable<autoPtr<multiphaseInter::phaseModel>>& phaseModels
        ) const;

        //- Generate pairs
        void generatePairs(const dictTable& modelDicts);

        //- Generate pair table
        void generatePairsTable();

        //- Generate pairs and sub-model tables using pair keys
        template<class modelType>
        void createSubModels
        (
            const dictTable& modelDicts,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and sub-model tables using mesh
        template<class modelType>
        void createSubModels
        (
            const dictTable& modelDicts,
            const fvMesh& mesh,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );


        //- Generate pairs and per-phase sub-model tables with mesh ref
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            const fvMesh& mesh,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );


        //- Generate pairs and per-phase sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                HashTable<autoPtr<modelType>>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );


public:

    //- Runtime type information
    TypeName("multiphaseInterSystem");

    //- Default name of the phase properties dictionary
    static const word phasePropertiesName;


    // Constructors

        //- Construct from fvMesh
        multiphaseInterSystem(const fvMesh& mesh);


    //- Destructor
    virtual ~multiphaseInterSystem();


    // Member Functions

    // Energy related thermo functionaliy functions

        //- Return access to the internal energy field [J/Kg]
        //  \note this mixture thermo is prepared to work with T
        virtual volScalarField& he()
        {
            NotImplemented;
            return const_cast<volScalarField&>(volScalarField::null());
        }

        //- Return access to the internal energy field [J/Kg]
        //  \note this mixture thermo is prepared to work with T
        virtual const volScalarField& he() const
        {
           NotImplemented;
           return volScalarField::null();
        }

        //- Enthalpy/Internal energy
        //  for given pressure and temperature [J/kg]
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

        //- Chemical enthalpy of the mixture [J/kg]
        virtual tmp<volScalarField> hc() const;

        //- Temperature from enthalpy/internal energy for cell-set
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& p,
            const scalarField& T0,
            const labelList& cells
        ) const;

        //- Temperature from enthalpy/internal energy for patch
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& p,
            const scalarField& T0,
            const label patchi
        ) const;


    // Thermo

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const;

        //- Return the mixture density on a patch
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
        ) const;

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

        //- Molecular weight [kg/kmol] of the mixture
        virtual tmp<volScalarField> W() const;


    // Transport

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        virtual tmp<scalarField> kappa
        (
            const label patchi
        ) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Effective thermal diffusivity for temperature
        //  of mixture [J/m/s/K]
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& kappat
        ) const;

        //- Effective thermal diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
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

        //- Return Prandt number
        const dimensionedScalar& Prt() const;


    // Access to transport state variables

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;


    // Turbulence

        //- Set turbulence model
        void setTurbulence(compressibleTurbulenceModel& turb)
        {
            turb_ = &turb;
        }

        //- Return pointer to turbulence model
        compressibleTurbulenceModel* turbulence() const
        {
            return turb_;
        }

        //- Return the turbulent dynamic viscosity
        tmp<volScalarField> mut() const;

        //- Return the effective dynamic viscosity
        tmp<volScalarField> muEff() const;

        //- Return the turbulent kinematic viscosity
        tmp<volScalarField> nut() const;

        //- Return the effective kinematic viscosity
        tmp<volScalarField> nuEff() const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        tmp<volScalarField> kappaEff() const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        tmp<scalarField> kappaEff(const label patchi) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        tmp<volScalarField> alphaEff() const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        tmp<scalarField> alphaEff(const label patchi) const;


    // Phase fluxes

        //- Constant access to the total flux
        const surfaceScalarField& phi() const;

        //- Access to the total mixture flux
        surfaceScalarField& phi();

        //- Constant access to the mixture mass flux
        const surfaceScalarField& rhoPhi() const;

        //- Access to the total mixture mass flux
        surfaceScalarField& rhoPhi();

        //- Mixture U
        tmp<volVectorField> U() const;


    // Surface tension

        //- Calculate surface tension of the mixture
        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Return the surface tension coefficient
        virtual tmp<volScalarField> surfaceTensionCoeff
        (
            const phasePairKey& key
        ) const;


        //- Return coefficients (1/rho)
        virtual tmp<volScalarField> coeffs(const word& key) const;


    // Interface porous between solid/fluid phases

        //- Add interface porosity on phasePair
        void addInterfacePorosity(fvVectorMatrix& UEqn);


    // Inter-Phase mass and heat transfe

        //- Return interfacial source mass rate per phase pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const = 0;

        //- Return the heat transfer matrices
        virtual tmp<fvScalarMatrix> heatTransfer
        (
            const volScalarField& T
        ) = 0;

        //- Return the volumetric rate transfer matrix
        virtual tmp<fvScalarMatrix> volTransfer
        (
            const volScalarField& p
        ) = 0;

        //- Calculate mass transfer for alpha's
        virtual void alphaTransfer(SuSpTable& Su, SuSpTable& Sp) = 0;

        //- Calculate mass transfer for species
        virtual void massSpeciesTransfer
        (
            const multiphaseInter::phaseModel& phase,
            volScalarField::Internal& Su,
            volScalarField::Internal& Sp,
            const word speciesName
        ) = 0;

        //- Add volume change in pEq
        virtual bool includeVolChange() = 0;


    // Solve phases and correct models

        //- Solve for the phase transport equations
        virtual void solve() = 0;

        //- Correct the mixture thermos
        virtual void correct();

        //- Correct mass sources
        virtual void correctMassSources(const volScalarField& T) = 0;

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            NotImplemented;
            return word();
        }

        //- Correct the turbulence
        //  \note Each phase could have its own turbulence
        virtual void correctTurbulence();

        //- Read base phaseProperties dictionary
        virtual bool read();


    // Access to phases models

        //- Constant access the total phase pairs
        const phasePairTable& totalPhasePairs() const;

        //- Non-constant access the total phase pairs
        phasePairTable& totalPhasePairs();

        //- Constant access the phases
        const phaseModelTable& phases() const;

        //- Access the phases
        phaseModelTable& phases();

        //- Access a sub model between a phase pair
        template <class modelType>
        const modelType& lookupSubModel(const phasePair& key) const;

        //- Access a sub model between two phases
        template <class modelType>
        const modelType& lookupSubModel
        (
            const multiphaseInter::phaseModel& from,
            const multiphaseInter::phaseModel& to
        ) const;


    // Query phase thermo information

        //- Return true if the equation of state is incompressible for all
        //  phases
        virtual bool incompressible() const;

        //- Return true if a phase is incompressible
        virtual bool incompressible(const word) const;

        //- Return true if the equation of state is isochoric for all phasses
        //  i.e. rho = const
        virtual bool isochoric() const;

        //- Return mesh
        const fvMesh& mesh() const;


    // Help functions for the interfaces

        //- Interface normal surface vector
        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Interface normal volField vector
        tmp<volVectorField> nVolHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Interface normal surface vector
        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Interface curvature
        tmp<volScalarField> K
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Near Interface of alpha1 and alpha2
        tmp<volScalarField> nearInterface
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        //- Near Interface of alpha'n
        tmp<volScalarField> nearInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiphaseInterSystemTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
