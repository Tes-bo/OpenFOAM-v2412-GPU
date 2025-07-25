/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::combustionModels::diffusionMulticomponent

Group
    grpCombustionModels

Description
    Diffusion based turbulent combustion model for multicomponent species.

    The model calculates the laminar finite rate source terms based on
    the kinetic for each reaction in order to begin the combustion and
    evaluates the minimum between this and the cross diffusion rate term
    defined as D*prob*muEff*mag(grad(Yi)*grad(Yj)) if laminarIgn is true.

    where:

    D :     is a model dynamic constant defined as C*f02 where:
            C is a model constant
            f02 = 1  + sqr(O2/oxidantRes), oxidantRes is an user input

    muEff : is the effective turbulent viscosity
    prob :  is a normalized Gaussian shaped distribution around the stoichiometric
            value of each reaction. The distribution is controlled by 'sigma'
            for standard deviation and ftCorr for correction of the stoichiometric
            value.

    In the combustion properties dictionary:

    diffusionMulticomponentCoeffs
    {
        Ci          (1.0 1.0);      // Default to 1
        fuels       (CH4 CO);
        oxidants    (O2 O2);
        YoxStream   (0.23 0.23);    // Default to 0.23
        YfStream    (1.0 1.0);      // Default to 1.0
        sigma       (0.02 0.02);    // Default to 0.02
        oxidantRes  (0.025 0.005);
        ftCorr      (0.0 0.0);      // Default to 0.0
        laminarIgn  false;          // Default false
    }


SourceFiles
    diffusionMulticomponent.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionMulticomponent_H
#define diffusionMulticomponent_H

#include "ChemistryCombustion.H"
#include "scalarList.H"
#include "tmp.H"
#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    Class diffusionMulticomponent Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class diffusionMulticomponent
:
    public ChemistryCombustion<ReactionThermo>
{
    // Private data

        //- Reactions
        const PtrList<Reaction<ThermoType>>& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Pointer list of source terms
        PtrList<volScalarField> RijPtr_;

        //- Model constants
        scalarList Ci_;

        //- List of fuels for each reaction
        wordList fuelNames_;

        //- List of oxidants for each reaction
        wordList oxidantNames_;

        //- Heat of combustion [J/Kg]
        scalarList qFuel_;

        //- Stoichiometric air-fuel mass ratio
        scalarList stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        scalarList s_;

        //- Oxidizer stream mass concentrations
        scalarList YoxStream_;

        //- Fuel stream mass concentrations
        scalarList YfStream_;

        //- Mean distribution for gaussian probability
        scalarList sigma_;

        //- Residual oxidizer
        scalarList oxidantRes_;

        //- ft stoichiometric correction
        scalarList ftCorr_;

        //- Relaxation factor on total source
        scalar alpha_;

        //- Switch on to laminar combustion for ignition
        bool laminarIgn_;


    // Private Member Functions

        //- Return the chemical time scale
        tmp<volScalarField> tc() const;

        //- Initialize
        void init();

        //- No copy construct
        diffusionMulticomponent(const diffusionMulticomponent&) = delete;

        //- No copy assignment
        void operator=(const diffusionMulticomponent&) = delete;


public:

    //- Runtime type information
    TypeName("diffusionMulticomponent");


    // Constructors

        //- Construct from components
        diffusionMulticomponent
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~diffusionMulticomponent() = default;


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix.
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> Qdot() const;


    // IO

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "diffusionMulticomponent.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
